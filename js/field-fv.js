/* ============================================================================
   field-fv.js — the live background, done compressibly.

   A drop-in replacement for js/field.js. Same public surface, same dust, same
   palette, same per-chapter presets; entirely different gas.

   field.js is incompressible: semi-Lagrangian advection, vorticity confinement,
   and an elliptic pressure solved by multigrid. It is a good-looking scheme and
   it has no shocks in it anywhere, by construction.

   This is the scheme from mini-ramses-ism's GPU hydro path (gpu_hydro.cuf),
   reduced to two dimensions and to its cheapest configuration:

     variables   : conserved (rho, rho u, rho v, rho Y), isothermal closure
                   P = rho cs^2 with cs = 1, so every velocity here is a Mach
                   number. Y is the dye: a passive scalar, which is exactly what
                   RAMSES's scalar_flux does with it in the HLL branch --
                   componentwise HLL on rho*Y with flux u*(rho*Y) -- so it is
                   simply the fourth conserved variable and costs nothing.
     reconstruct : none. Piecewise constant, i.e. slope_type = 0, which in the
                   reference sets dq = zero and with it every source term in the
                   trace step. First-order Godunov.
     flux        : HLL, with RAMSES's wave-speed estimate
                     SL = min(min(uL,uR) - cs, 0)
                     SR = max(max(uL,uR) + cs, 0)
                   Clamping through zero means SL <= 0 <= SR always, so the
                   central formula is always the right one, the supersonic
                   branches vanish, and SR - SL >= 2 cs needs no divide guard.
                   LLF is selectable and cheaper still.
     update      : unsplit and conservative, all four faces in one pass.
     timestep    : cmpdt, verbatim. ctot = sum over ndim of (|u_dim| + cs);
                   dt = C dx / ctot with a ceiling at C dx / smallc. C = 0.8, the
                   RAMSES default, and no margin -- because it is recomputed from
                   the current state before every step, reduced on the GPU into a
                   single texel so there is no readback in the loop.

   The dust is unchanged from field.js: a polydisperse population with stopping
   times drawn log-uniformly, advanced by a first-order implicit (backward Euler)
   drag on the gas velocity, lit by its drift |v - u_gas| and sized in proportion
   to its grain radius. What differs is what it is falling through. Grains in a
   compressible flow pile up against shock fronts, which is a different and more
   honest picture of preferential concentration than a divergence-free field can
   give.

   Single precision throughout, which is the only option: WebGL2's highp float is
   IEEE binary32. GLSL ES has no -ffast-math, so that intent is in the code
   instead -- reciprocal multiplies, no sqrt in the Riemann solver (an isothermal
   sound speed is a uniform), and a branch-free flux.

   Written from scratch. No libraries.
   ========================================================================= */

(function (global) {
  'use strict';

  // ---------------------------------------------------------------- shaders

  const VERT = `#version 300 es
precision highp float;
in vec2 aPos;
out vec2 vUv;
void main(){ vUv = aPos * 0.5 + 0.5; gl_Position = vec4(aPos, 0.0, 1.0); }`;

  const HEAD = `#version 300 es
precision highp float;
precision highp sampler2D;
in vec2 vUv;
out vec4 outColor;
`;

  // Isothermal equation of state. Primitives ride as (rho, u_normal,
  // u_transverse, Y) so one flux expression serves both sweep directions, the
  // same way the reference always puts the normal velocity in velocity_x.
  const EOS = `
uniform float cs, cs2, smallr;

vec4 toPrim(vec4 U){
  float r  = max(U.x, smallr);
  float ri = 1.0 / r;                        // one divide, then multiplies
  return vec4(r, U.y * ri, U.z * ri, U.w * ri);
}
vec4 toCons(vec4 q){
  return vec4(q.x, q.x * q.y, q.x * q.z, q.x * q.w);
}
`;

  // ---- Courant condition, mirroring cmpdt, reduced on the GPU --------------

  const F_CMAX_STATE = HEAD + EOS + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 8;
  float m = 0.0;
  for (int j = 0; j < 8; j++) {
    for (int i = 0; i < 8; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      vec4 q = toPrim(texelFetch(uSrc, p, 0));
      m = max(m, abs(q.y) + abs(q.z) + 2.0 * cs);
    }
  }
  outColor = vec4(m, 0.0, 0.0, 1.0);
}`;

  const F_CMAX_DOWN = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 8;
  float m = 0.0;
  for (int j = 0; j < 8; j++) {
    for (int i = 0; i < 8; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      m = max(m, texelFetch(uSrc, p, 0).x);
    }
  }
  outColor = vec4(m, 0.0, 0.0, 1.0);
}`;

  // dt = min(dt_default, dt_adaptive), which is how a production code is run: a
  // fixed timestep chosen well inside the stability limit, with the Courant
  // condition as a floor that only engages when the flow gets faster than the
  // state it was chosen for.
  //
  // Two things come out of that. The pace of the picture stops depending on the
  // flow's peak Mach, so it no longer wanders by fifteen per cent as the peak
  // wanders -- it is exactly fixed while nothing is disturbing the box. And an
  // interaction is free to raise the amplitude, because the margin absorbs it: at
  // the design state this runs at an effective Courant number of 0.32, so ctot can
  // rise by a factor of 2.5 before dt has to move at all, and beyond that the
  // adaptive branch takes over and the scheme stays stable. Defending the solver
  // by throttling the interaction had it backwards.
  const DTDX = `
uniform sampler2D uCmax;
uniform float cfl, smallc, dtdxMax;
float stepDtDx(){
  float ctot = max(texelFetch(uCmax, ivec2(0, 0), 0).x, 1e-20);
  return min(dtdxMax, min(cfl / smallc, cfl / ctot));
}
`;

  // ---- the solver: one unsplit pass ---------------------------------------
  //
  // Each face is solved twice, once for each cell that shares it. That is
  // wasteful on a CPU and right here: storing a flux would cost a second
  // full-screen pass, and on this GPU a pass costs far more than the arithmetic
  // it saves.
  const F_GODUNOV = HEAD + EOS + DTDX + `
uniform sampler2D uU;
uniform ivec2 size;
uniform float llf, dyeDiss, dxCell;

vec4 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0);
}

// HLL for a pair of primitive states whose .y is the face-normal velocity.
// No sqrt: an isothermal sound speed is a uniform.
vec4 hll(vec4 qL, vec4 qR){
  float SL, SR;
  if (llf > 0.5) {
    float sm = max(abs(qL.y), abs(qR.y)) + cs;
    SL = -sm; SR = sm;
  } else {
    SL = min(min(qL.y, qR.y) - cs, 0.0);
    SR = max(max(qL.y, qR.y) + cs, 0.0);
  }
  vec4 UL = toCons(qL), UR = toCons(qR);
  // .w is the passive scalar: flux u * (rho Y), which is what scalar_flux
  // reduces to for HLL -- componentwise HLL on rho*Y.
  vec4 FL = vec4(UL.y, qL.y * UL.y + qL.x * cs2, qL.y * UL.z, qL.y * UL.w);
  vec4 FR = vec4(UR.y, qR.y * UR.y + qR.x * cs2, qR.y * UR.z, qR.y * UR.w);
  return (SR * FL - SL * FR + SR * SL * (UR - UL)) / (SR - SL);
}

vec4 swapv(vec4 a){ return vec4(a.x, a.z, a.y, a.w); }

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4 U0 = cell(c);
  vec4 q0 = toPrim(U0);
  vec4 qw = toPrim(cell(c - ivec2(1, 0)));
  vec4 qe = toPrim(cell(c + ivec2(1, 0)));
  vec4 qs = toPrim(cell(c - ivec2(0, 1)));
  vec4 qn = toPrim(cell(c + ivec2(0, 1)));

  vec4 Fw = hll(qw, q0);
  vec4 Fe = hll(q0, qe);
  vec4 Fs = swapv(hll(swapv(qs), swapv(q0)));
  vec4 Fn = swapv(hll(swapv(q0), swapv(qn)));

  float dtdx = stepDtDx();
  vec4  U    = U0 + ((Fw - Fe) + (Fs - Fn)) * dtdx;
  U.x = max(U.x, smallr);
  // The dye needs a sink or it saturates; this is the dyeDiss knob the
  // incompressible build carried, applied to the scalar only.
  U.w = U.w / (1.0 + dyeDiss * dtdx * dxCell);
  outColor = U;
}`;

  // ---- second order: MC slopes, trace2d, HLLC -----------------------------
  //
  // The first-order scheme above cannot hold a shear layer. HLL has two waves, so
  // it has no contact: it averages the transverse velocity and the dye across the
  // interface, and a shear layer is nothing but a tangential discontinuity in
  // exactly those variables. Measured, a Kelvin-Helmholtz layer deposited at
  // delta_u = 1.3 cs decayed with an e-folding time of 2.5 s, about the same as its
  // own growth time, so it never completed a roll.
  //
  // This is the same scheme one configuration up, and the configuration is the
  // reference's own:
  //
  //   slope_type = 2   MonCen, the monotonised central difference, from uslope
  //   trace2d          the primitive variables advanced half a step by their own
  //                    Jacobian, transverse terms included, so the interface
  //                    states are centred in time as well as in space
  //   HLLC             three waves, from hllc_fluxes: the contact is resolved, so
  //                    v and Y are taken from whichever side of it the interface
  //                    lies on rather than blended across it
  //
  // One flux evaluation per step, as in the reference -- the half-step predictor is
  // what buys second-order time accuracy, not a second pass. The stencil is
  // thirteen cells: the five-cell cross, each of whose slopes needs its own
  // neighbours.
  const F_GODUNOV2 = HEAD + EOS + DTDX + `
uniform sampler2D uU;
uniform ivec2 size;
uniform float dyeDiss, dxCell, slopeType, useHllc;

vec4 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0);
}
vec4 prim(ivec2 c){ return toPrim(cell(c)); }

// uslope with slope_type = 2. dlft/drgt are the one-sided differences scaled by
// slope_type, dcen the central one; the limiter takes the smaller of the two
// one-sided slopes, and zero if they disagree in sign. slope_type = 0 returns zero
// and the whole scheme collapses to piecewise constant, which is what the
// reference does with it.
vec4 uslope(vec4 qm, vec4 q0, vec4 qp){
  vec4 dlft = slopeType * (q0 - qm);
  vec4 drgt = slopeType * (qp - q0);
  vec4 dcen = 0.5 * (dlft + drgt) / max(slopeType, 1.0);
  // step() is 1 where the product is >= 0; where it is exactly zero one of the
  // one-sided slopes is zero and the min below is zero anyway
  vec4 dlim = min(abs(dlft), abs(drgt)) * step(vec4(0.0), dlft * drgt);
  return sign(dcen) * min(dlim, abs(dcen));
}

// trace2d, isothermal. The pressure equation drops out -- p = rho cs^2 means
// dp = cs^2 drho -- and what is left is the reference's sr0/su0/sv0 with the
// passive scalar carried the same way. hx and hy are half slopes, as drx and dry
// are there.
vec4 trace(vec4 q, vec4 hx, vec4 hy, float dtdx){
  float r = max(q.x, smallr), u = q.y, v = q.z;
  float cs2r = cs2 / r;
  vec4 s;
  s.x = -u * hx.x - v * hy.x - (hx.y + hy.z) * r;
  s.y = -u * hx.y - v * hy.y - cs2r * hx.x;
  s.z = -u * hx.z - v * hy.z - cs2r * hy.x;
  s.w = -u * hx.w - v * hy.w;
  vec4 o = q + s * dtdx;
  // the reference falls back to the unpredicted value if the predictor undershoots
  // the density floor, rather than clamping to it
  o.x = o.x < smallr ? q.x : o.x;
  return o;
}

vec4 hll(vec4 qL, vec4 qR){
  float SL = min(min(qL.y, qR.y) - cs, 0.0);
  float SR = max(max(qL.y, qR.y) + cs, 0.0);
  vec4 UL = toCons(qL), UR = toCons(qR);
  vec4 FL = vec4(UL.y, qL.y * UL.y + qL.x * cs2, qL.y * UL.z, qL.y * UL.w);
  vec4 FR = vec4(UR.y, qR.y * UR.y + qR.x * cs2, qR.y * UR.z, qR.y * UR.w);
  return (SR * FL - SL * FR + SR * SL * (UR - UL)) / (SR - SL);
}

// hllc_fluxes. Wave speeds are Davis estimates and are deliberately *not* clamped
// through zero the way the HLL ones are: the four branches need their real signs.
vec4 hllc(vec4 qL, vec4 qR){
  float rl = max(qL.x, smallr), rr = max(qR.x, smallr);
  float ul = qL.y, ur = qR.y;
  float pl = rl * cs2, pr = rr * cs2;
  float sl = min(ul, ur) - cs;
  float sr = max(ul, ur) + cs;
  float rcl = rl * (ul - sl), rcr = rr * (sr - ur);   // both >= rho cs > 0
  float inv = 1.0 / (rcr + rcl);
  float ustar = (rcr * ur + rcl * ul + (pl - pr)) * inv;
  float pstar = (rcr * pl + rcl * pr + rcl * rcr * (ul - ur)) * inv;
  float ro, uo, vo, po, yo;
  if (sl > 0.0) {
    ro = rl; uo = ul; vo = qL.z; po = pl; yo = qL.w;
  } else if (ustar > 0.0) {
    // guards: the denominators vanish only where a wave speed coincides with the
    // contact, which is survivable in double precision and not in single
    ro = rl * (sl - ul) / min(sl - ustar, -1e-8); uo = ustar; vo = qL.z; po = pstar; yo = qL.w;
  } else if (sr > 0.0) {
    ro = rr * (sr - ur) / max(sr - ustar, 1e-8); uo = ustar; vo = qR.z; po = pstar; yo = qR.w;
  } else {
    ro = rr; uo = ur; vo = qR.z; po = pr; yo = qR.w;
  }
  float fd = ro * uo;
  // scalar_flux's HLLC branch: the dye rides on the mass flux, upwinded, which is
  // what the branch above has already selected
  return vec4(fd, fd * uo + po, fd * vo, fd * yo);
}

vec4 riem(vec4 qL, vec4 qR){ return useHllc > 0.5 ? hllc(qL, qR) : hll(qL, qR); }

vec4 swapv(vec4 a){ return vec4(a.x, a.z, a.y, a.w); }

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  const ivec2 ex = ivec2(1, 0), ey = ivec2(0, 1);

  vec4 q0 = prim(c);
  vec4 qE = prim(c + ex),      qW = prim(c - ex);
  vec4 qN = prim(c + ey),      qS = prim(c - ey);
  vec4 qEE = prim(c + 2 * ex), qWW = prim(c - 2 * ex);
  vec4 qNN = prim(c + 2 * ey), qSS = prim(c - 2 * ey);
  vec4 qNE = prim(c + ex + ey), qSE = prim(c + ex - ey);
  vec4 qNW = prim(c - ex + ey), qSW = prim(c - ex - ey);

  // half slopes on the five-cell cross, both directions for each: the predictor
  // couples them
  vec4 h0x = 0.5 * uslope(qW, q0, qE),   h0y = 0.5 * uslope(qS, q0, qN);
  vec4 hEx = 0.5 * uslope(q0, qE, qEE),  hEy = 0.5 * uslope(qSE, qE, qNE);
  vec4 hWx = 0.5 * uslope(qWW, qW, q0),  hWy = 0.5 * uslope(qSW, qW, qNW);
  vec4 hNx = 0.5 * uslope(qNW, qN, qNE), hNy = 0.5 * uslope(q0, qN, qNN);
  vec4 hSx = 0.5 * uslope(qSW, qS, qSE), hSy = 0.5 * uslope(qSS, qS, q0);

  float dtdx = stepDtDx();
  vec4 p0 = trace(q0, h0x, h0y, dtdx);
  vec4 pE = trace(qE, hEx, hEy, dtdx);
  vec4 pW = trace(qW, hWx, hWy, dtdx);
  vec4 pN = trace(qN, hNx, hNy, dtdx);
  vec4 pS = trace(qS, hSx, hSy, dtdx);

  vec4 Fe = riem(p0 + h0x, pE - hEx);
  vec4 Fw = riem(pW + hWx, p0 - h0x);
  vec4 Fn = swapv(riem(swapv(p0 + h0y), swapv(pN - hNy)));
  vec4 Fs = swapv(riem(swapv(pS + hSy), swapv(p0 - h0y)));

  vec4 U = cell(c) + ((Fw - Fe) + (Fs - Fn)) * dtdx;
  U.x = max(U.x, smallr);
  U.w = U.w / (1.0 + dyeDiss * dtdx * dxCell);
  outColor = U;
}`;

  // ---- third order in space: PPM ------------------------------------------
  //
  // Not from the reference: it offers slope_type 0 through 6 and every one of them
  // is piecewise linear, so this is Colella & Woodward 1984 directly, and the PLM
  // program above -- which *is* the reference's configuration -- stays in place as
  // the alternative.
  //
  //   1. limited slopes, the same MonCen ones PLM uses
  //   2. interface values from the fourth-order interpolation built on them,
  //      a(i+1/2) = (a_i + a_i+1)/2 + (dq_i - dq_i+1)/6
  //   3. monotonisation, CW84 (1.10): flatten at an extremum, and pull an
  //      overshooting edge back until the parabola is monotone across the cell
  //   4. the face state is the parabola averaged over the slab the flow crosses in
  //      half a step, which is where the normal direction's time centring comes
  //      from. What that average does *not* cover -- the transverse terms, the
  //      compression terms, and the pressure gradient -- is still trace2d's, so
  //      every -u*hx advection term is dropped there and nothing is counted twice.
  //
  // This is PPM reconstruction with an advective time centring, not full
  // characteristic tracing: one speed for all variables rather than one per wave
  // family. At an effective Courant number of 0.32 the difference between those is
  // small; the difference between a parabola and a line is not.
  const F_GODUNOV3 = HEAD + EOS + DTDX + `
uniform sampler2D uU;
uniform ivec2 size;
uniform float dyeDiss, dxCell, slopeType, useHllc;

vec4 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0);
}
vec4 prim(ivec2 c){ return toPrim(cell(c)); }

vec4 uslope(vec4 qm, vec4 q0, vec4 qp){
  vec4 dlft = slopeType * (q0 - qm);
  vec4 drgt = slopeType * (qp - q0);
  vec4 dcen = 0.5 * (dlft + drgt) / max(slopeType, 1.0);
  vec4 dlim = min(abs(dlft), abs(drgt)) * step(vec4(0.0), dlft * drgt);
  return sign(dcen) * min(dlim, abs(dcen));
}

// The monotonised parabola for one cell, from slopes computed by the caller.
// Adjacent cells' windows overlap by two, so computing the five slopes once per
// sweep and passing them in saves eight of eighteen uslope evaluations -- and
// uslope is the most expensive thing in the pass.
void ppmEdges(vec4 qm, vec4 q0, vec4 qp, vec4 dm, vec4 d0, vec4 dp, out vec4 aL, out vec4 aR){
  aL = 0.5 * (qm + q0) + (dm - d0) * (1.0 / 6.0);
  aR = 0.5 * (q0 + qp) + (d0 - dp) * (1.0 / 6.0);
  // an extremum: flatten to the cell average, or the parabola invents a new one
  // (flat is a reserved interpolation qualifier in GLSL ES, hence the name)
  vec4 ext = step((aR - q0) * (q0 - aL), vec4(0.0));
  aL = mix(aL, q0, ext);
  aR = mix(aR, q0, ext);
  // and the two overshoot cases, which are mutually exclusive, so the second mix
  // is a no-op wherever the first fired
  vec4 d = aR - aL;
  vec4 m = q0 - 0.5 * (aL + aR);
  vec4 c1 = step(d * d * (1.0 / 6.0), d * m);
  vec4 c2 = step(d * m, -d * d * (1.0 / 6.0));
  aL = mix(aL, 3.0 * q0 - 2.0 * aR, c1);
  aR = mix(aR, 3.0 * q0 - 2.0 * aL, c2);
}

// The parabola's average over the slab of width |u| dt/2 adjoining a face.
vec4 ppmFace(vec4 q0, vec4 aL, vec4 aR, float sig, float right){
  vec4 d = aR - aL;
  vec4 a6 = 6.0 * (q0 - 0.5 * (aL + aR));
  float k = 1.0 - (2.0 / 3.0) * sig;
  return right > 0.5 ? aR - 0.5 * sig * (d - k * a6)
                     : aL + 0.5 * sig * (d + k * a6);
}

vec4 floorRho(vec4 f, vec4 q){ return f.x < smallr ? vec4(q.x, f.yzw) : f; }

// trace2d minus its normal-advection terms, which the parabola average already
// carries. hx is the parabola's effective half slope, hy the transverse one.
vec4 src(vec4 q, vec4 hx, vec4 hy){
  float r = max(q.x, smallr), v = q.z;
  float cs2r = cs2 / r;
  vec4 s;
  s.x = -v * hy.x - (hx.y + hy.z) * r;
  s.y = -v * hy.y - cs2r * hx.x;
  s.z = -v * hy.z - cs2r * hy.x;
  s.w = -v * hy.w;
  return s;
}

vec4 hll(vec4 qL, vec4 qR){
  float SL = min(min(qL.y, qR.y) - cs, 0.0);
  float SR = max(max(qL.y, qR.y) + cs, 0.0);
  vec4 UL = toCons(qL), UR = toCons(qR);
  vec4 FL = vec4(UL.y, qL.y * UL.y + qL.x * cs2, qL.y * UL.z, qL.y * UL.w);
  vec4 FR = vec4(UR.y, qR.y * UR.y + qR.x * cs2, qR.y * UR.z, qR.y * UR.w);
  return (SR * FL - SL * FR + SR * SL * (UR - UL)) / (SR - SL);
}

vec4 hllc(vec4 qL, vec4 qR){
  float rl = max(qL.x, smallr), rr = max(qR.x, smallr);
  float ul = qL.y, ur = qR.y;
  float pl = rl * cs2, pr = rr * cs2;
  float sl = min(ul, ur) - cs;
  float sr = max(ul, ur) + cs;
  float rcl = rl * (ul - sl), rcr = rr * (sr - ur);
  float inv = 1.0 / (rcr + rcl);
  float ustar = (rcr * ur + rcl * ul + (pl - pr)) * inv;
  float pstar = (rcr * pl + rcl * pr + rcl * rcr * (ul - ur)) * inv;
  float ro, uo, vo, po, yo;
  if (sl > 0.0) {
    ro = rl; uo = ul; vo = qL.z; po = pl; yo = qL.w;
  } else if (ustar > 0.0) {
    ro = rl * (sl - ul) / min(sl - ustar, -1e-8); uo = ustar; vo = qL.z; po = pstar; yo = qL.w;
  } else if (sr > 0.0) {
    ro = rr * (sr - ur) / max(sr - ustar, 1e-8); uo = ustar; vo = qR.z; po = pstar; yo = qR.w;
  } else {
    ro = rr; uo = ur; vo = qR.z; po = pr; yo = qR.w;
  }
  float fd = ro * uo;
  return vec4(fd, fd * uo + po, fd * vo, fd * yo);
}

vec4 riem(vec4 qL, vec4 qR){ return useHllc > 0.5 ? hllc(qL, qR) : hll(qL, qR); }

vec4 swapv(vec4 a){ return vec4(a.x, a.z, a.y, a.w); }

// One sweep. a3..c..b3 run along the face normal; the six t* are the transverse
// neighbours of the three cells that own the two faces. Everything arrives with
// the normal velocity in .y, so the y sweep is this same code on swapped inputs.
void sweep(vec4 am3, vec4 am2, vec4 am1, vec4 a0, vec4 ap1, vec4 ap2, vec4 ap3,
           vec4 tmm, vec4 tpm, vec4 tm0, vec4 tp0, vec4 tmp, vec4 tpp,
           float dtdx, out vec4 Flo, out vec4 Fhi){
  // the five limited slopes this sweep needs, each computed once
  vec4 d2m = uslope(am3, am2, am1);
  vec4 d1m = uslope(am2, am1, a0);
  vec4 d00 = uslope(am1, a0,  ap1);
  vec4 d1p = uslope(a0,  ap1, ap2);
  vec4 d2p = uslope(ap1, ap2, ap3);

  vec4 lL, lR, cL, cR, rL, rR;
  ppmEdges(am2, am1, a0,  d2m, d1m, d00, lL, lR);   // the cell below
  ppmEdges(am1, a0,  ap1, d1m, d00, d1p, cL, cR);   // this cell
  ppmEdges(a0,  ap1, ap2, d00, d1p, d2p, rL, rR);   // the cell above

  vec4 hyL = 0.5 * uslope(tmm, am1, tpm);
  vec4 hy0 = 0.5 * uslope(tm0, a0,  tp0);
  vec4 hyR = 0.5 * uslope(tmp, ap1, tpp);

  float sL = clamp(abs(am1.y) * dtdx, 0.0, 1.0);
  float s0 = clamp(abs(a0.y)  * dtdx, 0.0, 1.0);
  float sR = clamp(abs(ap1.y) * dtdx, 0.0, 1.0);

  vec4 eL = src(am1, 0.5 * (lR - lL), hyL) * dtdx;
  vec4 e0 = src(a0,  0.5 * (cR - cL), hy0) * dtdx;
  vec4 eR = src(ap1, 0.5 * (rR - rL), hyR) * dtdx;

  // The reference floors its reconstructed states -- trace2d falls back to the
  // unpredicted cell value when the predictor undershoots smallr rather than
  // clamping to it -- and the PLM path here does the same. This one did not, and a
  // parabola through a freshly injected discontinuity can undershoot: a face
  // density at the floor next to a large velocity makes rcl tiny in HLLC, ustar
  // enormous, and the flux garbage. That is the instability mouse driving found.
  Flo = riem(floorRho(ppmFace(am1, lL, lR, sL, 1.0) + eL, am1),
             floorRho(ppmFace(a0,  cL, cR, s0, 0.0) + e0, a0));
  Fhi = riem(floorRho(ppmFace(a0,  cL, cR, s0, 1.0) + e0, a0),
             floorRho(ppmFace(ap1, rL, rR, sR, 0.0) + eR, ap1));
}

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  const ivec2 ex = ivec2(1, 0), ey = ivec2(0, 1);

  vec4 q0 = prim(c);
  vec4 x1p = prim(c + ex),     x1m = prim(c - ex);
  vec4 x2p = prim(c + 2 * ex), x2m = prim(c - 2 * ex);
  vec4 x3p = prim(c + 3 * ex), x3m = prim(c - 3 * ex);
  vec4 y1p = prim(c + ey),     y1m = prim(c - ey);
  vec4 y2p = prim(c + 2 * ey), y2m = prim(c - 2 * ey);
  vec4 y3p = prim(c + 3 * ey), y3m = prim(c - 3 * ey);
  vec4 dNE = prim(c + ex + ey), dSE = prim(c + ex - ey);
  vec4 dNW = prim(c - ex + ey), dSW = prim(c - ex - ey);

  float dtdx = stepDtDx();

  vec4 Fw, Fe;
  sweep(x3m, x2m, x1m, q0, x1p, x2p, x3p,
        dSW, dNW, y1m, y1p, dSE, dNE, dtdx, Fw, Fe);

  // the same sweep with the transverse velocity in .y, then swapped back
  vec4 Fs, Fn;
  sweep(swapv(y3m), swapv(y2m), swapv(y1m), swapv(q0), swapv(y1p), swapv(y2p), swapv(y3p),
        swapv(dSW), swapv(dSE), swapv(x1m), swapv(x1p), swapv(dNW), swapv(dNE),
        dtdx, Fs, Fn);
  Fs = swapv(Fs); Fn = swapv(Fn);

  vec4 U = cell(c) + ((Fw - Fe) + (Fs - Fn)) * dtdx;
  U.x = max(U.x, smallr);
  U.w = U.w / (1.0 + dyeDiss * dtdx * dxCell);
  outColor = U;
}`;

  // ---- forcing ------------------------------------------------------------
  //
  // Two counter-rotating stirrers wandering the box, which is what the
  // incompressible build used and what gives this background its character. The
  // kick is mostly tangential -- solenoidal, so it makes vortices rather than
  // sound -- with a small radial part, because a purely solenoidal drive in a
  // compressible gas takes a long time to build the shock network that is the
  // whole reason for being here.
  const F_STIR = HEAD + DTDX + `
uniform sampler2D uU;
uniform vec2  c0, c1;
uniform vec3  dye;
uniform float amp, radius, dyeRadius, zeta, aspect, dxCell, smallr;

vec2 kick(vec2 p, vec2 c, float dir){
  vec2 d = p - c;
  d.x *= aspect;
  float r2 = dot(d, d);
  float w  = exp(-r2 / radius);
  float r  = sqrt(r2) + 1e-5;
  vec2 t   = vec2(-d.y, d.x) / r;          // tangential: divergence-free
  vec2 g   = d / r;                        // radial: curl-free
  return dir * w * mix(t, g, zeta);
}

void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 a = kick(vUv, c0, 1.0) + kick(vUv, c1, -1.0);
  float dt = stepDtDx() * dxCell;
  float r  = max(U.x, smallr);
  // an acceleration, so the force is rho*a and heavy gas is harder to move
  vec2 m = U.yz + r * a * amp * dt;
  // The stirrers leave a little dye behind them, from a footprint narrower than
  // the one they push with: the ink has to be laid down in a track that the flow
  // then stretches into a filament. Inject it as broadly as the momentum and the
  // box fills with a uniform wash instead, which is what a wide footprint and a
  // slow dissipation rate produced here on the first attempt.
  float w = exp(-dot((vUv - c0) * vec2(aspect, 1.0), (vUv - c0) * vec2(aspect, 1.0)) / dyeRadius)
          + exp(-dot((vUv - c1) * vec2(aspect, 1.0), (vUv - c1) * vec2(aspect, 1.0)) / dyeRadius);
  outColor = vec4(U.x, m, U.w + r * dye.x * w * dt);
}`;

  // A pointer stroke: momentum where the cursor went, plus ink.
  const F_SPLAT = HEAD + `
uniform sampler2D uU;
uniform vec2  point, delta;
uniform float aspect, radius, ink, smallr;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float w = exp(-dot(p, p) / radius);
  float r = max(U.x, smallr);
  outColor = vec4(U.x, U.yz + delta * w * r, U.w + r * ink * w);
}`;

  // A hyperbolic-tangent shear layer across mid-screen with a sinusoidal
  // transverse seed: the canonical Kelvin-Helmholtz setup, so the roll-up is real.
  // Applied to momentum, because momentum is what is conserved here.
  //
  // Whether it actually rolls up is a race, and the layer thickness decides it.
  // Growth goes as delta_u / delta, so a thick layer grows slowly: at delta =
  // 0.055 the e-folding time was about five wall seconds, while the ambient
  // turbulence tears a layer that thick apart in half that -- the shear was being
  // deposited and then erased before the instability could do anything with it.
  // At delta = 0.028 it e-folds in a little over a wall second and wins the race.
  //
  // The seed sits at the wavelength that grows fastest. For a tanh layer that is
  // k*delta ~ 0.44, and five cycles across a box of aspect ~1.78 gives k*delta =
  // 0.49, which is why the seed is five and not some rounder number.
  //
  // The amplitude is bounded by physics at the other end: a vortex sheet stops
  // being unstable once the convective Mach number delta_u / 2 cs approaches one,
  // so a harder kick would suppress the very thing it is trying to excite -- and
  // would make strong shocks on the way. One gesture's worth of shear is delta_u ~
  // 1.4 cs, a convective Mach number of 0.7: vigorous rolls, weak fronts.
  const F_SHEAR = HEAD + `
uniform sampler2D uU;
uniform vec3  band;
uniform float amp, width, seed, phase, cycles, smallr;
void main(){
  vec4  U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  float r = max(U.x, smallr);
  float y = vUv.y - 0.5;
  float env = exp(-pow(y / (2.6 * width), 2.0));
  vec2  du = vec2(amp * tanh(y / width) * env,
                  seed * sin(6.2831853 * (vUv.x * cycles + phase)) * env);
  float bw = exp(-pow(y / (1.5 * width), 2.0));
  outColor = vec4(U.x, U.yz + r * du, U.w + r * band.x * bw);
}`;

  // A blast. In an isothermal gas P = rho cs^2, so piling up density *is* piling
  // up pressure: the Riemann solver turns it into an outward shock by itself,
  // which is rather the point of having one.
  const F_BLAST = HEAD + `
uniform sampler2D uU;
uniform vec2  point;
uniform float aspect, radius, dens, kick, smallr;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float rr = length(p) + 1e-5;
  float w  = exp(-dot(p, p) / radius);
  float r0 = max(U.x, smallr);
  float r1 = U.x + dens * w;
  // the new mass arrives with the local velocity, plus an outward kick
  vec2  m  = U.yz * (r1 / r0) + kick * w * (p / rr) * r0;
  outColor = vec4(r1, m, U.w * (r1 / r0));
}`;

  const F_INIT = HEAD + `
uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float h1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}
const float TAU = 6.28318530718;
void main(){
  vec2  uv  = vUv;
  float rho = 1.0 + 0.02 * (h1(uvec2(gl_FragCoord.xy), 17u) - 0.5);
  // A few large-scale modes and a smooth dye pattern, because this scheme cannot
  // afford to spin up from rest: at these timesteps a box crossing costs
  // thousands of steps, so the page would open on a still frame and stay there
  // for a minute. Starting from something already moving is the honest fix -- the
  // alternative is pretending the warm-up is free.
  float a = sin(TAU * uv.x * 2.0), b = cos(TAU * uv.y);
  float c = sin(TAU * (uv.x + uv.y)), d = cos(TAU * (uv.x * 2.0 - uv.y));
  vec2  u = 0.8 * vec2(a * b + 0.5 * d, -c + 0.4 * a);
  // Almost no dye: the ink belongs to the stirrers, which paint filaments. A
  // smooth box-scale dye pattern here reads as a single soft cloud for as long as
  // it takes the flow to stretch it, and at this clock that is half a minute.
  float Y = 0.04 + 0.03 * sin(TAU * (uv.x * 3.0 + uv.y)) * cos(TAU * (uv.y * 2.0 - uv.x));
  outColor = vec4(rho, rho * u, rho * Y);
}`;

  // Velocity, in a format the dust can sample with hardware filtering. RGBA32F
  // is not linearly filterable without OES_texture_float_linear, and asking for
  // LINEAR on a format that cannot filter makes the texture INCOMPLETE -- every
  // read returns zero, texelFetch included. RG16F filters everywhere.
  const F_VEL = HEAD + EOS + `
uniform sampler2D uU;
void main(){
  vec4 q = toPrim(texelFetch(uU, ivec2(gl_FragCoord.xy), 0));
  outColor = vec4(q.yz, 0.0, 1.0);
}`;

  // Everything the picture needs, evaluated once per step at grid resolution:
  // the dye, and the convergence of the velocity field. Deriving these here
  // rather than in the composite is what makes a smooth reconstruction
  // affordable -- the composite runs at display resolution, some twenty times
  // more pixels, and it was reading five cells of state and dividing five times
  // per pixel to get numbers that only ever vary on the grid. Now it interpolates
  // two scalars, in a filterable format, and can afford to do it properly.
  const F_DISP = HEAD + EOS + `
uniform sampler2D uU;
uniform ivec2 size;
vec4 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0);
}
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4  q = toPrim(cell(c));
  vec4 qe = toPrim(cell(c + ivec2(1, 0))), qw = toPrim(cell(c - ivec2(1, 0)));
  vec4 qn = toPrim(cell(c + ivec2(0, 1))), qs = toPrim(cell(c - ivec2(0, 1)));
  // centred, so this is the velocity difference per cell
  float div = 0.5 * ((qe.y - qw.y) + (qn.z - qs.z));
  outColor = vec4(q.w, max(-div, 0.0), 0.0, 1.0);
}`;

  // ---- reconstruction -----------------------------------------------------
  //
  // Two different interpolants, because the two consumers want different things.
  //
  // Catmull-Rom for the picture. It interpolates -- it passes through the sample
  // values -- so a shock front stays where the solver put it and keeps its
  // amplitude, but it is C1 across cell boundaries, which is what kills the
  // facets that plain bilinear leaves when a 176-cell grid is stretched over 800
  // pixels. Sixteen point samples collapse into nine hardware-bilinear taps: for
  // each axis the middle two weights are pooled into one fetch placed off-centre
  // at w2/(w1+w2), which is exact rather than an approximation.
  const CATROM = `
vec4 texCR(sampler2D tex, vec2 uv, vec2 size){
  vec2 p  = uv * size;
  vec2 c1 = floor(p - 0.5) + 0.5;
  vec2 f  = p - c1;
  vec2 f2 = f * f, f3 = f2 * f;
  vec2 w0 = -0.5 * f3 + f2 - 0.5 * f;
  vec2 w1 =  1.5 * f3 - 2.5 * f2 + 1.0;
  vec2 w2 = -1.5 * f3 + 2.0 * f2 + 0.5 * f;
  vec2 w3 =  0.5 * f3 - 0.5 * f2;
  vec2 w12 = w1 + w2;
  vec2 t0 = (c1 - 1.0) / size;
  vec2 t3 = (c1 + 2.0) / size;
  vec2 t12 = (c1 + w2 / w12) / size;
  return (texture(tex, vec2(t0.x,  t0.y))  * w0.x +
          texture(tex, vec2(t12.x, t0.y))  * w12.x +
          texture(tex, vec2(t3.x,  t0.y))  * w3.x) * w0.y
       + (texture(tex, vec2(t0.x,  t12.y)) * w0.x +
          texture(tex, vec2(t12.x, t12.y)) * w12.x +
          texture(tex, vec2(t3.x,  t12.y)) * w3.x) * w12.y
       + (texture(tex, vec2(t0.x,  t3.y))  * w0.x +
          texture(tex, vec2(t12.x, t3.y))  * w12.x +
          texture(tex, vec2(t3.x,  t3.y))  * w3.x) * w3.y;
}
`;

  // A cubic B-spline for the grains, which is the shape function a particle-in-
  // cell code gathers with, and for the same reason: it is C2, so a grain
  // crossing a cell boundary feels no kink in its acceleration and the population
  // cannot print the grid onto its own clustering. It approximates rather than
  // interpolates -- it is a mild low-pass on the gathered field, which is the
  // price of the smoothness and is exactly what PIC pays. Four bilinear taps.
  const BSPLINE = `
vec2 gatherVel(sampler2D tex, vec2 uv, vec2 size){
  vec2 p = uv * size - 0.5;
  vec2 i = floor(p), f = p - i;
  vec2 f2 = f * f, f3 = f2 * f;
  vec2 w0 = (1.0 - 3.0 * f + 3.0 * f2 - f3) * (1.0 / 6.0);
  vec2 w1 = (4.0 - 6.0 * f2 + 3.0 * f3) * (1.0 / 6.0);
  vec2 w2 = (1.0 + 3.0 * f + 3.0 * f2 - 3.0 * f3) * (1.0 / 6.0);
  vec2 w3 = f3 * (1.0 / 6.0);
  vec2 g0 = w0 + w1, g1 = w2 + w3;
  vec2 h0 = (i - 0.5 + w1 / g0) / size;
  vec2 h1 = (i + 1.5 + w3 / g1) / size;
  return g0.y * (g0.x * texture(tex, vec2(h0.x, h0.y)).xy +
                 g1.x * texture(tex, vec2(h1.x, h0.y)).xy)
       + g1.y * (g0.x * texture(tex, vec2(h0.x, h1.y)).xy +
                 g1.x * texture(tex, vec2(h1.x, h1.y)).xy);
}
`;

  // ---- dust ---------------------------------------------------------------

  // The grain size spectrum, identical to the incompressible build: the rank is
  // hashed out of the grain's own texel, so it costs no storage and cannot fall
  // out of step with the particle state.
  const SPECTRUM = `
uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float hash1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}
float stRank(ivec2 tc){ return hash1(uvec2(tc), 7u); }
float tauMul(float s){ return exp((s * 2.0 - 1.0) * 1.6118); }
// In the Epstein regime tau ~ a, so the tau multiplier IS the size multiplier,
// normalised to 1 at the top of the distribution.
float aRel(float s){ return exp((s - 1.0) * 2.0 * 1.6118); }
`;

  const F_PART = HEAD + SPECTRUM + BSPLINE + `
uniform sampler2D uPart, uVel;
uniform vec2  stream, gridSize;
uniform float dragA, dt, aspect, brown;
uniform uint  uFrame;

void main(){
  vec4  P   = texture(uPart, vUv);
  vec2  pos = P.xy;
  vec2  vel = P.zw;
  uvec2 id  = uvec2(gl_FragCoord.xy);

  vec2 ug = gatherVel(uVel, pos, gridSize);
  // Backward Euler on a per-grain stopping time. tau varies across the
  // population but not in time, so the median grain's dt/tau is the only scalar
  // the CPU supplies; a grain rescales it by its place in the size spectrum.
  float a = dragA / tauMul(stRank(ivec2(id)));
  vel = (vel + a * ug) / (1.0 + a);

  if (brown > 0.0) {
    vel += brown * (vec2(hash1(id, uFrame + 11u), hash1(id, uFrame + 523u)) - 0.5);
  }

  // physical velocity -> uv: the box is aspect wide and 1 tall
  pos = fract(pos + vec2((vel.x + stream.x) / aspect, vel.y + stream.y) * dt);
  outColor = vec4(pos, vel);
}`;

  const V_PARTICLE = `#version 300 es
precision highp float;
uniform sampler2D uPart, uVel;
uniform int   uPW;
uniform vec2  uStream;
uniform vec2  gridSize;
uniform float uPointSize, uDriftNorm, uAlpha, uVignette, uAspect;
out float vBright;
out float vDrift;
` + SPECTRUM + BSPLINE + `
void main(){
  int   id  = gl_VertexID;
  ivec2 tc  = ivec2(id % uPW, id / uPW);
  vec4  P   = texelFetch(uPart, tc, 0);
  vec2  pos = P.xy;
  float st  = stRank(tc);

  // Drift, not speed. A grain swept along with the gas is doing nothing
  // interesting however fast it is going; a grain slipping through it is the
  // whole subject -- and in a compressible flow the slipping happens hardest at
  // the shock fronts.
  // The same gather the drag used, so the drift a grain is lit by is the drift it
  // actually felt -- reading the gas through a different interpolant here would
  // light grains for a slip they never had.
  vec2 ug = gatherVel(uVel, pos, gridSize);
  vDrift  = clamp(length(P.zw + uStream - ug) * uDriftNorm, 0.0, 1.0);

  gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);

  vec2  q  = (pos - 0.5) * vec2(uAspect, 1.0);
  float vg = 1.0 - uVignette * smoothstep(0.30, 1.05, length(q));
  vBright  = uAlpha * vg * (0.30 + 0.70 * vDrift);
  gl_PointSize = max(1.0, uPointSize * aRel(st));
}`;

  const F_PARTICLE = `#version 300 es
precision highp float;
in float vBright;
in float vDrift;
out vec4 outColor;
uniform vec3 uColor;
void main(){
  vec2  d  = gl_PointCoord * 2.0 - 1.0;
  float r2 = dot(d, d);
  if (r2 > 1.0) discard;
  float a = exp(-r2 * 2.7) * (1.0 - r2);
  vec3  c = mix(uColor, mix(uColor, vec3(1.0), 0.42), vDrift * vDrift);
  float b = a * vBright;
  outColor = vec4(c * b, b);
}`;

  // ---- composite ----------------------------------------------------------

  const F_COMPOSITE = HEAD + CATROM + `
uniform sampler2D uDisp;
uniform vec2  gridSize;
uniform vec2  uRes;
uniform vec3  uBg, uTint, uAccent;
uniform float uTime, uGrain, uDyeGain, uVignette, uShockGain, uShockLift;

uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float h1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}

void main(){
  // One smooth reconstruction, two fields: .x is the dye, .y the convergence.
  vec2 d2 = texCR(uDisp, vUv, gridSize).xy;

  // Dye, through a squared density response: faint gas stays dark and only
  // genuinely dense wisps register, so the gas sits behind the dust.
  float dens = clamp(d2.x * uDyeGain, 0.0, 1.0);
  vec3  dye  = uAccent * dens * dens;

  // faint cold gradient so the page never reads as flat black
  float g = smoothstep(1.15, -0.15, vUv.y + vUv.x * 0.22);
  vec3  c = uBg + uTint * g * 0.5;
  c += dye;

  // Convergence, which is the thing an incompressible background cannot show at
  // all: div u < 0 is a shock, and here it is a real discontinuity the Riemann
  // solver captured rather than a gradient in a projected field.
  //
  // Drawn in the chapter's own accent, lifted towards white on the strongest
  // fronts exactly as a fast-drifting grain is. An amber highlight did read as
  // "shock", but it introduced a hue that appears nowhere else on the page, so
  // the fronts stopped being part of the picture and became an annotation on it.
  float sh = clamp(d2.y * uShockGain, 0.0, 1.0);
  c += mix(uAccent, vec3(1.0), 0.30) * sh * sh * uShockLift;

  float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
  c *= 1.0 / (1.0 + 0.62 * l);

  float d = length((vUv - 0.5) * vec2(uRes.x / uRes.y, 1.0));
  c *= 1.0 - uVignette * smoothstep(0.30, 1.05, d);

  c += (h1(uvec2(gl_FragCoord.xy), uint(uTime * 60.0)) - 0.5) * uGrain;
  outColor = vec4(max(c, 0.0), 1.0);
}`;

  // Diagnostics: two maxima and two sums, carried through one reduction so the
  // whole measurement is a single readback. It used to take two -- one here and one
  // off the tip of the CFL chain for ctot -- and a readback is not a cheap thing to
  // do twice: it drains a queue that is now a hundred and eighty steps a second
  // deep, and the pair of them cost about six frames per second.
  //
  //   x = ctot, max-reduced (the Courant sum the solver steps on)
  //   y = |u|^2, summed      (rms Mach)
  //   z = rho, summed        (mass, and so a conservation check)
  //   w = |u|, max-reduced   (peak Mach, for the readout)
  const F_METRICS = HEAD + EOS + `
uniform sampler2D uU;
void main(){
  vec4  q = toPrim(texelFetch(uU, ivec2(gl_FragCoord.xy), 0));
  outColor = vec4(abs(q.y) + abs(q.z) + 2.0 * cs, dot(q.yz, q.yz), q.x, length(q.yz));
}`;

  const F_REDUCE = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 4;
  float mx = 0.0, mw = 0.0;
  vec2  s  = vec2(0.0);
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      vec4 v = texelFetch(uSrc, p, 0);
      mx = max(mx, v.x);
      mw = max(mw, v.w);
      s += v.yz;
    }
  }
  outColor = vec4(mx, s, mw);
}`;

  // ------------------------------------------------------------- gl helpers

  function compile(gl, type, src) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
      throw new Error('shader: ' + gl.getShaderInfoLog(s));
    }
    return s;
  }

  function program(gl, vsrc, fsrc) {
    const p = gl.createProgram();
    gl.attachShader(p, compile(gl, gl.VERTEX_SHADER, vsrc));
    gl.attachShader(p, compile(gl, gl.FRAGMENT_SHADER, fsrc));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
      throw new Error('link: ' + gl.getProgramInfoLog(p));
    }
    const u = {};
    const n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
    for (let i = 0; i < n; i++) {
      const nm = gl.getActiveUniform(p, i).name.replace(/\[0\]$/, '');
      u[nm] = gl.getUniformLocation(p, nm);
    }
    return { p, u };
  }

  function makeFBO(gl, w, h, internal, format, type, filter, wrap) {
    const tex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, wrap);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, wrap);
    gl.texImage2D(gl.TEXTURE_2D, 0, internal, w, h, 0, format, type, null);
    const fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) {
      throw new Error('incomplete framebuffer ' + internal);
    }
    return {
      tex, fbo, w, h,
      bind(unit) { gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, tex); return unit; }
    };
  }

  // ------------------------------------------------------------- presets
  //
  // The same chapters, the same accents. What changed is the meaning of the
  // velocity-like numbers: the incompressible build worked in cells per second
  // with |u| ~ 100, this one works in Mach numbers with cs = 1, so anything with
  // units of speed is about forty times smaller. tau is unchanged, because it
  // was already in units of an eddy turnover.

  const BASE = {
    accent: [0.62, 0.72, 0.98],
    tint: [0.020, 0.024, 0.040],
    tau: 0.10,
    // Target rms Mach, scaled by stirGain. Transonic on purpose: at M ~ 0.8 the
    // rms flow is subsonic, so vortices are not shredded as fast as they form and
    // the field keeps the swirl of the incompressible build, while the tail of the
    // distribution -- peaks run 2-3x the rms -- is supersonic and steepens into
    // real shocks. Whirls and shocks. The earlier M ~ 2.75 was a shock-dominated
    // flow, which is why it read as violent: at that Mach the density contrast is
    // of order M^2 and eddies are gone within a crossing time of forming.
    // 0.55 rather than 0.80: at rms Mach 0.7 the density contrast is of order a
    // half, the fronts read as structure rather than as events, and there is
    // enough margin left that a scroll cannot push the flow supersonic on average.
    mach: 0.55,
    dustGain: 1.00,
    dyeGain: 0.75,
    dyeDiss: 0.75,
    pointSize: 1.5,
    // Coefficient on tau * M^2 -- see the drift reference in drawDust.
    driftRef: 2.0,
    brown: 0.0,
    bhat: [1.0, 0.0],
    stream: 0.0,        // grain bulk streaming, in units of cs
    stirGain: 1.00,
    // Compressive fraction of the stirring. Low, because a solenoidal drive is
    // what makes eddies; the compressive part of the flow can come from the
    // nonlinear steepening of those eddies instead of being injected directly,
    // which is both cheaper visually and closer to how it happens.
    zeta: 0.14,
    grain: 0.030,
    vignette: 0.62
  };

  const PRESETS = {
    hero:      { accent: [0.50, 0.64, 1.00], tau: 0.14, dustGain: 1.48, stirGain: 1.25 },
    bio:       { accent: [0.96, 0.80, 0.50], tau: 0.10, dustGain: 1.05 },
    picdust:   { accent: [1.00, 0.58, 0.18], tau: 0.030, dustGain: 1.45, pointSize: 1.2, stirGain: 1.15 },
    dfmm:      { accent: [0.26, 0.95, 0.72], tau: 0.85,  dustGain: 2.24, pointSize: 1.75, brown: 0.08 },
    mhd:       { accent: [0.60, 0.46, 1.00], tau: 0.12,  dustGain: 1.81, bhat: [0.94, 0.34], zeta: 0.08 },
    cosmicray: { accent: [0.22, 0.86, 1.00], tau: 0.020, dustGain: 2.40, stream: 1.5, bhat: [0.94, 0.34], pointSize: 1.15 },
    phrike:    { accent: [1.00, 0.24, 0.48], tau: 0.06,  dustGain: 1.38, stirGain: 1.35, zeta: 0.26 },
    papers:    { accent: [0.94, 0.82, 0.52], tau: 0.09,  dustGain: 0.95, stirGain: 0.85 },
    blog:      { accent: [0.56, 0.72, 1.00], tau: 0.11,  dustGain: 0.95, stirGain: 0.8 },
    iron:      { accent: [1.00, 0.13, 0.10], tau: 0.55,  dustGain: 1.00, pointSize: 1.6, stirGain: 0.65, vignette: 0.74, zeta: 0.28 },
    cactus:    { accent: [0.54, 0.84, 0.32], tau: 0.95,  dustGain: 2.06, pointSize: 1.7, brown: 0.16, stirGain: 0.42, grain: 0.040 },
    skate:     { accent: [1.00, 0.48, 0.06], tau: 0.20,  dustGain: 1.19, pointSize: 1.3, stirGain: 1.1 },
    contact:   { accent: [0.92, 0.86, 0.58], tau: 0.10,  dustGain: 1.05, stirGain: 0.9 }
  };

  function resolvePreset(key) {
    return Object.assign({}, BASE, PRESETS[key] || PRESETS.hero);
  }

  // ---------------------------------------------------------------- engine

  function Field(canvas) {
    const gl = canvas.getContext('webgl2', {
      alpha: false, depth: false, stencil: false, antialias: false,
      premultipliedAlpha: false, preserveDrawingBuffer: false,
      powerPreference: 'default'
    });
    if (!gl) throw new Error('no webgl2');
    if (!gl.getExtension('EXT_color_buffer_float')) throw new Error('no float render targets');

    const reduced = global.matchMedia && global.matchMedia('(prefers-reduced-motion: reduce)').matches;
    const mobile = Math.min(global.innerWidth, global.innerHeight) < 720;

    const PART_SIDE = 320;
    const CBLOCK = 8;
    // 1e-4 was reachable: at Mach 5 the rarefactions span five decades and pin
    // the floor. Three decades below the mean is still far outside the physics.
    const SMALLR = 1e-3;
    const SMALLC = 1e-3;
    const CFL = 0.8;                 // RAMSES courant_factor
    const CS = 1.0;
    // Solver configuration, in the reference's own terms. slope_type 0 makes the
    // second-order program collapse to piecewise constant, so the four combinations
    // of (slope_type, riemann solver) can be compared through one code path.
    const SLOPE_TYPE = 2;            // 0 = piecewise constant, 2 = MonCen
    const USE_HLLC = true;
    // Reconstruction: 'ppm' or 'plm'. Backend only, no control -- PLM is the
    // reference's own configuration and is kept as the fallback.
    const RECON = 'ppm';
    // ...and it is also what runs while the box is being driven, even when PPM is
    // selected. An interaction drops a discontinuity into the middle of the flow,
    // and a parabola fitted across a fresh one is exactly where PPM is fragile: the
    // monotonisation acts on each primitive variable independently, so a state that
    // is monotone in every variable separately can still be unphysical taken
    // together. MonCen slopes are the more robust thing to have under a mouse. The
    // window covers the injection and the early roll-up, then PPM comes back.
    const PLM_WINDOW = 1.6;          // wall seconds after any interaction
    // The default timestep, as dt/dx so it is independent of the tier's grid. At
    // the design state (rms Mach 0.7, ctot ~ 4.3) the Courant condition would allow
    // 0.186, so this is a factor of 2.5 inside it.
    const DTDX_MAX = 0.075;

    // [gridH, gasStepsPerFrame, grainFrac, dprCap].
    //
    // The work per frame is a fixed integer number of Godunov steps, never a
    // while-loop chasing a sim-time target: a fixed cost is what makes the frame
    // rate stable, and the varying quantity is then how much sim time a frame
    // buys. The tier ladder below is an emergency brake, not a working range --
    // at one step per frame the top rung has better than 2x headroom on this
    // class of GPU, so it should essentially never engage.
    //
    // How fast the picture moves is set here rather than left to chance. In box
    // heights per wall second it is
    //
    //   rate = M * fps * n * DTDX_MAX / gridH   ~   3.3 * n / gridH
    //
    // at the rms Mach 0.7 this page holds. The ladder keeps n/gridH constant, so
    // every rung moves at the same apparent speed. The incompressible build,
    // measured over a minute, runs at 0.057 box heights per wall second; this lands
    // at 0.060, and because the default timestep is fixed rather than tied to the
    // flow's peak, the residual wander is only the wander in M itself -- about
    // eight per cent, against twenty-five when dt came from the Courant condition
    // on every step.
    //
    // Note what this does to the earlier version's conclusion. Five steps per
    // frame at gridH = 112 was chosen to cover ground fast, and it did: 0.70 box
    // heights per second, twelve times the main page, which is what made it look
    // violent. Slowing down to the main page's pace hands back the entire budget,
    // so the gas can be *finer* than the version that was rushing. The acoustic
    // CFL penalty is only crippling if you insist on fast advective motion.
    const TIERS = [
      [176, 3.08, 1.00, 1.35],
      [144, 2.52, 0.85, 1.20],
      [112, 1.96, 0.65, 1.05],
      [ 88, 1.54, 0.45, 1.00]
    ];
    let tier = mobile ? 2 : 0;
    let ceiling = 0;

    const P = {
      god:    program(gl, VERT, F_GODUNOV),
      god2:   program(gl, VERT, F_GODUNOV2),
      god3:   program(gl, VERT, F_GODUNOV3),
      cmax0:  program(gl, VERT, F_CMAX_STATE),
      cmaxN:  program(gl, VERT, F_CMAX_DOWN),
      stir:   program(gl, VERT, F_STIR),
      splat:  program(gl, VERT, F_SPLAT),
      shear:  program(gl, VERT, F_SHEAR),
      blast:  program(gl, VERT, F_BLAST),
      init:   program(gl, VERT, F_INIT),
      vel:    program(gl, VERT, F_VEL),
      disp:   program(gl, VERT, F_DISP),
      part:   program(gl, VERT, F_PART),
      pdraw:  program(gl, V_PARTICLE, F_PARTICLE),
      comp:   program(gl, VERT, F_COMPOSITE),
      metric: program(gl, VERT, F_METRICS),
      reduce: program(gl, VERT, F_REDUCE)
    };

    const quad = gl.createVertexArray();
    gl.bindVertexArray(quad);
    const vb = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vb);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    const emptyVAO = gl.createVertexArray();

    let U = null, vel = null, disp = null, part = null, cmax = [], met = null, red = [];
    let grid = { w: 0, h: 0 }, pSide = 0, nPart = 0, nDraw = 0, dpr = 1, aspect = 1;
    let owned = [], persistent = [];

    const state = {
      preset: resolvePreset('hero'),
      target: resolvePreset('hero'),
      live: null,
      blend: 1,
      wall: 0, time: 0, steps: 0,
      running: true, fps: 60,
      amp: 4.0,
      // ctot is the Courant sum the shader uses, max over the grid: it is read
      // back with the metrics so the clock JS keeps matches the one the solver
      // actually integrated with.
      ctot: 2 * CS, maxSig: 1, machRms: 0, machMax: 0, dens: 1, div: 0
    };
    const pending = { splats: [], shear: 0, blasts: [] };

    function activePreset() { return state.live || state.preset; }
    // The scale every interaction is measured against: the flow's own rms Mach,
    // floored so a page that has just opened still responds to a click.
    function machNow() { return Math.max(state.machRms, 0.5); }

    function drawQuad(f) {
      if (f) { gl.bindFramebuffer(gl.FRAMEBUFFER, f.fbo); gl.viewport(0, 0, f.w, f.h); }
      else { gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.viewport(0, 0, canvas.width, canvas.height); }
      gl.bindVertexArray(quad);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function mk(w, h, i, f, t, fil, wr) {
      const o = makeFBO(gl, w, h, i, f, t, fil, wr);
      owned.push(o);
      return o;
    }
    function mkDouble(w, h, i, f, t, fil, wr) {
      let a = mk(w, h, i, f, t, fil, wr), b = mk(w, h, i, f, t, fil, wr);
      return { w, h, get read() { return a; }, get write() { return b; }, swap() { const s = a; a = b; b = s; } };
    }

    function eos(pr) {
      gl.uniform1f(pr.u.cs, CS);
      gl.uniform1f(pr.u.cs2, CS * CS);
      gl.uniform1f(pr.u.smallr, SMALLR);
    }
    // The spin-up runs at the Courant limit rather than the default step. Nothing
    // is on screen yet, so the fixed pace buys nothing there, and it would cost 2.5
    // times as many steps to reach the same developed state -- the whole stall,
    // multiplied by two and a half.
    let warming = false;
    function dtdxCap() { return warming ? 1e9 : DTDX_MAX; }
    function dtUniforms(pr) {
      gl.uniform1f(pr.u.cfl, CFL);
      gl.uniform1f(pr.u.smallc, SMALLC);
      gl.uniform1f(pr.u.dtdxMax, dtdxCap());
    }

    let allocW = -1, allocH = -1, allocTier = -1;

    function allocate(force) {
      const t = TIERS[tier];
      dpr = Math.min(global.devicePixelRatio || 1, t[3]);
      const cw = Math.max(1, Math.round(canvas.clientWidth * dpr));
      const ch = Math.max(1, Math.round(canvas.clientHeight * dpr));
      if (!force && cw === allocW && ch === allocH && tier === allocTier) return;
      allocW = cw; allocH = ch; allocTier = tier;

      const prev = owned;
      owned = [];
      canvas.width = cw; canvas.height = ch;
      aspect = cw / ch;

      grid.h = t[0];
      grid.w = Math.max(48, Math.round(t[0] * aspect / 4) * 4);

      const R = gl.REPEAT, C = gl.CLAMP_TO_EDGE, L = gl.LINEAR, N = gl.NEAREST;
      U = mkDouble(grid.w, grid.h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R);
      vel = mk(grid.w, grid.h, gl.RG16F, gl.RG, gl.HALF_FLOAT, L, R);
      // Half float, and LINEAR: both reconstructions below are built out of
      // hardware bilinear taps, so the format has to be one that filters. Asking
      // for LINEAR on RGBA32F makes the texture incomplete and every read -- even
      // texelFetch -- returns zero, silently.
      disp = mk(grid.w, grid.h, gl.RG16F, gl.RG, gl.HALF_FLOAT, L, R);
      met = mk(grid.w, grid.h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);

      cmax = [];
      let a = grid.w, b = grid.h;
      do {
        a = Math.max(1, Math.ceil(a / CBLOCK));
        b = Math.max(1, Math.ceil(b / CBLOCK));
        cmax.push(mk(a, b, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C));
      } while (a > 1 || b > 1);

      red = [];
      let rw = grid.w, rh = grid.h;
      while (rw > 4 || rh > 4) {
        rw = Math.max(1, Math.ceil(rw / 4));
        rh = Math.max(1, Math.ceil(rh / 4));
        red.push(mk(rw, rh, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C));
      }

      // Grains are allocated once for the life of the context and never
      // re-scattered; a tier change draws a shorter prefix of the same
      // population, which is a uniform random subsample because grain i is
      // scattered independently of i.
      if (!part) {
        pSide = PART_SIDE;
        nPart = pSide * pSide;
        part = mkDouble(pSide, pSide, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);
        persistent = owned.splice(owned.length - 2, 2);
        seedParticles();
      }
      nDraw = Math.max(4096, Math.round(nPart * t[2]));

      for (let i = 0; i < prev.length; i++) { gl.deleteFramebuffer(prev[i].fbo); gl.deleteTexture(prev[i].tex); }

      reset();
      warmUp();
      writeVel();
      paint(activePreset());
    }

    // Spin the seeded initial condition up into a developed flow before anything
    // is shown. A box crossing costs about ctot/(C M) * gridH ~ 1000 steps, and at
    // one step per frame that is sixteen seconds of watching a smooth sinusoid
    // shear itself apart -- so it happens here instead, in under a second, and the
    // page opens on turbulence rather than on an initial condition.
    //
    // The servo runs during the spin-up too, with a gain it would ring at while
    // anyone was watching. Without that the warm-up ends wherever the arbitrary
    // starting amplitude took it, and the first ten seconds on screen are the
    // driving settling down rather than the flow doing anything.
    // 700 steps, about two thirds of a box crossing, and seven servo updates.
    // This is a synchronous stall on the main thread, so its length is a real
    // cost: 1200 steps with a servo update every 50 measured 2.1 s, of which the
    // readbacks were the larger half -- each one drains the queue. What is left
    // finishes developing on screen, from a field that is already turbulent.
    function warmUp() {
      warming = true;
      for (let i = 0; i < 700; i++) {
        stepGas();
        if (i % 100 === 99) { measure(); servo(0.7); }
      }
      warming = false;
      cflReduce();
    }

    function reset() {
      gl.useProgram(P.init.p);
      drawQuad(U.read);
      drawQuad(U.write);
      state.time = 0; state.steps = 0;
      cflReduce();
      const probe = new Float32Array(4);
      gl.bindFramebuffer(gl.FRAMEBUFFER, U.read.fbo);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, probe);
      if (!(probe[0] > 0.1)) throw new Error('state reads back as ' + probe[0]);
    }

    function seedParticles() {
      const data = new Float32Array(nPart * 4);
      for (let i = 0; i < nPart; i++) {
        data[i * 4] = Math.random();
        data[i * 4 + 1] = Math.random();
      }
      for (const f of [part.read, part.write]) {
        gl.bindTexture(gl.TEXTURE_2D, f.tex);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, pSide, pSide, 0, gl.RGBA, gl.FLOAT, data);
      }
    }

    // ------------------------------------------------------------- stepping

    function cflReduce() {
      gl.useProgram(P.cmax0.p);
      eos(P.cmax0);
      gl.uniform1i(P.cmax0.u.uSrc, U.read.bind(0));
      gl.uniform2i(P.cmax0.u.srcSize, grid.w, grid.h);
      drawQuad(cmax[0]);
      for (let i = 1; i < cmax.length; i++) {
        gl.useProgram(P.cmaxN.p);
        gl.uniform1i(P.cmaxN.u.uSrc, cmax[i - 1].bind(0));
        gl.uniform2i(P.cmaxN.u.srcSize, cmax[i - 1].w, cmax[i - 1].h);
        drawQuad(cmax[i]);
      }
    }
    function cmaxTex() { return cmax[cmax.length - 1]; }
    // The same expression the shader evaluates, on the last ctot read back from
    // the reduction. Reconstructing it from the rms Mach instead was wrong by
    // tens of percent -- |ux|+|uy| is not |u| -- and every bit of that error went
    // straight into state.time and into the dust's timestep.
    function dtdxNow() {
      return Math.min(dtdxCap(), Math.min(CFL / SMALLC, CFL / Math.max(state.ctot, 1e-20)));
    }
    function dtNow() { return dtdxNow() / grid.h; }

    function applyForcing(pr) {
      const dx = 1 / grid.h;

      // two counter-rotating stirrers, wandering as they did before
      const t = state.time;
      const c0 = [0.5 + 0.31 * Math.cos(t * 1.7 * 0.11 + 0), 0.5 + 0.27 * Math.sin(t * 1.31 * 0.11)];
      const c1 = [0.5 + 0.31 * Math.cos(t * 1.7 * 0.162 + 1), 0.5 + 0.27 * Math.sin(t * 1.31 * 0.162 + 2)];
      const a = pr.accent;
      gl.useProgram(P.stir.p);
      dtUniforms(P.stir);
      gl.uniform1i(P.stir.u.uU, U.read.bind(0));
      gl.uniform1i(P.stir.u.uCmax, cmaxTex().bind(1));
      gl.uniform2f(P.stir.u.c0, c0[0], c0[1]);
      gl.uniform2f(P.stir.u.c1, c1[0], c1[1]);
      gl.uniform1f(P.stir.u.amp, state.amp);
      gl.uniform1f(P.stir.u.radius, 0.055);
      gl.uniform1f(P.stir.u.dyeRadius, 0.020);
      gl.uniform1f(P.stir.u.zeta, pr.zeta);
      gl.uniform1f(P.stir.u.aspect, aspect);
      gl.uniform1f(P.stir.u.dxCell, dx);
      gl.uniform1f(P.stir.u.smallr, SMALLR);
      // Ink rate. What sets the concentration in a plume is not the dissipation
      // rate but the residence time -- the dye is carried out of the footprint in
      // about (footprint / u), so the peak is rate x that, and 0.85 left the plumes
      // three times too faint to see. Dissipation sets the *length* of the wisp
      // instead: a trail survives u / dyeDiss ~ half a box before it fades, which
      // is what keeps the mean level dark while the plumes stay bright.
      gl.uniform3f(P.stir.u.dye, 2.1, 0.0, 0.0);
      drawQuad(U.write); U.swap();

      for (const s of pending.splats) {
        gl.useProgram(P.splat.p);
        gl.uniform1i(P.splat.u.uU, U.read.bind(0));
        gl.uniform2f(P.splat.u.point, s.x, s.y);
        gl.uniform2f(P.splat.u.delta, s.dx, s.dy);
        gl.uniform1f(P.splat.u.aspect, aspect);
        gl.uniform1f(P.splat.u.radius, s.radius);
        gl.uniform1f(P.splat.u.ink, s.ink);
        gl.uniform1f(P.splat.u.smallr, SMALLR);
        drawQuad(U.write); U.swap();
      }
      pending.splats.length = 0;

      if (Math.abs(pending.shear) > 1e-4) {
        // The entry point already bounds this to 3 M; the clamp is a backstop.
        const amp = Math.max(-6, Math.min(6, pending.shear));
        const full = 1.5 * machNow();
        gl.useProgram(P.shear.p);
        gl.uniform1i(P.shear.u.uU, U.read.bind(0));
        gl.uniform1f(P.shear.u.amp, amp);
        // Three cells. At first order with HLL this would have been pointless --
        // the scheme's own diffusion smeared a five-cell layer faster than the
        // instability could grow in it -- but MonCen slopes and a resolved contact
        // support a layer this thin, and growth goes as delta_u / delta. The seed
        // then sits at the wavelength that grows fastest for it: k*delta ~ 0.44 for
        // a tanh layer, which at delta = 3/gridH is eight cycles across the box.
        const w = Math.max(3 / grid.h, 0.012);
        gl.uniform1f(P.shear.u.width, w);
        gl.uniform1f(P.shear.u.cycles, Math.round(0.44 * aspect / (6.2831853 * w)));
        gl.uniform1f(P.shear.u.seed, amp * 0.50);
        gl.uniform1f(P.shear.u.phase, state.time * 0.21);
        gl.uniform1f(P.shear.u.smallr, SMALLR);
        // The dye band is not decoration: it is the tracer that makes a roll-up
        // visible at all. A velocity field has no colour, and the fronts overlay
        // only shows convergence -- without ink in the layer, a textbook KH spiral
        // is invisible on screen.
        const s = Math.min(0.9, Math.abs(amp) / full) * 0.60;
        gl.uniform3f(P.shear.u.band, a[0] * s, a[1] * s, a[2] * s);
        drawQuad(U.write); U.swap();
        pending.shear *= 0.55;
        if (Math.abs(pending.shear) < 1e-3) pending.shear = 0;
      }

      for (const b of pending.blasts) {
        gl.useProgram(P.blast.p);
        gl.uniform1i(P.blast.u.uU, U.read.bind(0));
        gl.uniform2f(P.blast.u.point, b.x, b.y);
        gl.uniform1f(P.blast.u.aspect, aspect);
        gl.uniform1f(P.blast.u.radius, b.radius);
        gl.uniform1f(P.blast.u.dens, b.dens);
        gl.uniform1f(P.blast.u.kick, b.kick);
        gl.uniform1f(P.blast.u.smallr, SMALLR);
        drawQuad(U.write); U.swap();
      }
      pending.blasts.length = 0;
    }

    function stepGas() {
      const pr = activePreset();
      // Once per step, ahead of the source terms -- the order courant_fine and
      // amr_step run in. Reducing again after the forcing cost three more passes
      // for a correction the servo keeps small anyway.
      cflReduce();
      applyForcing(pr);
      const G = (RECON === 'ppm' && agitated <= 0) ? P.god3 : P.god2;
      gl.useProgram(G.p);
      eos(G);
      dtUniforms(G);
      gl.uniform1i(G.u.uU, U.read.bind(0));
      gl.uniform1i(G.u.uCmax, cmaxTex().bind(1));
      gl.uniform2i(G.u.size, grid.w, grid.h);
      gl.uniform1f(G.u.slopeType, SLOPE_TYPE);
      gl.uniform1f(G.u.useHllc, USE_HLLC ? 1 : 0);
      gl.uniform1f(G.u.dyeDiss, pr.dyeDiss);
      gl.uniform1f(G.u.dxCell, 1 / grid.h);
      drawQuad(U.write); U.swap();
      state.time += dtNow();
      state.steps++;
    }

    function writeVel() {
      gl.useProgram(P.vel.p);
      eos(P.vel);
      gl.uniform1i(P.vel.u.uU, U.read.bind(0));
      drawQuad(vel);

      gl.useProgram(P.disp.p);
      eos(P.disp);
      gl.uniform1i(P.disp.u.uU, U.read.bind(0));
      gl.uniform2i(P.disp.u.size, grid.w, grid.h);
      drawQuad(disp);
    }

    function stepDust(pr, dt) {
      const a = dt / Math.max(pr.tau, 1e-4);
      const b = pr.bhat, bn = Math.hypot(b[0], b[1]) || 1;
      gl.useProgram(P.part.p);
      gl.uniform1i(P.part.u.uPart, part.read.bind(0));
      gl.uniform1i(P.part.u.uVel, vel.bind(1));
      gl.uniform1f(P.part.u.dragA, a);
      gl.uniform1f(P.part.u.dt, dt);
      gl.uniform1f(P.part.u.aspect, aspect);
      gl.uniform2f(P.part.u.gridSize, grid.w, grid.h);
      gl.uniform1f(P.part.u.brown, pr.brown * (dt * 60));
      gl.uniform1ui(P.part.u.uFrame, state.steps >>> 0);
      gl.uniform2f(P.part.u.stream, pr.stream * b[0] / bn, pr.stream * b[1] / bn);
      drawQuad(part.write); part.swap();
    }

    // -------------------------------------------------------------- painting

    function composite(pr) {
      gl.useProgram(P.comp.p);
      gl.uniform1i(P.comp.u.uDisp, disp.bind(0));
      gl.uniform2f(P.comp.u.gridSize, grid.w, grid.h);
      gl.uniform2f(P.comp.u.uRes, canvas.width, canvas.height);
      gl.uniform3f(P.comp.u.uBg, 0.014, 0.015, 0.021);
      gl.uniform3f(P.comp.u.uTint, pr.tint[0], pr.tint[1], pr.tint[2]);
      gl.uniform3f(P.comp.u.uAccent, pr.accent[0], pr.accent[1], pr.accent[2]);
      gl.uniform1f(P.comp.u.uTime, state.wall);
      gl.uniform1f(P.comp.u.uGrain, pr.grain);
      gl.uniform1f(P.comp.u.uDyeGain, pr.dyeGain);
      gl.uniform1f(P.comp.u.uVignette, pr.vignette);
      // The convergence scale, in velocity difference per cell. It has to track
      // the driving: a front at rms Mach 0.7 is a jump of a few tenths of cs, not
      // several times it, and a fixed scale either saturates or shows nothing.
      gl.uniform1f(P.comp.u.uShockGain, 1.0 / Math.max(0.22 * machNow(), 0.08));
      gl.uniform1f(P.comp.u.uShockLift, 0.30);
      drawQuad(null);
    }

    function drawDust(pr) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, canvas.width, canvas.height);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE);
      gl.useProgram(P.pdraw.p);
      gl.uniform1i(P.pdraw.u.uPart, part.read.bind(0));
      gl.uniform1i(P.pdraw.u.uVel, vel.bind(1));
      gl.uniform1i(P.pdraw.u.uPW, pSide);
      gl.uniform2f(P.pdraw.u.gridSize, grid.w, grid.h);
      gl.uniform1f(P.pdraw.u.uPointSize, pr.pointSize * dpr * 3.9);
      // The drift scale rides on the measured Mach number so the lighting keeps
      // its meaning as the driving changes -- but it has to ride on M^2, not M. A
      // grain of stopping time tau reaches a drift of order tau times the gas
      // acceleration, and the acceleration of a flow of rms Mach M across a box of
      // order unity is M^2 cs^2. Normalising by M alone (which is what this did,
      // tuned at M ~ 2.75) makes the dust lose all its contrast the moment the
      // driving is turned down: every grain sits at the floor of the ramp.
      const M = Math.max(state.machRms, 0.35);
      const ref = pr.driftRef * pr.tau * M * M + 1.6 * pr.stream;
      gl.uniform1f(P.pdraw.u.uDriftNorm, 1 / Math.max(ref, 0.05));
      gl.uniform1f(P.pdraw.u.uAlpha, pr.dustGain * 0.30);
      gl.uniform1f(P.pdraw.u.uVignette, pr.vignette);
      gl.uniform1f(P.pdraw.u.uAspect, aspect);
      const b = pr.bhat, bn = Math.hypot(b[0], b[1]) || 1;
      gl.uniform2f(P.pdraw.u.uStream, pr.stream * b[0] / bn, pr.stream * b[1] / bn);
      gl.uniform3f(P.pdraw.u.uColor, pr.accent[0], pr.accent[1], pr.accent[2]);
      gl.bindVertexArray(emptyVAO);
      gl.drawArrays(gl.POINTS, 0, nDraw);
      gl.disable(gl.BLEND);
    }

    function paint(pr) { composite(pr); drawDust(pr); }

    // ----------------------------------------------------------- diagnostics

    const buf = new Float32Array(4 * 64);
    let measureCountdown = 2;

    function measure() {
      gl.useProgram(P.metric.p);
      eos(P.metric);
      gl.uniform1i(P.metric.u.uU, U.read.bind(0));
      drawQuad(met);
      let src = met;
      for (const d of red) {
        gl.useProgram(P.reduce.p);
        gl.uniform1i(P.reduce.u.uSrc, src.bind(0));
        gl.uniform2i(P.reduce.u.srcSize, src.w, src.h);
        drawQuad(d);
        src = d;
      }
      gl.bindFramebuffer(gl.FRAMEBUFFER, src.fbo);
      const n = Math.min(64, src.w * src.h);
      try { gl.readPixels(0, 0, src.w, src.h, gl.RGBA, gl.FLOAT, buf); }
      catch (e) { return; }
      let mx = 0, sq = 0, sm = 0, mw = 0;
      for (let i = 0; i < n; i++) {
        mx = Math.max(mx, buf[i * 4]);
        sq += buf[i * 4 + 1];
        sm += buf[i * 4 + 2];
        mw = Math.max(mw, buf[i * 4 + 3]);
      }
      const cells = grid.w * grid.h;
      if (isFinite(mw)) state.machMax = mw / CS;
      if (isFinite(sq)) state.machRms = Math.sqrt(Math.max(0, sq / cells)) / CS;
      if (isFinite(sm) && sm > 0) state.dens = sm / cells;
      if (isFinite(mx) && mx > 0) state.ctot = mx;
      state.maxSig = state.machMax * CS + 2 * CS;

      if (!isFinite(state.machRms) || state.machMax > 400) { state.amp = 2.0; reset(); }
    }

    // Hold the target Mach. Rate-limited, and only ever called once per
    // measurement -- per step it is a runaway. The running gain is small because
    // the flow now responds on a timescale of order ten wall seconds, and a servo
    // faster than the system it drives only rings; the warm-up passes a larger
    // gain, where there is no viewer to see it hunt.
    function servo(gain) {
      const pr = activePreset();
      const want = pr.mach * pr.stirGain;
      const rel = (want - state.machRms) / Math.max(want, 0.4);
      const g = gain * Math.max(-0.25, Math.min(0.25, rel));
      state.amp = Math.max(0.05, Math.min(400, state.amp * (1 + g)));
    }

    // ------------------------------------------------------------------ loop

    function lerpPreset(dt) {
      if (state.blend >= 1) return;
      state.blend = Math.min(1, state.blend + dt * 1.35);
      const k = state.blend, e = k * k * (3 - 2 * k);
      const cur = state.preset, tgt = state.target, out = {};
      for (const key in tgt) {
        const a = cur[key], b = tgt[key];
        if (Array.isArray(b)) out[key] = b.map((v, i) => a[i] + (v - a[i]) * e);
        else if (typeof b === 'number') out[key] = a + (b - a) * e;
        else out[key] = b;
      }
      state.live = out;
    }

    const TARGET_FPS = 60;
    const FRAME_MIN = 1 / TARGET_FPS - 0.002;
    const FPS_FLOOR = 30;
    let last = 0, frames = 0, fpsT = 0, slow = 0, fast = 0, tierChanges = 0;
    let frozen = false, lost = false, dead = false;

    function retune(d) {
      if (tierChanges >= 4) return;
      const next = tier + d;
      if (next < 0 || next > TIERS.length - 1) return;
      tier = next;
      if (d > 0) ceiling = tier;
      tierChanges++;
      credit = 0;
      try { allocate(true); } catch (e) { frozen = true; dead = true; }
    }

    // One frame's worth of work, factored out so the profiler below can run it
    // off the clock.
    //
    // The steps-per-frame figure is fractional, so the coarse rungs take a step
    // only every other frame or so; that is what keeps n/gridH, and with it the
    // apparent speed, the same on every rung. The cost swing is one Godunov pass,
    // far below a frame's budget, so the frame rate does not notice.
    // The interaction budget: how much velocity recent impulses have injected, and
    // how fast that allowance comes back. In units of the rms Mach.
    const IMP_CAP = 1.0, IMP_LEAK = 0.7;
    // The pointer's own impulse budget. A splat is a velocity increment, not a
    // force, and site.js can emit one every 28 ms -- so a ceiling on each one does
    // not bound a stroke, which is the same mistake the scroll shear had. Measured
    // before this: thirty-six splats a second at the old ceiling of 1.0 M is 86
    // sound speeds a second deposited into a spot a twentieth of the box across. It
    // took the peak Mach to 8 here and 22 on the magnetised page, drove ctot past
    // the point where the Courant condition takes over from the default timestep,
    // and stayed there for seconds afterwards because nothing had removed the
    // energy: the picture visibly slowed down and did not speed back up.
    //
    // So the momentum is metered and the ink is not. The ink is a passive scalar
    // that costs the solver nothing and it is what makes a stroke visible, so it is
    // scaled by the gesture; the momentum draws on a budget that refills at a fixed
    // rate. A flick lands at a bit under half of what it used to, and a held scribble is bounded
    // by the leak rather than by how long it is held: a three-second scribble now
    // delivers a fifth of what it did, and a flick about half, with the timestep
    // measured unthrottled through both -- including a scribble held in one spot,
    // which is the case that concentrates it worst.
    const SPLAT_CAP = 0.45;      // per event, in units of the rms Mach
    const SPLAT_BUDGET = 1.60;   // how much may be outstanding
    const SPLAT_LEAK = 7.0;      // and how fast that comes back, per second
    let pimp = 0;
    let impulse = 0;
    // Wall seconds of PLM left to run. Any interaction refreshes it.
    let agitated = 0;
    function agitate() { agitated = PLM_WINDOW; }

    let credit = 0;
    function tick(pr) {
      const n = TIERS[tier][1];
      credit += n;
      let stepped = false;
      while (credit >= 1) { stepGas(); credit -= 1; stepped = true; }
      if (stepped) writeVel();
      // The dust gets the sim time a frame is *worth* rather than the time that
      // happened to be integrated during it, so it does not stutter on a frame
      // that skipped the gas. Backward Euler is L-stable for any step.
      stepDust(pr, Math.max(n * dtNow(), 1e-9));
      paint(pr);
    }

    // Offline profiling hook, reached as field.__bench(n). rAF caps at 60 Hz, so
    // fps alone cannot tell a frame that costs 4 ms from one that costs 16; and
    // gl.finish() does not serialise on ANGLE/Metal, nor are timer queries
    // coherent there. A batch of whole frames bracketed by a readPixels -- which
    // does block -- divided by the batch size, is the only honest number.
    function bench(n) {
      const pr = activePreset();
      const px = new Float32Array(4);
      const sync = function () {
        gl.bindFramebuffer(gl.FRAMEBUFFER, U.read.fbo);
        gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, px);
      };
      for (let i = 0; i < 4; i++) tick(pr);
      sync();
      const t0 = performance.now();
      for (let i = 0; i < n; i++) tick(pr);
      sync();
      return (performance.now() - t0) / n;
    }

    function frame(now) {
      if (dead) return;
      requestAnimationFrame(frame);
      if (lost || gl.isContextLost()) return;
      if (!state.running) { last = 0; return; }
      if (!last) { last = now; return; }
      const el = (now - last) / 1000;
      if (el < FRAME_MIN) return;
      last = now;
      const dtWall = Math.min(el, 1 / 12);
      state.wall += dtWall;

      frames++; fpsT += dtWall;
      if (state.wall < 2.5) { frames = 0; fpsT = 0; }
      else if (fpsT > 0.75) {
        state.fps = frames / fpsT; frames = 0; fpsT = 0;
        if (state.fps < FPS_FLOOR - 6 && tier < TIERS.length - 1) {
          fast = 0; if (++slow >= 3) { slow = 0; retune(+1); }
        } else if (state.fps > FPS_FLOOR + 12 && tier > ceiling) {
          slow = 0; if (++fast >= 5) { fast = 0; retune(-1); }
        } else if (state.fps > FPS_FLOOR) { slow = 0; fast = 0; }
      }

      lerpPreset(dtWall);
      impulse = Math.max(0, impulse - IMP_LEAK * machNow() * dtWall);
      pimp = Math.max(0, pimp - SPLAT_LEAK * machNow() * dtWall);
      agitated = Math.max(0, agitated - dtWall);
      const pr = activePreset();

      if (!frozen) {
        tick(pr);
        if (--measureCountdown <= 0) { measureCountdown = 40; measure(); servo(0.10); }
      } else {
        paint(pr);
      }
      if (reduced && state.wall > 2.5) frozen = true;
    }

    canvas.addEventListener('webglcontextlost', function (e) {
      e.preventDefault(); lost = true; owned = []; persistent = []; part = null; U = null;
    }, false);

    let resizeTimer = null;
    global.addEventListener('resize', function () {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(function () { try { allocate(false); } catch (e) { dead = true; } }, 250);
    });
    document.addEventListener('visibilitychange', function () {
      state.running = !document.hidden; last = 0;
    });

    allocate(true);
    requestAnimationFrame(frame);

    return {
      kind: 'webgl2',
      // The interaction surface, in the units site.js speaks: cells per second,
      // hundreds of them, because that is what the incompressible build used.
      //
      // Every one of these is scaled by the *current* rms Mach rather than capped
      // at some absolute multiple of cs. That is not a detail. Turning a chapter
      // fires a shear of 190, which the first version of this file clamped to Mach
      // 6 -- and then applied again every frame while the pending value decayed,
      // for a cumulative twenty-odd cs of momentum dumped into a band across a
      // flow whose rms is 0.9. The page opened on that impulse and spent the next
      // twenty seconds recovering from it, which read as a much more violent
      // background than the solver actually produces. An interaction should be a
      // perturbation on the flow, so it is written as a multiple of the flow.
      splat(x, y, dx, dy, radius) {
        const M = machNow();
        const sp = Math.hypot(dx, dy);
        if (sp < 1e-6) return;
        const cap = SPLAT_CAP * M;
        // the ink follows the gesture, so the stroke stays visible even when the
        // budget has nothing left to give the gas
        const m = Math.min(1, sp / cap);
        const give = Math.min(sp, cap, Math.max(0, SPLAT_BUDGET * M - pimp));
        const k = give / sp;
        pimp += give;
        pending.splats.push({
          x, y, dx: dx * k, dy: dy * k,
          ink: 0.06 + m * 0.30,
          radius: radius || 0.0030
        });
        if (pending.splats.length > 6) pending.splats.splice(0, pending.splats.length - 6);
        frozen = false;
        agitate();
      },
      // A chapter turn sends 190. A scroll gesture sends a stream of them -- one
      // per scroll event -- and that is the thing to defend against: an impulse
      // has to be a velocity increment rather than a force, because it must
      // register within a frame or two of the gesture, and a ceiling on the
      // amplitude does not bound a sequence that keeps refilling it. Fourteen
      // events, which is one flick of a trackpad, tripled the rms Mach and drove
      // the peak to 5.2 through a 0.45 M ceiling.
      //
      // So the impulses draw on a budget that refills at a fixed rate. One turn
      // lands in full; a held scroll gets the budget and then the leak, which is
      // about one M of shear per second however hard it is scrolled. The band is a
      // tenth of the box, so that is a weak front, not a supersonic one.
      shear(amount) {
        const M = machNow();
        // Normalised on 70, not 190. 190 is what a chapter turn sends, but a scroll
        // frame sends dy*3 and a trackpad frame is dy ~ 20, so the real path was
        // arriving with about 60 -- a third of the assumed scale, giving a layer of
        // delta_u ~ 0.2 cs against an rms of 0.75. Invisible, and nowhere near able
        // to roll up. The budget below still bounds where a sustained gesture ends
        // up; this only sets how fast it gets there, and a page turn now deposits
        // most of a layer in one event.
        const want = Math.max(-1.5, Math.min(1.5, amount / 70)) * 0.22 * M;
        // the pending value is re-applied with a 0.55 decay, so what finally lands
        // is 1/(1-0.55) = 2.22 times it
        const room = Math.max(0, IMP_CAP * M - impulse) / 2.22;
        const s = Math.sign(want) * Math.min(Math.abs(want), room);
        if (Math.abs(s) < 1e-5) return;
        impulse += Math.abs(s) * 2.22;
        pending.shear += s;
        frozen = false;
        agitate();
      },
      // In an isothermal gas piling up density *is* piling up pressure, so the
      // density bump does most of the work and the kick only aims it.
      blast(x, y, amp, radius) {
        const M = machNow();
        const kick = Math.min(1.8 * M, ((amp || 300) / 300) * 0.9 * M,
                              Math.max(0, SPLAT_BUDGET * M - pimp));
        pimp += kick;
        pending.blasts.push({
          x, y,
          dens: 1.0,
          kick: kick,
          radius: radius || 0.010
        });
        frozen = false;
        agitate();
      },
      setPreset(key) {
        const next = resolvePreset(key);
        state.preset = Object.assign({}, activePreset());
        state.target = next;
        state.blend = 0;
        frozen = false;
      },
      accentOf(key) { return resolvePreset(key).accent; },
      __bench: bench,
      // A column of primitives, for measuring a profile off the clock -- how a
      // deposited shear layer decays, for instance, which is not something a
      // screenshot can tell you.
      __col(ix) {
        const px = new Float32Array(grid.h * 4);
        gl.bindFramebuffer(gl.FRAMEBUFFER, U.read.fbo);
        gl.readPixels(ix | 0, 0, 1, grid.h, gl.RGBA, gl.FLOAT, px);
        const out = [];
        for (let j = 0; j < grid.h; j++) {
          const r = Math.max(px[j * 4], 1e-6);
          out.push([px[j * 4 + 1] / r, px[j * 4 + 2] / r]);
        }
        return out;
      },
      stats() {
        const pr = activePreset();
        return {
          gridW: grid.w, gridH: grid.h, nPart, nDraw, fps: state.fps, tier,
          tau: pr.tau, dpr,
          rms: state.machRms, max: state.machMax, div: state.div,
          mach: state.machRms, dens: state.dens,
          solver: (USE_HLLC ? 'HLLC' : 'HLL') + ' · ' +
            (SLOPE_TYPE === 0 ? 'piecewise constant'
              : RECON === 'ppm' ? 'PPM (CW84) · MonCen slopes'
              : 'MonCen slopes (slope_type 2) · trace2d') + ' · unsplit',
          mgLevels: null,
          // sim time and substeps per frame: between them they fix how fast the
          // picture moves, which is a design parameter here and not an accident
          ch: state.ctot, dtdx: dtdxNow(),
          time: state.time, sub: TIERS[tier][1], steps: state.steps, amp: state.amp,
          recon: (RECON === 'ppm' && agitated <= 0) ? 'ppm' : 'plm',
          timeScale: 1, targetFps: TARGET_FPS, cfl: CFL,
          drag: 'backward Euler'
        };
      }
    };
  }

  global.Field = {
    create(canvas) {
      try { return Field(canvas); }
      catch (e) {
        if (global.console) console.warn('[field-fv] ' + e.message);
        return null;
      }
    }
  };
})(window);

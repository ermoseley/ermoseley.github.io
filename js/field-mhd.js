/* ============================================================================
   field-mhd.js — the same background, magnetised.

   A third gas for the same page. field.js is incompressible (projection);
   field-fv.js is compressible hydrodynamics (the live site); this is compressible
   *magneto*hydrodynamics at plasma beta = 1, with charged dust.

   Everything above the gas -- the palette, the chapters, the dye, the grain size
   spectrum, the interaction surface site.js drives -- is unchanged. What differs
   is what the grains are falling through, and that they now feel it twice: once
   through drag, once through the Lorentz force.

   ---------------------------------------------------------------- the gas

     variables   : (rho, rho u, rho v, rho Y) and (Bx, By, psi), two RGBA32F
                   targets written in one pass through MRT. Isothermal,
                   P = rho cs^2, cs = 1, so every velocity is a Mach number and
                   every field strength is in units of rho^1/2 cs.
     div B       : Dedner et al. (2002) generalised Lagrange multiplier. B is
                   cell-centred and psi is an eighth variable that carries
                   divergence error out of the box at speed ch and is damped as it
                   goes. This is the alternative to the reference's constrained
                   transport, and it is the right one here: CT needs a staggered
                   field and a corner EMF solver, which is two more passes and a
                   second stencil, where GLM is two extra numbers in a texture that
                   the same 1D solver handles as part of its own flux.
     flux        : HLLD (Miyoshi & Kusano), the reference's hlld_mhd_fluxes reduced
                   to isothermal and to two dimensions -- five waves, so the
                   Alfven waves and the contact are all resolved and a shear layer
                   in the transverse field survives.
     reconstruct : PLM. slope_type = 2, the MonCen limiter from uslope, on all
                   seven primitives, plus a half-step predictor in the primitive
                   variables. No PPM on this page: a parabola fitted through a
                   fresh discontinuity is monotonised variable by variable, which
                   is fragile enough in hydrodynamics and worse when two of the
                   variables are a magnetic field whose divergence is supposed to
                   vanish.
     update      : unsplit and conservative, all four faces in one pass, plus the
                   one non-conservative term GLM needs (below).
     timestep    : dt = min(dt_default, dt_Courant), as on the live site. ctot is
                   now the sum over directions of |u_dim| + cf_dim with cf the fast
                   magnetosonic speed, so it is roughly twice the hydrodynamic
                   value at the same Mach number and the default step is set lower
                   to keep the fixed branch the binding one.

   Three places where this departs from the reference, all deliberate:

     1. Isothermal HLLD. hlld_mhd_fluxes is adiabatic and carries an energy
        equation; there is none here. P = rho cs^2 goes in wherever gamma*P/rho
        appeared, the etot star states and the v.B terms drop out entirely, and
        what is left is the same wave structure and the same algebra. The star
        region keeps a density on each side of the contact, exactly as the
        reference's does. Strictly, isothermal MHD has no entropy wave and the two
        should be equal; the same approximation is in the isothermal HLLC the live
        site runs, where it costs nothing visible and buys a resolved tangential
        discontinuity.

     2. GLM rather than CT, as above. The normal-field flux is not zero here
        (hlld_mhd_fluxes sets flux%Bx = 0 because CT owns that component): it is
        psi, and psi's own flux is ch^2 Bn. That pair is a linear 2x2 hyperbolic
        system with eigenvalues +-ch, decoupled from the rest, so it is solved
        exactly and added to the HLLD flux rather than approximated by it.

     3. The Powell source term -(div B) B in the momentum equation, which GLM needs
        and a CT scheme does not. Without it a residual divergence exerts a force
        along B that no physics put there. It is evaluated from the same
        face-averaged normal fields the psi flux uses, so the divergence the
        momentum equation feels is the one the cleaning is damping.

   ---------------------------------------------------------------- the dust

   Charged grains, charge-to-mass ratio 100 for the whole population. The update is
   mini-ramses's, from pm/move_fine.f90:

     - the force acts on the *drift* velocity w = v - u_gas, which is correct and
       not a shortcut: in ideal MHD the electric field in the lab frame is -u x B,
       so q(E + v x B) = q(w x B). A grain moving with the gas feels nothing.
     - compute_lorentz, the implicit Cayley form
         w' = Inverse[e - (dteff/2) m] . [e + (dteff/2) m] . w,   dteff = -dt*charge
       with m the cross-product matrix of B. That is an exact rotation about B for
       any step -- it preserves |w| to roundoff and cannot heat the population --
       and it is the reference's default (analytic_dust_force = .false.). The
       matrix is theirs, with Bz = 0 substituted; the sign of dteff is theirs too.
     - then drag. Here is the one exception the brief asked for: a single backward
       Euler stage after the Lorentz update, w <- w/(1 + dt/tau), rather than the
       reference's half-drag / Lorentz / half-drag splitting. One stage, first
       order, L-stable at any stopping time.

   The grains are 2.5-dimensional and the gas is not. B lies in the plane, so
   w x B for an in-plane drift points out of it: a charged grain cannot stay in the
   plane even though the flow it lives in is planar. So the population carries
   (vx, vy, vz) -- a second particle target -- and gyrates in the plane containing
   the field normal. The gas needs no such thing: with Bz = vz = 0 initially, their
   equations are homogeneous and they stay zero exactly.

   At charge 100 and |B| ~ 1.4 the gyrofrequency is 140 and the gyration is
   resolved at about thirty steps a turn (the reference's gyro_factor limits the
   same angle to 0.63 rad; this runs at 0.2). The gyroradius is a twelfth of a
   cell, so the grains are strongly magnetised: they slide along field lines and
   are held across them, and the dust draws the field rather than the flow.

   That anisotropy is also why the Kelvin-Helmholtz layer a scroll deposits behaves
   differently here. The mean field is vertical and the layer is horizontal, so the
   field threads it at right angles and its tension resists the roll-up -- the
   classical result is that a layer is stable to KH once delta_u < 2 v_A, and at
   beta = 1 that is delta_u < 2.8 cs. A scroll here makes Alfven waves and current
   sheets where the live site makes vortices. That is the physics of beta = 1 with a
   mean field, not a bug; BETA below is the knob if you want the rolls back.

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

  // Two attachments: the hydrodynamic state and the field. Explicit locations,
  // because ES 3.0 has no bindFragDataLocation.
  const HEAD2 = `#version 300 es
precision highp float;
precision highp sampler2D;
in vec2 vUv;
layout(location = 0) out vec4 outColor;
layout(location = 1) out vec4 outColor1;
`;

  // The primitive state, in two vec4s so one limiter serves both. The field slot
  // carries psi as its third component: it is limited and predicted like any other
  // primitive, and its fourth component is unused.
  const EOS = `
uniform float cs, cs2, smallr;

struct Prim { vec4 h; vec4 f; };      // h = (rho,u,v,Y), f = (Bx,By,psi,-)
struct Flux { vec4 h; float bn; float bt; float ps; };

vec4 toCons(vec4 q){ return vec4(q.x, q.x * q.y, q.x * q.z, q.x * q.w); }
`;

  // ---- Courant condition ---------------------------------------------------
  //
  // cmpdt, with the fast magnetosonic speed in place of cs. Two numbers come off
  // this reduction, not one:
  //
  //   x  ctot = sum over dim of (|u_dim| + cf_dim), which sets dt
  //   y  ch   = max over dim of (|u_dim| + cf_dim), the largest signal speed in
  //             the cell, which is what the GLM cleaning wave should travel at
  //
  // ch has to be the real maximum rather than ctot: it is a wave speed the
  // timestep must respect, and using the sum would either waste half the step or
  // -- worse, if the step were fixed regardless -- put the psi waves outside the
  // stability limit while every physical wave stayed inside it.
  const F_CMAX_STATE = HEAD + EOS + `
uniform sampler2D uSrc, uSrcB;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 8;
  vec2 m = vec2(0.0);
  for (int j = 0; j < 8; j++) {
    for (int i = 0; i < 8; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      vec4 U = texelFetch(uSrc, p, 0);
      vec4 B = texelFetch(uSrcB, p, 0);
      float ri = 1.0 / max(U.x, smallr);
      vec2  u  = U.yz * ri;
      // find_speed_fast: d2 = (b^2/rho + c^2)/2, cf^2 = d2 + sqrt(d2^2 - c^2 bn^2/rho)
      float d2 = 0.5 * (dot(B.xy, B.xy) * ri + cs2);
      float cx = sqrt(d2 + sqrt(max(d2 * d2 - cs2 * B.x * B.x * ri, 0.0)));
      float cy = sqrt(d2 + sqrt(max(d2 * d2 - cs2 * B.y * B.y * ri, 0.0)));
      float sx = abs(u.x) + cx, sy = abs(u.y) + cy;
      m = max(m, vec2(sx + sy, max(sx, sy)));
    }
  }
  outColor = vec4(m, 0.0, 1.0);
}`;

  const F_CMAX_DOWN = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 8;
  vec2 m = vec2(0.0);
  for (int j = 0; j < 8; j++) {
    for (int i = 0; i < 8; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      m = max(m, texelFetch(uSrc, p, 0).xy);
    }
  }
  outColor = vec4(m, 0.0, 1.0);
}`;

  // dt = min(dt_default, dt_adaptive), as on the live site: a fixed step chosen
  // inside the stability limit, with the Courant condition as a floor that engages
  // only when the flow gets faster than the state it was chosen for. The default is
  // lower here than in the hydrodynamic build because ctot is about twice as large
  // at the same Mach number -- at beta = 1 the fast speed is sqrt(cs^2 + vA^2) =
  // 1.73 across the field and vA = 1.41 along it, against cs = 1.
  const DTDX = `
uniform sampler2D uCmax;
uniform float cfl, smallc, dtdxMax;
vec2 waves(){ return max(texelFetch(uCmax, ivec2(0, 0), 0).xy, vec2(1e-20)); }
float stepDtDx(){ return min(dtdxMax, min(cfl / smallc, cfl / waves().x)); }
float chSpeed(){ return waves().y; }
`;

  // ---- the solver: one unsplit pass, two targets ---------------------------
  const F_GODUNOV = HEAD2 + EOS + DTDX + `
uniform sampler2D uU, uB;
uniform ivec2 size;
uniform float dyeDiss, dxCell, slopeType, psiDamp, powell, fric;

Prim prim(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  vec4 U = texelFetch(uU, p, 0);
  vec4 B = texelFetch(uB, p, 0);
  float r = max(U.x, smallr), ri = 1.0 / r;
  return Prim(vec4(r, U.y * ri, U.z * ri, U.w * ri), B);
}

// uslope, slope_type = 2: MonCen. dlft/drgt are the one-sided differences scaled
// by slope_type, dcen the central one; the limiter takes the smaller one-sided
// slope and zero where they disagree in sign.
vec4 uslope(vec4 qm, vec4 q0, vec4 qp){
  vec4 dlft = slopeType * (q0 - qm);
  vec4 drgt = slopeType * (qp - q0);
  vec4 dcen = 0.5 * (dlft + drgt) / max(slopeType, 1.0);
  vec4 dlim = min(abs(dlft), abs(drgt)) * step(vec4(0.0), dlft * drgt);
  return sign(dcen) * min(dlim, abs(dcen));
}

// The half-step predictor, trace2d's job, for isothermal GLM-MHD. Derived rather
// than transcribed: the reference's MHD trace predicts a staggered field through
// EMFs, which a cell-centred scheme has no use for. In primitive form,
//
//   D rho/Dt = -rho div u
//   D u  /Dt = (-cs^2 rho_x - By Jz) / rho        Jz = By_x - Bx_y
//   D v  /Dt = (-cs^2 rho_y + Bx Jz) / rho        so the force is (J x B)/rho
//   D Y  /Dt = 0
//   dBx/dt   = u_y By + u By_y - v_y Bx - v Bx_y - psi_x
//   dBy/dt   = -u_x By - u By_x + v_x Bx + v Bx_x - psi_y
//   dpsi/dt  = -ch^2 div B
//
// The magnetic pressure gradient and the field-line tension have been combined
// into J x B, which is exact and cheaper: the two cancel along B. hx/hy and gx/gy
// are half slopes, as drx/dry are in the reference, so the whole source is
// multiplied by dt/dx once.
void trace(inout Prim q, vec4 hx, vec4 hy, vec4 gx, vec4 gy, float dtdx, float ch2){
  float r = max(q.h.x, smallr), u = q.h.y, v = q.h.z, ri = 1.0 / r;
  float bx = q.f.x, by = q.f.y;
  float jz = gx.y - gy.x;
  vec4 sh, sf;
  sh.x = -u * hx.x - v * hy.x - (hx.y + hy.z) * r;
  sh.y = -u * hx.y - v * hy.y + (-cs2 * hx.x - by * jz) * ri;
  sh.z = -u * hx.z - v * hy.z + (-cs2 * hy.x + bx * jz) * ri;
  sh.w = -u * hx.w - v * hy.w;
  sf.x =  hy.y * by + u * gy.y - hy.z * bx - v * gy.x - gx.z;
  sf.y = -hx.y * by - u * gx.y + hx.z * bx + v * gx.x - gy.z;
  sf.z = -ch2 * (gx.x + gy.y);
  sf.w = 0.0;
  vec4 oh = q.h + sh * dtdx;
  // the reference falls back to the unpredicted value if the predictor undershoots
  // the density floor, rather than clamping to it
  oh.x = oh.x < smallr ? q.h.x : oh.x;
  q.h = oh;
  q.f = q.f + sf * dtdx;
}

// hlld_mhd_fluxes, isothermal and two-dimensional. Five waves: two fast, two
// Alfven, and the contact. The transverse pair (v_t, B_t) is what the Alfven waves
// rotate and what a two-wave solver would average away, which is the whole reason
// to pay for this one.
Flux hlld(Prim L, Prim R, float ch){
  float rl = max(L.h.x, smallr), rr = max(R.h.x, smallr);
  float ul = L.h.y, ur = R.h.y;
  float vl = L.h.z, vr = R.h.z;
  float bl = L.f.y, br = R.f.y;                 // transverse field
  float A  = 0.5 * (L.f.x + R.f.x);             // normal field: continuity enforced
  float A2 = A * A;
  float sgnm = A >= 0.0 ? 1.0 : -1.0;
  float ptl = rl * cs2 + 0.5 * (A2 + bl * bl);  // total pressure, gas + magnetic
  float ptr = rr * cs2 + 0.5 * (A2 + br * br);

  // fast magnetosonic speeds, with c^2 = cs^2 a uniform
  float ril = 1.0 / rl, rir = 1.0 / rr;
  float d2l = 0.5 * ((A2 + bl * bl) * ril + cs2);
  float d2r = 0.5 * ((A2 + br * br) * rir + cs2);
  float cfL = sqrt(d2l + sqrt(max(d2l * d2l - cs2 * A2 * ril, 0.0)));
  float cfR = sqrt(d2r + sqrt(max(d2r * d2r - cs2 * A2 * rir, 0.0)));
  float cm  = max(cfL, cfR);

  // HLL speeds and the star state, which follow from mass and normal-momentum
  // conservation across the whole fan and so do not involve the energy equation
  float SL = min(ul, ur) - cm;
  float SR = max(ul, ur) + cm;
  float rcl = rl * (ul - SL), rcr = rr * (SR - ur);
  float inv = 1.0 / (rcr + rcl);
  float ustar = (rcr * ur + rcl * ul + (ptl - ptr)) * inv;
  float pstar = (rcr * ptl + rcl * ptr + rcl * rcr * (ul - ur)) * inv;

  // left star region. The guard is the reference's: where estar vanishes the
  // Alfven wave has met the contact and the transverse state is simply continuous.
  float eps = max(1e-4 * A2, 1e-9);
  float rsl = max(rl * (SL - ul) / min(SL - ustar, -1e-8), smallr);
  float esl = rl * (SL - ul) * (SL - ustar) - A2;
  float ell = rl * (SL - ul) * (SL - ul)    - A2;
  float vsl = vl, bsl = bl;
  if (abs(esl) >= eps) {
    float ie = 1.0 / esl;
    vsl = vl - A * bl * (ustar - ul) * ie;
    bsl = bl * ell * ie;
  }
  float sql = sqrt(rsl), SAL = ustar - abs(A) / sql;

  // right star region
  float rsr = max(rr * (SR - ur) / max(SR - ustar, 1e-8), smallr);
  float esr = rr * (SR - ur) * (SR - ustar) - A2;
  float err = rr * (SR - ur) * (SR - ur)    - A2;
  float vsr = vr, bsr = br;
  if (abs(esr) >= eps) {
    float ie = 1.0 / esr;
    vsr = vr - A * br * (ustar - ur) * ie;
    bsr = br * err * ie;
  }
  float sqr = sqrt(rsr), SAR = ustar + abs(A) / sqr;

  // between the Alfven waves the transverse state is a single rotated one
  float den = 1.0 / (sql + sqr);
  float vss = (sql * vsl + sqr * vsr + sgnm * (bsr - bsl)) * den;
  float bss = (sql * bsr + sqr * bsl + sgnm * sql * sqr * (vsr - vsl)) * den;

  // sample at x/t = 0
  float ro, uo, vo, bo, pto;
  if (SL > 0.0)         { ro = rl;  uo = ul;    vo = vl;  bo = bl;  pto = ptl; }
  else if (SAL > 0.0)   { ro = rsl; uo = ustar; vo = vsl; bo = bsl; pto = pstar; }
  else if (ustar > 0.0) { ro = rsl; uo = ustar; vo = vss; bo = bss; pto = pstar; }
  else if (SAR > 0.0)   { ro = rsr; uo = ustar; vo = vss; bo = bss; pto = pstar; }
  else if (SR > 0.0)    { ro = rsr; uo = ustar; vo = vsr; bo = bsr; pto = pstar; }
  else                  { ro = rr;  uo = ur;    vo = vr;  bo = br;  pto = ptr; }

  float fd = ro * uo;
  Flux f;
  // scalar_flux: the dye rides on the mass flux, upwinded on its sign, which in
  // the star region is the side of the contact the interface lies on
  f.h  = vec4(fd, fd * uo + pto - A2, fd * vo - A * bo, fd * (fd >= 0.0 ? L.h.w : R.h.w));
  f.bt = bo * uo - A * vo;
  // GLM. (Bn, psi) is a linear system with eigenvalues +-ch and no coupling to
  // anything else, so its Riemann problem has an exact solution -- an upwind
  // average, which is what these two lines are -- and it is added to the flux
  // above rather than approximated by it. This is the term that makes flux%Bx
  // nonzero, where a constrained-transport scheme sets it to zero and lets the EMF
  // own that component.
  f.bn = 0.5 * (L.f.z + R.f.z) - 0.5 * ch * (R.f.x - L.f.x);
  f.ps = ch * (0.5 * ch * (L.f.x + R.f.x) - 0.5 * (R.f.z - L.f.z));
  return f;
}

vec4 swapv(vec4 a){ return vec4(a.x, a.z, a.y, a.w); }
// the y sweep, in the same frame the x sweep uses: normal velocity in .y, normal
// field in .x, exactly as the reference always puts the normal in velocity_x
Prim rot(Prim q){ return Prim(swapv(q.h), vec4(q.f.y, q.f.x, q.f.z, 0.0)); }

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  const ivec2 ex = ivec2(1, 0), ey = ivec2(0, 1);

  Prim q0 = prim(c);
  Prim qE = prim(c + ex),      qW = prim(c - ex);
  Prim qN = prim(c + ey),      qS = prim(c - ey);
  Prim qEE = prim(c + 2 * ex), qWW = prim(c - 2 * ex);
  Prim qNN = prim(c + 2 * ey), qSS = prim(c - 2 * ey);
  Prim qNE = prim(c + ex + ey), qSE = prim(c + ex - ey);
  Prim qNW = prim(c - ex + ey), qSW = prim(c - ex - ey);

  // half slopes on the five-cell cross, both directions and both halves of the
  // state: the predictor couples all of it
  vec4 h0x = 0.5 * uslope(qW.h, q0.h, qE.h),   h0y = 0.5 * uslope(qS.h, q0.h, qN.h);
  vec4 g0x = 0.5 * uslope(qW.f, q0.f, qE.f),   g0y = 0.5 * uslope(qS.f, q0.f, qN.f);
  vec4 hEx = 0.5 * uslope(q0.h, qE.h, qEE.h),  hEy = 0.5 * uslope(qSE.h, qE.h, qNE.h);
  vec4 gEx = 0.5 * uslope(q0.f, qE.f, qEE.f),  gEy = 0.5 * uslope(qSE.f, qE.f, qNE.f);
  vec4 hWx = 0.5 * uslope(qWW.h, qW.h, q0.h),  hWy = 0.5 * uslope(qSW.h, qW.h, qNW.h);
  vec4 gWx = 0.5 * uslope(qWW.f, qW.f, q0.f),  gWy = 0.5 * uslope(qSW.f, qW.f, qNW.f);
  vec4 hNx = 0.5 * uslope(qNW.h, qN.h, qNE.h), hNy = 0.5 * uslope(q0.h, qN.h, qNN.h);
  vec4 gNx = 0.5 * uslope(qNW.f, qN.f, qNE.f), gNy = 0.5 * uslope(q0.f, qN.f, qNN.f);
  vec4 hSx = 0.5 * uslope(qSW.h, qS.h, qSE.h), hSy = 0.5 * uslope(qSS.h, qS.h, q0.h);
  vec4 gSx = 0.5 * uslope(qSW.f, qS.f, qSE.f), gSy = 0.5 * uslope(qSS.f, qS.f, q0.f);

  float dtdx = stepDtDx();
  float ch = chSpeed(), ch2 = ch * ch;

  trace(q0, h0x, h0y, g0x, g0y, dtdx, ch2);
  trace(qE, hEx, hEy, gEx, gEy, dtdx, ch2);
  trace(qW, hWx, hWy, gWx, gWy, dtdx, ch2);
  trace(qN, hNx, hNy, gNx, gNy, dtdx, ch2);
  trace(qS, hSx, hSy, gSx, gSy, dtdx, ch2);

  Prim e0 = Prim(q0.h + h0x, q0.f + g0x), e1 = Prim(qE.h - hEx, qE.f - gEx);
  Prim w0 = Prim(qW.h + hWx, qW.f + gWx), w1 = Prim(q0.h - h0x, q0.f - g0x);
  Prim n0 = Prim(q0.h + h0y, q0.f + g0y), n1 = Prim(qN.h - hNy, qN.f - gNy);
  Prim s0 = Prim(qS.h + hSy, qS.f + gSy), s1 = Prim(q0.h - h0y, q0.f - g0y);

  Flux Fe = hlld(e0, e1, ch);
  Flux Fw = hlld(w0, w1, ch);
  Flux Fn = hlld(rot(n0), rot(n1), ch);
  Flux Fs = hlld(rot(s0), rot(s1), ch);

  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);

  vec4 dU = (Fw.h - Fe.h) + swapv(Fs.h - Fn.h);
  // the transverse flux of one component is the normal (psi) flux of the other
  vec2 dB = (vec2(Fw.bn, Fw.bt) - vec2(Fe.bn, Fe.bt))
          + (vec2(Fs.bt, Fs.bn) - vec2(Fn.bt, Fn.bn));
  float dpsi = (Fw.ps - Fe.ps) + (Fs.ps - Fn.ps);

  vec4 Un = U + dU * dtdx;
  // Powell. The divergence here is the one the psi flux sees -- the same
  // face-averaged normal fields -- so the force the momentum equation loses is
  // exactly the one the cleaning is working on.
  float divb = (0.5 * (e0.f.x + e1.f.x) - 0.5 * (w0.f.x + w1.f.x))
             + (0.5 * (n0.f.y + n1.f.y) - 0.5 * (s0.f.y + s1.f.y));
  Un.yz -= powell * divb * B.xy * dtdx;
  // A large-scale drag, backward Euler so it is stable at any step. Two dimensions
  // do not have three's energy cascade: enstrophy goes to the grid and energy goes
  // the other way, into box-scale coherent structures that nothing dissipates.
  // Measured, that is not a small effect here -- with the drive cut by a factor of a
  // hundred the rms Mach held at 0.85 for half a minute while the peaks decayed, so
  // the box was storing everything it had ever been given at the largest scale and
  // the servo had no plant to act on. This is the sink that gives it one, and it is
  // what a two-dimensional driven-turbulence run does for the same reason. The
  // timescale is longer than a box-scale turnover, so it sets the equilibrium
  // without touching the eddies.
  Un.yz = Un.yz / (1.0 + fric * dtdx * dxCell);
  Un.x = max(Un.x, smallr);
  Un.w = Un.w / (1.0 + dyeDiss * dtdx * dxCell);

  // and the parabolic half of Dedner's operator: psi is damped as it travels, so
  // the divergence error is not merely moved somewhere else in the box
  float psi = (B.z + dpsi * dtdx) * exp(-psiDamp * ch * dtdx);

  outColor  = Un;
  outColor1 = vec4(B.xy + dB * dtdx, psi, 0.0);
}`;

  // ---- third order in space: PPM for MHD ----------------------------------
  //
  // Colella & Woodward 1984 on all seven primitives, which is not in the reference
  // -- its slope_type 0 through 6 are every one of them piecewise linear -- so the
  // PLM program above, which *is* the reference's configuration, stays as the
  // alternative and as what runs while the box is being driven.
  //
  // Steps 1 to 3 are the hydrodynamic page's, componentwise, and the field rides
  // through them unchanged: limited MonCen slopes, the fourth-order interface value
  // a(i+1/2) = (a_i + a_i+1)/2 + (dq_i - dq_i+1)/6 built on them, and CW84 (1.10)
  // monotonisation -- flatten at an extremum, pull an overshooting edge back until
  // the parabola is monotone across the cell.
  //
  // Step 4 is where MHD differs, and getting it wrong would quietly undo the
  // divergence cleaning. The face state in PPM is the parabola averaged over the
  // slab the flow crosses in half a step, which is an advective time centring: it
  // carries the -u dq/dx term for you, at one speed for every variable. That is an
  // approximation for the hydrodynamic variables, all of which *are* advected along
  // the normal. It is simply wrong for two of the three field variables, because
  //
  //   dBn/dt  has no  -u dBn/dx  term at all -- the normal component is advected
  //           only transversely, which is the same fact constrained transport is
  //           built around
  //   dpsi/dt = -ch^2 div B  is not advected by anything
  //   dBt/dt  does have -u dBt/dx, so the transverse component is the one that
  //           belongs in the slab average
  //
  // Shifting Bn and psi by u dt/2 anyway would inject a spurious -u dBn/dx of
  // exactly the order of the terms being added, straight into the divergence the
  // cleaning is trying to damp. So the field takes the slab average in its
  // transverse component and the plain parabola edge in the other two -- one mix
  // against a constant mask, no extra cost -- and its source term keeps every term
  // the predictor has except the one the average now carries.
  const F_GODUNOV3 = HEAD2 + EOS + DTDX + `
uniform sampler2D uU, uB;
uniform ivec2 size;
uniform float dyeDiss, dxCell, slopeType, psiDamp, powell, fric;

Prim prim(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  vec4 U = texelFetch(uU, p, 0);
  vec4 B = texelFetch(uB, p, 0);
  float r = max(U.x, smallr), ri = 1.0 / r;
  return Prim(vec4(r, U.y * ri, U.z * ri, U.w * ri), B);
}

vec4 uslope(vec4 qm, vec4 q0, vec4 qp){
  vec4 dlft = slopeType * (q0 - qm);
  vec4 drgt = slopeType * (qp - q0);
  vec4 dcen = 0.5 * (dlft + drgt) / max(slopeType, 1.0);
  vec4 dlim = min(abs(dlft), abs(drgt)) * step(vec4(0.0), dlft * drgt);
  return sign(dcen) * min(dlim, abs(dcen));
}

// The monotonised parabola for one cell, from slopes computed by the caller.
// Adjacent cells' windows overlap by two, so the five slopes are computed once per
// sweep and passed in: that saves eight of eighteen uslope evaluations per variable
// group, and with seven variables uslope is the most expensive thing in the pass.
void ppmEdges(vec4 qm, vec4 q0, vec4 qp, vec4 dm, vec4 d0, vec4 dp, out vec4 aL, out vec4 aR){
  aL = 0.5 * (qm + q0) + (dm - d0) * (1.0 / 6.0);
  aR = 0.5 * (q0 + qp) + (d0 - dp) * (1.0 / 6.0);
  vec4 ext = step((aR - q0) * (q0 - aL), vec4(0.0));
  aL = mix(aL, q0, ext);
  aR = mix(aR, q0, ext);
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

// Which field components are advected along the sweep normal, and so belong in the
// slab average: the transverse one, and only that.
const vec4 ADVECTS = vec4(0.0, 1.0, 0.0, 0.0);

vec4 floorRho(vec4 f, vec4 q){ return f.x < smallr ? vec4(q.x, f.yzw) : f; }

// The predictor, minus whatever the parabola average already carries. hx is the
// parabola's effective half slope along the normal, hy the transverse one.
vec4 srcH(vec4 h, vec4 f, vec4 hx, vec4 hy, vec4 gx, vec4 gy){
  float r = max(h.x, smallr), v = h.z, ri = 1.0 / r;
  float jz = gx.y - gy.x;
  vec4 s;
  s.x = -v * hy.x - (hx.y + hy.z) * r;
  s.y = -v * hy.y + (-cs2 * hx.x - f.y * jz) * ri;
  s.z = -v * hy.z + (-cs2 * hy.x + f.x * jz) * ri;
  s.w = -v * hy.w;
  return s;
}
vec4 srcF(vec4 h, vec4 f, vec4 hx, vec4 hy, vec4 gx, vec4 gy, float ch2){
  float u = h.y, v = h.z;
  vec4 s;
  s.x =  hy.y * f.y + u * gy.y - hy.z * f.x - v * gy.x - gx.z;
  s.y = -hx.y * f.y             + hx.z * f.x + v * gx.x - gy.z;   // -u*gx.y dropped
  s.z = -ch2 * (gx.x + gy.y);
  s.w = 0.0;
  return s;
}

Flux hlld(Prim L, Prim R, float ch){
  float rl = max(L.h.x, smallr), rr = max(R.h.x, smallr);
  float ul = L.h.y, ur = R.h.y;
  float vl = L.h.z, vr = R.h.z;
  float bl = L.f.y, br = R.f.y;
  float A  = 0.5 * (L.f.x + R.f.x);
  float A2 = A * A;
  float sgnm = A >= 0.0 ? 1.0 : -1.0;
  float ptl = rl * cs2 + 0.5 * (A2 + bl * bl);
  float ptr = rr * cs2 + 0.5 * (A2 + br * br);
  float ril = 1.0 / rl, rir = 1.0 / rr;
  float d2l = 0.5 * ((A2 + bl * bl) * ril + cs2);
  float d2r = 0.5 * ((A2 + br * br) * rir + cs2);
  float cfL = sqrt(d2l + sqrt(max(d2l * d2l - cs2 * A2 * ril, 0.0)));
  float cfR = sqrt(d2r + sqrt(max(d2r * d2r - cs2 * A2 * rir, 0.0)));
  float cm  = max(cfL, cfR);
  float SL = min(ul, ur) - cm;
  float SR = max(ul, ur) + cm;
  float rcl = rl * (ul - SL), rcr = rr * (SR - ur);
  float inv = 1.0 / (rcr + rcl);
  float ustar = (rcr * ur + rcl * ul + (ptl - ptr)) * inv;
  float pstar = (rcr * ptl + rcl * ptr + rcl * rcr * (ul - ur)) * inv;
  float eps = max(1e-4 * A2, 1e-9);
  float rsl = max(rl * (SL - ul) / min(SL - ustar, -1e-8), smallr);
  float esl = rl * (SL - ul) * (SL - ustar) - A2;
  float ell = rl * (SL - ul) * (SL - ul)    - A2;
  float vsl = vl, bsl = bl;
  if (abs(esl) >= eps) {
    float ie = 1.0 / esl;
    vsl = vl - A * bl * (ustar - ul) * ie;
    bsl = bl * ell * ie;
  }
  float sql = sqrt(rsl), SAL = ustar - abs(A) / sql;
  float rsr = max(rr * (SR - ur) / max(SR - ustar, 1e-8), smallr);
  float esr = rr * (SR - ur) * (SR - ustar) - A2;
  float err = rr * (SR - ur) * (SR - ur)    - A2;
  float vsr = vr, bsr = br;
  if (abs(esr) >= eps) {
    float ie = 1.0 / esr;
    vsr = vr - A * br * (ustar - ur) * ie;
    bsr = br * err * ie;
  }
  float sqr = sqrt(rsr), SAR = ustar + abs(A) / sqr;
  float den = 1.0 / (sql + sqr);
  float vss = (sql * vsl + sqr * vsr + sgnm * (bsr - bsl)) * den;
  float bss = (sql * bsr + sqr * bsl + sgnm * sql * sqr * (vsr - vsl)) * den;
  float ro, uo, vo, bo, pto;
  if (SL > 0.0)         { ro = rl;  uo = ul;    vo = vl;  bo = bl;  pto = ptl; }
  else if (SAL > 0.0)   { ro = rsl; uo = ustar; vo = vsl; bo = bsl; pto = pstar; }
  else if (ustar > 0.0) { ro = rsl; uo = ustar; vo = vss; bo = bss; pto = pstar; }
  else if (SAR > 0.0)   { ro = rsr; uo = ustar; vo = vss; bo = bss; pto = pstar; }
  else if (SR > 0.0)    { ro = rsr; uo = ustar; vo = vsr; bo = bsr; pto = pstar; }
  else                  { ro = rr;  uo = ur;    vo = vr;  bo = br;  pto = ptr; }
  float fd = ro * uo;
  Flux fx;
  fx.h  = vec4(fd, fd * uo + pto - A2, fd * vo - A * bo, fd * (fd >= 0.0 ? L.h.w : R.h.w));
  fx.bt = bo * uo - A * vo;
  fx.bn = 0.5 * (L.f.z + R.f.z) - 0.5 * ch * (R.f.x - L.f.x);
  fx.ps = ch * (0.5 * ch * (L.f.x + R.f.x) - 0.5 * (R.f.z - L.f.z));
  return fx;
}

vec4 swapv(vec4 a){ return vec4(a.x, a.z, a.y, a.w); }
Prim rot(Prim q){ return Prim(swapv(q.h), vec4(q.f.y, q.f.x, q.f.z, 0.0)); }

// One sweep, in the frame where the normal velocity is in .y and the normal field
// in .x. a3..b3 run along the normal; the six t* are the transverse neighbours of
// the three cells that own the two faces. An* come back so the caller can build the
// same face-averaged divergence the psi flux sees.
void sweep(Prim am3, Prim am2, Prim am1, Prim a0, Prim ap1, Prim ap2, Prim ap3,
           Prim tmm, Prim tpm, Prim tm0, Prim tp0, Prim tmp, Prim tpp,
           float dtdx, float ch, out Flux Flo, out Flux Fhi,
           out float Anlo, out float Anhi){
  vec4 h2m = uslope(am3.h, am2.h, am1.h), f2m = uslope(am3.f, am2.f, am1.f);
  vec4 h1m = uslope(am2.h, am1.h, a0.h),  f1m = uslope(am2.f, am1.f, a0.f);
  vec4 h00 = uslope(am1.h, a0.h,  ap1.h), f00 = uslope(am1.f, a0.f,  ap1.f);
  vec4 h1p = uslope(a0.h,  ap1.h, ap2.h), f1p = uslope(a0.f,  ap1.f, ap2.f);
  vec4 h2p = uslope(ap1.h, ap2.h, ap3.h), f2p = uslope(ap1.f, ap2.f, ap3.f);

  vec4 hlL, hlR, hcL, hcR, hrL, hrR, flL, flR, fcL, fcR, frL, frR;
  ppmEdges(am2.h, am1.h, a0.h,  h2m, h1m, h00, hlL, hlR);
  ppmEdges(am1.h, a0.h,  ap1.h, h1m, h00, h1p, hcL, hcR);
  ppmEdges(a0.h,  ap1.h, ap2.h, h00, h1p, h2p, hrL, hrR);
  ppmEdges(am2.f, am1.f, a0.f,  f2m, f1m, f00, flL, flR);
  ppmEdges(am1.f, a0.f,  ap1.f, f1m, f00, f1p, fcL, fcR);
  ppmEdges(a0.f,  ap1.f, ap2.f, f00, f1p, f2p, frL, frR);

  vec4 hyL = 0.5 * uslope(tmm.h, am1.h, tpm.h), gyL = 0.5 * uslope(tmm.f, am1.f, tpm.f);
  vec4 hy0 = 0.5 * uslope(tm0.h, a0.h,  tp0.h), gy0 = 0.5 * uslope(tm0.f, a0.f,  tp0.f);
  vec4 hyR = 0.5 * uslope(tmp.h, ap1.h, tpp.h), gyR = 0.5 * uslope(tmp.f, ap1.f, tpp.f);

  float sL = clamp(abs(am1.h.y) * dtdx, 0.0, 1.0);
  float s0 = clamp(abs(a0.h.y)  * dtdx, 0.0, 1.0);
  float sR = clamp(abs(ap1.h.y) * dtdx, 0.0, 1.0);

  float ch2 = ch * ch;
  vec4 hxL = 0.5 * (hlR - hlL), gxL = 0.5 * (flR - flL);
  vec4 hx0 = 0.5 * (hcR - hcL), gx0 = 0.5 * (fcR - fcL);
  vec4 hxR = 0.5 * (hrR - hrL), gxR = 0.5 * (frR - frL);

  vec4 ehL = srcH(am1.h, am1.f, hxL, hyL, gxL, gyL) * dtdx;
  vec4 ehC = srcH(a0.h,  a0.f,  hx0, hy0, gx0, gy0) * dtdx;
  vec4 ehR = srcH(ap1.h, ap1.f, hxR, hyR, gxR, gyR) * dtdx;
  vec4 efL = srcF(am1.h, am1.f, hxL, hyL, gxL, gyL, ch2) * dtdx;
  vec4 efC = srcF(a0.h,  a0.f,  hx0, hy0, gx0, gy0, ch2) * dtdx;
  vec4 efR = srcF(ap1.h, ap1.f, hxR, hyR, gxR, gyR, ch2) * dtdx;

  // The reference floors its reconstructed states, and a parabola through a fresh
  // discontinuity can undershoot: a face density at the floor beside a large
  // velocity makes rcl tiny, ustar enormous and the flux garbage.
  Prim lo0 = Prim(floorRho(ppmFace(am1.h, hlL, hlR, sL, 1.0) + ehL, am1.h),
                  mix(flR, ppmFace(am1.f, flL, flR, sL, 1.0), ADVECTS) + efL);
  Prim lo1 = Prim(floorRho(ppmFace(a0.h,  hcL, hcR, s0, 0.0) + ehC, a0.h),
                  mix(fcL, ppmFace(a0.f,  fcL, fcR, s0, 0.0), ADVECTS) + efC);
  Prim hi0 = Prim(floorRho(ppmFace(a0.h,  hcL, hcR, s0, 1.0) + ehC, a0.h),
                  mix(fcR, ppmFace(a0.f,  fcL, fcR, s0, 1.0), ADVECTS) + efC);
  Prim hi1 = Prim(floorRho(ppmFace(ap1.h, hrL, hrR, sR, 0.0) + ehR, ap1.h),
                  mix(frL, ppmFace(ap1.f, frL, frR, sR, 0.0), ADVECTS) + efR);

  Anlo = 0.5 * (lo0.f.x + lo1.f.x);
  Anhi = 0.5 * (hi0.f.x + hi1.f.x);
  Flo = hlld(lo0, lo1, ch);
  Fhi = hlld(hi0, hi1, ch);
}

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  const ivec2 ex = ivec2(1, 0), ey = ivec2(0, 1);

  Prim q0 = prim(c);
  Prim x1p = prim(c + ex),     x1m = prim(c - ex);
  Prim x2p = prim(c + 2 * ex), x2m = prim(c - 2 * ex);
  Prim x3p = prim(c + 3 * ex), x3m = prim(c - 3 * ex);
  Prim y1p = prim(c + ey),     y1m = prim(c - ey);
  Prim y2p = prim(c + 2 * ey), y2m = prim(c - 2 * ey);
  Prim y3p = prim(c + 3 * ey), y3m = prim(c - 3 * ey);
  Prim dNE = prim(c + ex + ey), dSE = prim(c + ex - ey);
  Prim dNW = prim(c - ex + ey), dSW = prim(c - ex - ey);

  float dtdx = stepDtDx();
  float ch = chSpeed();

  Flux Fw, Fe; float Aw, Ae;
  sweep(x3m, x2m, x1m, q0, x1p, x2p, x3p,
        dSW, dNW, y1m, y1p, dSE, dNE, dtdx, ch, Fw, Fe, Aw, Ae);

  Flux Fs, Fn; float As, An;
  sweep(rot(y3m), rot(y2m), rot(y1m), rot(q0), rot(y1p), rot(y2p), rot(y3p),
        rot(dSW), rot(dSE), rot(x1m), rot(x1p), rot(dNW), rot(dNE),
        dtdx, ch, Fs, Fn, As, An);

  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);

  vec4 dU = (Fw.h - Fe.h) + swapv(Fs.h - Fn.h);
  vec2 dB = (vec2(Fw.bn, Fw.bt) - vec2(Fe.bn, Fe.bt))
          + (vec2(Fs.bt, Fs.bn) - vec2(Fn.bt, Fn.bn));
  float dpsi = (Fw.ps - Fe.ps) + (Fs.ps - Fn.ps);

  vec4 Un = U + dU * dtdx;
  Un.yz -= powell * ((Ae - Aw) + (An - As)) * B.xy * dtdx;
  Un.yz = Un.yz / (1.0 + fric * dtdx * dxCell);
  Un.x = max(Un.x, smallr);
  Un.w = Un.w / (1.0 + dyeDiss * dtdx * dxCell);
  float psi = (B.z + dpsi * dtdx) * exp(-psiDamp * ch * dtdx);

  outColor  = Un;
  outColor1 = vec4(B.xy + dB * dtdx, psi, 0.0);
}`;

  // ---- forcing ------------------------------------------------------------
  //
  // Every pass that writes the state writes both targets, so the field rides
  // through untouched rather than going stale behind a ping-pong swap. It costs one
  // extra texel read and write on a grid of fifty thousand, which is nothing beside
  // the solver.
  const PASS = `
uniform sampler2D uB;
vec4 passField(){ return texelFetch(uB, ivec2(gl_FragCoord.xy), 0); }
`;

  const F_STIR = HEAD2 + PASS + DTDX + `
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
  vec2 m = U.yz + r * a * amp * dt;
  float w = exp(-dot((vUv - c0) * vec2(aspect, 1.0), (vUv - c0) * vec2(aspect, 1.0)) / dyeRadius)
          + exp(-dot((vUv - c1) * vec2(aspect, 1.0), (vUv - c1) * vec2(aspect, 1.0)) / dyeRadius);
  outColor  = vec4(U.x, m, U.w + r * dye.x * w * dt);
  outColor1 = passField();
}`;

  const F_SPLAT = HEAD2 + PASS + `
uniform sampler2D uU;
uniform vec2  point, delta;
uniform float aspect, radius, ink, smallr;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float w = exp(-dot(p, p) / radius);
  float r = max(U.x, smallr);
  outColor  = vec4(U.x, U.yz + delta * w * r, U.w + r * ink * w);
  outColor1 = passField();
}`;

  // A tanh shear layer with a sinusoidal seed, as on the live site. Whether it
  // rolls up is now a race with the field as well as with the ambient strain --
  // see the note at the top of this file.
  const F_SHEAR = HEAD2 + PASS + `
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
  outColor  = vec4(U.x, U.yz + r * du, U.w + r * band.x * bw);
  outColor1 = passField();
}`;

  const F_BLAST = HEAD2 + PASS + `
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
  vec2  m  = U.yz * (r1 / r0) + kick * w * (p / rr) * r0;
  outColor  = vec4(r1, m, U.w * (r1 / r0));
  outColor1 = passField();
}`;

  // The initial field is uniform and vertical, which is exactly divergence-free:
  // the cleaning diagnostic then measures what the scheme does to a clean field
  // rather than what it inherited. |B| = sqrt(2 P / beta), and with rho = cs = 1
  // that is sqrt(2/beta) -- 1.414 at beta = 1, so vA = 1.41 against cs = 1.
  const F_INIT = HEAD2 + `
uniform float b0;
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
  float a = sin(TAU * uv.x * 2.0), b = cos(TAU * uv.y);
  float c = sin(TAU * (uv.x + uv.y)), d = cos(TAU * (uv.x * 2.0 - uv.y));
  vec2  u = 0.8 * vec2(a * b + 0.5 * d, -c + 0.4 * a);
  float Y = 0.04 + 0.03 * sin(TAU * (uv.x * 3.0 + uv.y)) * cos(TAU * (uv.y * 2.0 - uv.x));
  outColor  = vec4(rho, rho * u, rho * Y);
  outColor1 = vec4(0.0, b0, 0.0, 0.0);
}`;

  // What the dust samples: velocity and field in one filterable texture, so a
  // grain gathers both in four bilinear taps instead of eight. RGBA16F, because
  // asking for LINEAR on RGBA32F makes the texture incomplete and every read --
  // texelFetch included -- silently returns zero.
  const F_FLOW = HEAD + `
uniform sampler2D uU, uB;
uniform float smallr;
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);
  outColor = vec4(U.yz / max(U.x, smallr), B.xy);
}`;

  // Everything the picture needs, at grid resolution, so the composite can afford
  // a proper reconstruction: the dye, and the convergence of the velocity field.
  const F_DISP = HEAD + `
uniform sampler2D uU;
uniform ivec2 size;
uniform float smallr;
vec4 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0);
}
vec3 q(ivec2 c){ vec4 U = cell(c); float ri = 1.0 / max(U.x, smallr); return vec3(U.w * ri, U.yz * ri); }
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec3 q0 = q(c);
  vec3 qe = q(c + ivec2(1, 0)), qw = q(c - ivec2(1, 0));
  vec3 qn = q(c + ivec2(0, 1)), qs = q(c - ivec2(0, 1));
  float div = 0.5 * ((qe.y - qw.y) + (qn.z - qs.z));
  outColor = vec4(q0.x, max(-div, 0.0), 0.0, 1.0);
}`;


  // ---- the display field, smoothed ----------------------------------------
  //
  // Why this pass exists. The convergence channel is a centred difference of the
  // velocity, so a captured front is one cell wide and effectively all of its power
  // sits at the grid's Nyquist frequency. No interpolant can fix that: Catmull-Rom
  // is C1 and passes through the samples, which is exactly right for the dye, but
  // fed a one-cell spike it has no sub-cell information to work with and reproduces
  // the lattice -- lozenges along cell boundaries, and ringing from its negative
  // lobes. Stretched over eight screen pixels per cell that reads as grid artifacts
  // on every shock, which is what it was.
  //
  // So the smoothing happens before the reconstruction, not in it: a 5x5 binomial on
  // the convergence channel at grid resolution, which turns the spike into a
  // resolved two-and-a-half-cell feature the interpolant can actually represent. It
  // costs one pass over seventeen thousand pixels, against doing anything at all at
  // display resolution, where there are seventy times as many. The dye is passed
  // through untouched -- it is a smooth field already and blurring it would cost the
  // wisps their contrast.
  const F_DSMOOTH = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 size;
uniform float sigma;
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  float s2 = 2.0 * max(sigma * sigma, 1e-6);
  float acc = 0.0, wsum = 0.0;
  for (int j = -2; j <= 2; j++) {
    for (int i = -2; i <= 2; i++) {
      ivec2 p = c + ivec2(i, j);
      p = ivec2((p.x + size.x) % size.x, (p.y + size.y) % size.y);
      float w = exp(-float(i * i + j * j) / s2);
      acc += w * texelFetch(uSrc, p, 0).y;
      wsum += w;
    }
  }
  outColor = vec4(texelFetch(uSrc, c, 0).x, acc / wsum, 0.0, 1.0);
}`;

  // ---- reconstruction -----------------------------------------------------

  // Catmull-Rom for the picture: it interpolates, so a front stays where the
  // solver put it, and it is C1 across cell boundaries, which is what kills the
  // facets bilinear leaves when a 128-cell grid is stretched over 800 pixels.
  // Sixteen point samples in nine hardware taps.
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

  // A cubic B-spline for the display field, alongside the Catmull-Rom above --
  // because the two channels want opposite things.
  //
  // Catmull-Rom interpolates, which is what the dye wants: it passes through the
  // samples, so a wisp keeps its amplitude. But it is built from a kernel with
  // negative lobes, and negative lobes on a one-cell feature ring: an undershoot
  // ridge on each side, aligned with the axes because the kernel is separable.
  // That is what the facets on the fronts actually were -- not the interpolation
  // being too coarse, but it overshooting a feature it could not represent.
  //
  // The B-spline is C2 and strictly positive, so it cannot ring at all. It
  // approximates instead of interpolating, which costs a front about half a cell of
  // width, and that is exactly the trade to take here: nobody is reading a
  // convergence value off the screen. Four bilinear taps against nine.
  const BSPLINE4 = `
vec4 texBS(sampler2D tex, vec2 uv, vec2 size){
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
  return g0.y * (g0.x * texture(tex, vec2(h0.x, h0.y)) +
                 g1.x * texture(tex, vec2(h1.x, h0.y)))
       + g1.y * (g0.x * texture(tex, vec2(h0.x, h1.y)) +
                 g1.x * texture(tex, vec2(h1.x, h1.y)));
}
`;

  // A cubic B-spline for the grains: the shape function a particle-in-cell code
  // gathers with, and for the same reason -- it is C2, so a grain crossing a cell
  // boundary feels no kink in its acceleration and the population cannot print the
  // grid onto its own clustering. Four bilinear taps, and here they carry the
  // field as well as the velocity.
  const BSPLINE = `
vec4 gather4(sampler2D tex, vec2 uv, vec2 size){
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
  return g0.y * (g0.x * texture(tex, vec2(h0.x, h0.y)) +
                 g1.x * texture(tex, vec2(h1.x, h0.y)))
       + g1.y * (g0.x * texture(tex, vec2(h0.x, h1.y)) +
                 g1.x * texture(tex, vec2(h1.x, h1.y)));
}
`;

  // ---- dust ---------------------------------------------------------------

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
float aRel(float s){ return exp((s - 1.0) * 2.0 * 1.6118); }
`;

  // compute_lorentz, from mini-ramses pm/move_fine.f90, with Bz = 0 substituted
  // and the sign of dteff as the reference carries it. The matrix is the Cayley
  // transform of the cross-product matrix,
  //
  //   w' = Inverse[e - (dteff/2) m] . [e + (dteff/2) m] . w
  //
  // which is an exact rotation about B through 2 atan(|B| dteff / 2) -- so |w| is
  // preserved for any step and the drift along B is untouched, both of which are
  // visible in the algebra below: the .y row divides by det where the diagonal is
  // det. That is what makes this safe at charge 100, where an explicit push would
  // need the gyration substepped.
  const LORENTZ = `
vec3 lorentz(vec3 w, vec2 B, float dt, float q){
  float d  = -dt * q;                       // dteff
  float d2 = d * d;
  float bx = B.x, by = B.y;
  float b2 = bx * bx + by * by;
  float det = 1.0 + 0.25 * b2 * d2;
  float k = 0.25 * (bx * bx - by * by) * d2;
  float s = 0.5 * d2 * bx * by;
  vec3 o;
  o.x = (1.0 + k) * w.x + s * w.y           - d * by * w.z;
  o.y = s * w.x         + (1.0 - k) * w.y   + d * bx * w.z;
  o.z = d * by * w.x    - d * bx * w.y      + (1.0 - 0.25 * b2 * d2) * w.z;
  return o / det;
}
`;

  // The grain push: gather, rotate, drag, move. Two targets, because a charged
  // grain in a planar field does not stay planar.
  const F_PART = HEAD2 + SPECTRUM + BSPLINE + LORENTZ + `
uniform sampler2D uPart, uPartZ, uFlow;
uniform vec2  gridSize;
uniform float dragA, dt, aspect, brown, charge, stream;
uniform uint  uFrame;

void main(){
  ivec2 id  = ivec2(gl_FragCoord.xy);
  vec4  P   = texture(uPart, vUv);
  vec4  Z   = texture(uPartZ, vUv);
  vec2  pos = P.xy;
  vec3  v   = vec3(P.zw, Z.x);

  vec4 g = gather4(uFlow, pos, gridSize);
  vec3 u = vec3(g.xy, 0.0);
  vec2 B = g.zw;

  // In ideal MHD the lab-frame electric field is -u x B, so the whole Lorentz
  // force on a grain is q (v - u) x B: it acts on the drift and nothing else.
  vec3 w = v - u;
  w = lorentz(w, B, dt, charge);

  // then one backward Euler drag stage, on a per-grain stopping time. tau varies
  // across the population but not in time, so the median grain's dt/tau is the
  // only scalar the CPU supplies and a grain rescales it by its place in the size
  // spectrum.
  float a = dragA / tauMul(stRank(id));
  w = w / (1.0 + a);

  if (brown > 0.0) {
    w.xy += brown * (vec2(hash1(uvec2(id), uFrame + 11u),
                          hash1(uvec2(id), uFrame + 523u)) - 0.5);
  }
  v = u + w;

  // bulk streaming, along the local field now that there is one to follow
  vec2 sdir = stream * B / max(length(B), 1e-4);
  vec3 drift = w + vec3(sdir, 0.0);

  pos = fract(pos + vec2((v.x + sdir.x) / aspect, v.y + sdir.y) * dt);
  outColor  = vec4(pos, v.xy);
  outColor1 = vec4(v.z, length(drift), 0.0, 1.0);
}`;

  // The drift scale, measured rather than predicted. On the live site the
  // lighting was normalised by tau * M^2 -- the drift a grain of stopping time tau
  // reaches in a flow of rms Mach M -- and that estimate is simply wrong once the
  // grains are magnetised: the cross-field drift is suppressed by 1/(1 + (omega
  // tau)^2), which at charge 100 is a factor of two hundred, so every grain would
  // sit at the floor of the ramp. The population's own mean drift is reduced on the
  // GPU and read straight out of the vertex shader, so nothing has to be guessed
  // and nothing has to come back to the CPU.
  const F_PDRIFT = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
uniform float first;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 8;
  float sum = 0.0, mx = 0.0;
  for (int j = 0; j < 8; j++) {
    for (int i = 0; i < 8; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      vec4 v = texelFetch(uSrc, p, 0);
      // The first pass reads the particle state, where the drift is .y; later
      // passes read their own output, where the running sum is .x. The maximum
      // lives in .y in both layouts, so it needs no branch.
      sum += first > 0.5 ? v.y : v.x;
      mx = max(mx, v.y);
    }
  }
  outColor = vec4(sum, mx, 0.0, 1.0);
}`;

  const V_PARTICLE = `#version 300 es
precision highp float;
uniform sampler2D uPart, uPartZ, uDrift;
uniform int   uPW;
uniform float uPointSize, uDriftRef, uInvN, uAlpha, uVignette, uAspect;
out float vBright;
out float vDrift;
` + SPECTRUM + `
void main(){
  int   id  = gl_VertexID;
  ivec2 tc  = ivec2(id % uPW, id / uPW);
  vec2  pos = texelFetch(uPart, tc, 0).xy;
  float w   = texelFetch(uPartZ, tc, 0).y;
  float st  = stRank(tc);

  // Drift, not speed: a grain swept along with the gas is doing nothing
  // interesting however fast it is going. The drift is the one the grain actually
  // felt, computed in the push from the gather the drag used, and it includes the
  // out-of-plane component the field put there.
  float mean = texelFetch(uDrift, ivec2(0, 0), 0).x * uInvN;
  vDrift = clamp(w / max(uDriftRef * mean, 1e-5), 0.0, 1.0);

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

  // ---- the field, drawn ---------------------------------------------------
  //
  // Line integral convolution: walk a few cells along the field line through both
  // ends of each pixel, average a static noise field over the walk, and what comes
  // back is a texture whose grain lies along B. It is the standard way to draw a
  // vector field without drawing arrows on it, and it is the right one here for a
  // reason that is not only aesthetic: the dust at charge 100 is already stuck to
  // field lines, so the streaks it produces run the same way the grains do. The
  // field reinforces the thing the dust is telling you instead of annotating it.
  //
  // Three choices that keep it inside the palette rather than on top of it:
  //
  //   - the contribution is *signed* around zero, so the field grains the medium
  //     rather than lighting it, and the mean brightness of the frame does not move
  //   - it is drawn in the chapter's own accent, like the dye and the fronts. No
  //     hue appears here that appears nowhere else
  //   - it is scaled by |B| / B0, so a stretched or compressed field shows itself
  //     and an undisturbed one stays quiet
  //
  // The kernel travels along the walk, which makes bright packets run along the
  // field lines at a fixed rate. That is not decoration either: it is what an Alfven
  // wave does, and it is the only motion in the picture that is not advective.
  //
  // Run at half the display resolution into one RG16F, because the pattern is smooth
  // by construction and eleven taps per pixel at full resolution is not free.
  const F_LIC = HEAD + `
uniform sampler2D uFlow, uNoise;
uniform float aspect, time, b0, stepUv, noiseFreq, waveLen, waveSpeed;

void main(){
  vec2 uv = vUv;
  vec4 g0 = texture(uFlow, uv);
  float bmag = length(g0.zw);
  if (bmag < 1e-5) { outColor = vec4(0.0, 0.0, 0.0, 1.0); return; }

  // Both the noise grain and the walk step are set in pixels of *this* target, not
  // in box units, and that is not a detail -- it is the whole thing working or not.
  // Written in box units first, the noise repeated a hundred and ten times across
  // the box, which put its texels thirty-five times below the pixel Nyquist: every
  // sample along the walk drew an independent random number, the average went to a
  // flat half everywhere, and there was no streak to see. The step has the same
  // requirement from the other side -- step further than the noise is correlated and
  // consecutive samples are again independent. So one noise texel is a few target
  // pixels, and the step is the same few pixels.
  const int N = 8;                       // 2N+1 samples along the field line
  float h = stepUv;                      // step, in uv, sized from the target
  float acc = 0.0, wsum = 0.0;

  // Both directions from the pixel, following the local field as it turns. The
  // walk re-reads B at every step, so a curved field line stays followed rather
  // than being approximated by its tangent at the centre.
  for (int side = 0; side < 2; side++) {
    vec2 p = uv;
    float sgn = side == 0 ? 1.0 : -1.0;
    vec2 b = g0.zw;
    for (int i = 0; i <= N; i++) {
      if (i > 0) {
        b = texture(uFlow, p).zw;
        float m = length(b);
        if (m < 1e-5) break;
        b /= m;
        // physical -> uv: x spans aspect box widths, y spans one
        p += sgn * h * vec2(b.x / aspect, b.y);
      }
      // arc length from the centre, and a kernel that travels along it
      float s = float(i) * h;
      float w = 0.55 + 0.45 * cos(6.2831853 * (s / waveLen - time * waveSpeed));
      if (side == 1 && i == 0) continue;  // the centre sample belongs to one side only
      float n = texture(uNoise, vec2(p.x * aspect, p.y) * noiseFreq).x;
      acc += n * w;
      wsum += w;
    }
  }
  // Averaging eleven samples of white noise leaves a standard deviation of
  // 0.29/sqrt(11) ~ 0.09, not 0.5, so the signed value is normalised here rather
  // than leaving the gain below to absorb a factor of three it cannot be read as.
  float lic = wsum > 0.0 ? acc / wsum : 0.5;
  outColor = vec4((lic - 0.5) * 3.3, bmag / b0, 0.0, 1.0);
}`;

  // ---- composite ----------------------------------------------------------
  //
  // Otherwise identical to the live site's, deliberately: this page is a change of
  // physics, not a change of picture.
  const F_COMPOSITE = HEAD + CATROM + BSPLINE4 + `
uniform sampler2D uDisp, uLic;
uniform float uFieldGain;
uniform vec2  gridSize;
uniform vec2  uRes;
uniform vec3  uBg, uTint, uAccent;
uniform float uTime, uGrain, uDyeGain, uVignette, uShockGain, uShockLift, uShockPow;

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
  // One reconstruction per channel, because they want opposite things: the dye
  // interpolated so a wisp keeps its amplitude, the convergence approximated so a
  // front cannot ring. See BSPLINE4.
  vec2 d2 = vec2(texCR(uDisp, vUv, gridSize).x, texBS(uDisp, vUv, gridSize).y);

  float dens = clamp(d2.x * uDyeGain, 0.0, 1.0);
  vec3  dye  = uAccent * dens * dens;

  float g = smoothstep(1.15, -0.15, vUv.y + vUv.x * 0.22);
  vec3  c = uBg + uTint * g * 0.5;
  c += dye;

  // A soft shoulder rather than a clamp. Where the old expression saturated, sh was
  // exactly 1 over a whole region and the edge of that region was a level set of the
  // interpolant -- an angular, grid-following contour that read as a facet on the
  // brightest fronts, which are the ones you look at. 1 - exp(-x) is linear at the
  // onset, bounded by one, and has no contour anywhere.
  float x = max(d2.y, 0.0) * uShockGain;
  float sh = 1.0 - exp(-x);
  // The exponent, and it is a display knob only -- nothing here touches the solver.
  // A front weakens as it propagates, so its convergence falls, so a steeper
  // exponent makes an old front fade faster while a fresh one stays where it was.
  // The live page keeps the square; this page is a little steeper, because at beta =
  // 1 the fronts are longer-lived than the eye wants them to be.
  c += mix(uAccent, vec3(1.0), 0.30) * pow(sh, uShockPow) * uShockLift;

  // The field. Signed about zero so it grains rather than lights, in the chapter's
  // own accent, and scaled by its own strength so only a field that has been done
  // something to shows up. uFieldGain = 0 removes it entirely, at the cost of one
  // texture fetch.
  if (uFieldGain > 0.0) {
    vec2 L = texture(uLic, vUv).xy;
    c += uAccent * L.x * clamp(L.y, 0.0, 1.8) * uFieldGain;
  }

  float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
  c *= 1.0 / (1.0 + 0.62 * l);

  float d = length((vUv - 0.5) * vec2(uRes.x / uRes.y, 1.0));
  c *= 1.0 - uVignette * smoothstep(0.30, 1.05, d);

  c += (h1(uvec2(gl_FragCoord.xy), uint(uTime * 60.0)) - 0.5) * uGrain;
  outColor = vec4(max(c, 0.0), 1.0);
}`;

  // Diagnostics, two maxima and two sums through one reduction so the whole
  // measurement is a single readback:
  //   x = ctot, max      y = |u|^2, sum      z = rho, sum      w = |u|, max
  const F_METRICS = HEAD + EOS + `
uniform sampler2D uU, uB;
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);
  float ri = 1.0 / max(U.x, smallr);
  vec2  u  = U.yz * ri;
  float d2 = 0.5 * (dot(B.xy, B.xy) * ri + cs2);
  float cx = sqrt(d2 + sqrt(max(d2 * d2 - cs2 * B.x * B.x * ri, 0.0)));
  float cy = sqrt(d2 + sqrt(max(d2 * d2 - cs2 * B.y * B.y * ri, 0.0)));
  outColor = vec4(abs(u.x) + cx + abs(u.y) + cy, dot(u, u), U.x, length(u));
}`;

  // On demand only, so it costs nothing per frame: how well the cleaning is
  // working, and how strong the field and its current sheets are.
  //   x = |div B| dx / |B|, max     y = (div B dx)^2, sum
  //   z = |B|, sum                  w = |J| dx / |B|, max
  const F_DIVB = HEAD + `
uniform sampler2D uB;
uniform ivec2 size;
vec2 bf(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uB, p, 0).xy;
}
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec2 b0 = bf(c);
  vec2 be = bf(c + ivec2(1, 0)), bw = bf(c - ivec2(1, 0));
  vec2 bn = bf(c + ivec2(0, 1)), bs = bf(c - ivec2(0, 1));
  float div = 0.5 * ((be.x - bw.x) + (bn.y - bs.y));
  float jz  = 0.5 * ((be.y - bw.y) - (bn.x - bs.x));
  float mag = max(length(b0), 1e-6);
  outColor = vec4(abs(div) / mag, div * div, length(b0), abs(jz) / mag);
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

  // Static value noise for the LIC to average along. A texture rather than a hash
  // in the shader: eleven samples a pixel is eleven bilinear taps this way and
  // forty-four integer hashes the other.
  function makeNoise(gl, n) {
    const data = new Uint8Array(n * n);
    for (let i = 0; i < data.length; i++) data[i] = Math.random() * 256;
    const tex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, n, n, 0, gl.RED, gl.UNSIGNED_BYTE, data);
    return {
      tex,
      bind(unit) { gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, tex); return unit; }
    };
  }

  // Two attachments on one framebuffer, so the state and the field are written by
  // the same pass. The draw-buffer list is per-framebuffer state in ES 3.0, so it
  // is set once here and every single-target framebuffer keeps its default.
  function makeMRT(gl, a, b) {
    const fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, a.tex, 0);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, b.tex, 0);
    gl.drawBuffers([gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1]);
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) {
      throw new Error('incomplete MRT framebuffer');
    }
    return { fbo, w: a.w, h: a.h, u: a, f: b };
  }

  // ------------------------------------------------------------- presets

  const BASE = {
    accent: [0.62, 0.72, 0.98],
    tint: [0.020, 0.024, 0.040],
    tau: 0.10,
    mach: 0.55,
    dustGain: 1.00,
    dyeGain: 0.75,
    dyeDiss: 0.75,
    pointSize: 1.5,
    // Saturation point of the dust brightness ramp, as a multiple of the
    // population's measured mean drift. A relative number now, not the absolute
    // tau*M^2 coefficient the unmagnetised build used -- see F_PDRIFT.
    driftRef: 2.0,
    brown: 0.0,
    bhat: [1.0, 0.0],
    stream: 0.0,
    stirGain: 1.00,
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
    if (gl.getParameter(gl.MAX_DRAW_BUFFERS) < 2) throw new Error('no MRT');

    const reduced = global.matchMedia && global.matchMedia('(prefers-reduced-motion: reduce)').matches;
    const mobile = Math.min(global.innerWidth, global.innerHeight) < 720;

    const PART_SIDE = 320;
    const CBLOCK = 8;
    const SMALLR = 1e-3;
    const SMALLC = 1e-3;
    const CFL = 0.8;                 // RAMSES courant_factor
    const CS = 1.0;
    const SLOPE_TYPE = 2;            // 0 = piecewise constant, 2 = MonCen
    // Reconstruction: 'ppm' or 'plm'. Backend only, no control. PLM is the
    // reference's own configuration and is what runs while the box is being driven
    // either way -- a parabola fitted across a freshly injected discontinuity is
    // monotonised variable by variable, so a state that is monotone in every
    // variable separately can still be unphysical taken together, and with a
    // magnetic field among those variables there is more to get wrong.
    const RECON = 'ppm';
    const PLM_WINDOW = 1.6;          // wall seconds of PLM after any interaction

    // Plasma beta = P / (B^2/2). At beta = 1 with rho = cs = 1 the field is
    // sqrt(2), so vA = 1.41 and the fast speed is 1.73 across the field -- the
    // gas is stiffer than the hydrodynamic build in every direction.
    const BETA = 1.0;
    const B0 = Math.sqrt(2 * CS * CS / BETA);
    // Dedner's parabolic term, as the factor exp(-PSI_DAMP * ch * dt/dx) applied
    // once per step. Larger damps the divergence error faster and closer to where
    // it was made; too large and psi cannot carry it out of the box at all.
    const PSI_DAMP = 0.4;
    const POWELL = 1.0;              // the -(div B) B momentum source, 0 disables
    // Large-scale drag, as an inverse timescale: see the note in the solver. 1/0.625
    // sim time units, a little under a box-scale turnover at the Mach number this
    // page holds.
    const FRIC = 1.6;

    // Grain charge-to-mass ratio: charge_parameter in the reference, one value for
    // the whole population. At 100 the gyrofrequency is q|B| = 141, the gyration
    // is resolved thirty times over per step, and the gyroradius is a twelfth of a
    // cell -- the grains are stuck to field lines.
    const CHARGE = 100.0;

    // The field visualisation. FIELD_VIS = false removes the pass and the fetch;
    // setFieldVis(false) does the same at runtime, and the page binds it to a key.
    let fieldVis = true;
    const FIELD_GAIN = 0.14;     // how strongly the grain reads, at |B| = B0
    const NOISE_PX = 3.0;        // noise grain, in pixels of the LIC target
    const STEP_PX = 3.0;         // and the walk step, which has to match it
    const WAVE_LEN = 0.055;      // travelling kernel wavelength, in box heights
    const WAVE_SPEED = 0.35;     // and its rate, wavelengths per wall second
    const LIC_SCALE = 0.5;       // resolution, relative to the canvas

    // The interactions are scaled up on this page, and here is why. At beta = 1 the
    // field carries as much energy as the gas and resists being pushed around, so an
    // impulse that reads as a perturbation on the live site barely registers here --
    // measured, a scroll flick moved the rms Mach by two per cent. What it cannot do
    // at any amplitude is roll up: the mean field is vertical, the layer a scroll
    // deposits is horizontal, and field-line tension holds a layer against
    // Kelvin-Helmholtz until delta_u > 2 vA, which at beta = 1 is 2.83 sound speeds.
    // One gesture's worth here is 1.4 -- half the threshold, so the layer bows the
    // field hard and launches Alfven waves and current sheets, and the roll-up stays
    // suppressed. That is the physics of a magnetised shear layer, and BETA above is
    // the knob if you would rather have the vortices back.
    const DRIVE_GAIN = 2.0;

    // Lower than the hydrodynamic build's 0.075, because ctot is about twice as
    // large at the same Mach number: this keeps the fixed branch of
    // min(dt_default, dt_Courant) the binding one in the developed state, which is
    // what makes the pace of the picture independent of the flow's peaks.
    const DTDX_MAX = 0.060;

    // Front rendering. The 5x5 binomial spreads a one-cell spike over about two and
    // a half cells, which takes roughly a factor of two off its peak, and the soft
    // shoulder takes about another two off the mid-range where most of the visible
    // front is -- so the gain carries the first and the lift carries the second.
    // Grid-level smoothing of the convergence, in cells. The B-spline reconstruction
    // is what removes the facets; this is only here to take the front from one cell
    // wide to something a display interpolant has any information about at all, so it
    // is deliberately small -- at 1.0 the fronts went soft and stopped reading as
    // fronts.
    const SHOCK_BLUR = 0.40;
    const SHOCK_GAIN = 2.2;
    const SHOCK_LIFT = 0.36;
    const SHOCK_POW = 2.6;   // 2.0 is the live page; higher fades a propagating front faster

    // The servo's gains. The plant here is slower and less dissipative than the
    // unmagnetised one: at beta = 1 the field resists the stirring, so it takes a
    // large amplitude to get the flow moving, and then the Alfvenic cascade carries
    // that energy to the grid more slowly than a hydrodynamic one does, so the
    // amplitude that got the box going is far above the one that holds it.
    // Proportional control alone loses that race -- measured, the rms was still
    // climbing through 1.4 while the drive was being cut at its maximum rate -- so
    // there is a derivative term, which stops the servo pushing while an earlier
    // push is still arriving.
    const SERVO_D = 2.0;
    const SERVO_CLAMP = 0.35;
    const SERVO_GAIN = 0.30;
    // and a bounded output, which is all the anti-windup a multiplicative servo
    // needs: the drive is what makes the vortices and the dye plumes, so it is not
    // allowed to reach zero however far off the rms is.
    const AMP_MIN = 0.30, AMP_MAX = 60;

    // [gridH, gasStepsPerFrame, grainFrac, dprCap]. n/gridH is held constant, so
    // every rung moves the picture at the same apparent speed:
    //   rate = M * fps * n * DTDX_MAX / gridH  ~  0.058 box heights per wall
    // second at rms Mach 0.69, which is what the live site measures.
    //
    // The top rung is coarser than the live site's 176 because an MHD step costs
    // more: seven primitives through the limiter instead of four, two textures on a
    // thirteen-cell stencil instead of one, and HLLD's five waves and four square
    // roots instead of HLLC's three waves and none. At 96 cells and 2.24 steps a
    // frame the cost per frame comes out level with the live site's, which is the
    // budget this page has to fit in.
    const TIERS = [
      [ 96, 2.24, 0.90, 1.30],
      [ 80, 1.86, 0.75, 1.20],
      [ 64, 1.49, 0.55, 1.05],
      [ 56, 1.30, 0.40, 1.00]
    ];
    let tier = mobile ? 1 : 0;
    let ceiling = 0;

    const P = {
      god:    program(gl, VERT, F_GODUNOV),
      god3:   program(gl, VERT, F_GODUNOV3),
      cmax0:  program(gl, VERT, F_CMAX_STATE),
      cmaxN:  program(gl, VERT, F_CMAX_DOWN),
      stir:   program(gl, VERT, F_STIR),
      splat:  program(gl, VERT, F_SPLAT),
      shear:  program(gl, VERT, F_SHEAR),
      blast:  program(gl, VERT, F_BLAST),
      init:   program(gl, VERT, F_INIT),
      flow:   program(gl, VERT, F_FLOW),
      disp:   program(gl, VERT, F_DISP),
      part:   program(gl, VERT, F_PART),
      pdrift: program(gl, VERT, F_PDRIFT),
      pdraw:  program(gl, V_PARTICLE, F_PARTICLE),
      lic:    program(gl, VERT, F_LIC),
      dsm:    program(gl, VERT, F_DSMOOTH),
      comp:   program(gl, VERT, F_COMPOSITE),
      metric: program(gl, VERT, F_METRICS),
      divb:   program(gl, VERT, F_DIVB),
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

    let S = null, flow = null, disp = null, disp2 = null, part = null, lic = null, noise = null;
    let cmax = [], met = null, red = [], pdr = [];
    let grid = { w: 0, h: 0 }, pSide = 0, nPart = 0, nDraw = 0, dpr = 1, aspect = 1;
    let owned = [];

    const state = {
      preset: resolvePreset('hero'),
      target: resolvePreset('hero'),
      live: null,
      blend: 1,
      wall: 0, time: 0, steps: 0,
      running: true, fps: 60,
      amp: 4.0,
      ctot: 2 * CS, maxSig: 1, machRms: 0, machMax: 0, dens: 1, div: 0,
      divb: 0, divbRms: 0, bmean: B0, jmax: 0
    };
    const pending = { splats: [], shear: 0, blasts: [] };

    function activePreset() { return state.live || state.preset; }
    function machNow() { return Math.max(state.machRms, 0.5); }

    function drawQuad(f) {
      if (f) { gl.bindFramebuffer(gl.FRAMEBUFFER, f.fbo); gl.viewport(0, 0, f.w, f.h); }
      else { gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.viewport(0, 0, canvas.width, canvas.height); }
      gl.bindVertexArray(quad);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function mk(w, h, i, f, t, fil, wr, keep) {
      const o = makeFBO(gl, w, h, i, f, t, fil, wr);
      if (!keep) owned.push(o);
      return o;
    }
    // A double-buffered pair of pairs: (state, field) written together through MRT,
    // swapped together, so the two can never fall out of step.
    function mkDoubleState(w, h, keep) {
      const N = gl.NEAREST, R = gl.REPEAT;
      const a0 = mk(w, h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R, keep);
      const a1 = mk(w, h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R, keep);
      const b0 = mk(w, h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R, keep);
      const b1 = mk(w, h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R, keep);
      let A = makeMRT(gl, a0, b0), B = makeMRT(gl, a1, b1);
      return {
        w, h,
        get read() { return A; }, get write() { return B; },
        swap() { const s = A; A = B; B = s; }
      };
    }

    function eos(pr) {
      gl.uniform1f(pr.u.cs, CS);
      gl.uniform1f(pr.u.cs2, CS * CS);
      gl.uniform1f(pr.u.smallr, SMALLR);
    }
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
      grid.w = Math.max(32, Math.round(t[0] * aspect / 4) * 4);

      const R = gl.REPEAT, C = gl.CLAMP_TO_EDGE, L = gl.LINEAR, N = gl.NEAREST;
      S = mkDoubleState(grid.w, grid.h);
      // velocity and field together, filterable: the grains gather both at once
      flow = mk(grid.w, grid.h, gl.RGBA16F, gl.RGBA, gl.HALF_FLOAT, L, R);
      disp = mk(grid.w, grid.h, gl.RG16F, gl.RG, gl.HALF_FLOAT, L, R);
      disp2 = mk(grid.w, grid.h, gl.RG16F, gl.RG, gl.HALF_FLOAT, L, R);
      lic = mk(Math.max(2, Math.round(cw * LIC_SCALE)), Math.max(2, Math.round(ch * LIC_SCALE)),
               gl.RG16F, gl.RG, gl.HALF_FLOAT, L, C);
      if (!noise) noise = makeNoise(gl, 256);
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

      // Grains live for the life of the context: a tier change draws a shorter
      // prefix of the same population, which is a uniform random subsample.
      if (!part) {
        pSide = PART_SIDE;
        nPart = pSide * pSide;
        part = mkDoubleState(pSide, pSide, true);
        let pw = pSide, ph = pSide;
        do {
          pw = Math.max(1, Math.ceil(pw / CBLOCK));
          ph = Math.max(1, Math.ceil(ph / CBLOCK));
          pdr.push(mk(pw, ph, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C, true));
        } while (pw > 1 || ph > 1);
        seedParticles();
      }
      nDraw = Math.max(4096, Math.round(nPart * t[2]));

      for (let i = 0; i < prev.length; i++) { gl.deleteFramebuffer(prev[i].fbo); gl.deleteTexture(prev[i].tex); }

      reset();
      warmUp();
      writeFlow();
      driftReduce();
      paint(activePreset());
    }

    // Spin the initial condition up into developed turbulence before anything is
    // shown. The step count scales with the grid, because a box crossing costs
    // ctot/(C M) * gridH steps: 448 at gridH = 112 is the same half-crossing the
    // live site buys with 700 at 176.
    function warmUp() {
      warming = true;
      const n = Math.round(4 * grid.h);
      for (let i = 0; i < n; i++) {
        stepGas();
        if (i % 48 === 47) { measure(); servo(0.8); }
      }
      warming = false;
      cflReduce();
    }

    function reset() {
      gl.useProgram(P.init.p);
      gl.uniform1f(P.init.u.b0, B0);
      drawQuad(S.read);
      drawQuad(S.write);
      state.time = 0; state.steps = 0;
      cflReduce();
      const probe = new Float32Array(4);
      gl.bindFramebuffer(gl.FRAMEBUFFER, S.read.fbo);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, probe);
      if (!(probe[0] > 0.1)) throw new Error('state reads back as ' + probe[0]);
    }

    function seedParticles() {
      const data = new Float32Array(nPart * 4);
      for (let i = 0; i < nPart; i++) {
        data[i * 4] = Math.random();
        data[i * 4 + 1] = Math.random();
      }
      const zero = new Float32Array(nPart * 4);
      for (const f of [part.read.u, part.write.u]) {
        gl.bindTexture(gl.TEXTURE_2D, f.tex);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, pSide, pSide, 0, gl.RGBA, gl.FLOAT, data);
      }
      for (const f of [part.read.f, part.write.f]) {
        gl.bindTexture(gl.TEXTURE_2D, f.tex);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, pSide, pSide, 0, gl.RGBA, gl.FLOAT, zero);
      }
    }

    // ------------------------------------------------------------- stepping

    function cflReduce() {
      gl.useProgram(P.cmax0.p);
      eos(P.cmax0);
      gl.uniform1i(P.cmax0.u.uSrc, S.read.u.bind(0));
      gl.uniform1i(P.cmax0.u.uSrcB, S.read.f.bind(1));
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
    function dtdxNow() {
      return Math.min(dtdxCap(), Math.min(CFL / SMALLC, CFL / Math.max(state.ctot, 1e-20)));
    }
    function dtNow() { return dtdxNow() / grid.h; }

    function applyForcing(pr) {
      const dx = 1 / grid.h;
      const t = state.time;
      const c0 = [0.5 + 0.31 * Math.cos(t * 1.7 * 0.11 + 0), 0.5 + 0.27 * Math.sin(t * 1.31 * 0.11)];
      const c1 = [0.5 + 0.31 * Math.cos(t * 1.7 * 0.162 + 1), 0.5 + 0.27 * Math.sin(t * 1.31 * 0.162 + 2)];
      const a = pr.accent;
      gl.useProgram(P.stir.p);
      dtUniforms(P.stir);
      gl.uniform1i(P.stir.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.stir.u.uB, S.read.f.bind(1));
      gl.uniform1i(P.stir.u.uCmax, cmaxTex().bind(2));
      gl.uniform2f(P.stir.u.c0, c0[0], c0[1]);
      gl.uniform2f(P.stir.u.c1, c1[0], c1[1]);
      gl.uniform1f(P.stir.u.amp, state.amp);
      gl.uniform1f(P.stir.u.radius, 0.055);
      gl.uniform1f(P.stir.u.dyeRadius, 0.020);
      gl.uniform1f(P.stir.u.zeta, pr.zeta);
      gl.uniform1f(P.stir.u.aspect, aspect);
      gl.uniform1f(P.stir.u.dxCell, dx);
      gl.uniform1f(P.stir.u.smallr, SMALLR);
      gl.uniform3f(P.stir.u.dye, 2.1, 0.0, 0.0);
      drawQuad(S.write); S.swap();

      for (const s of pending.splats) {
        gl.useProgram(P.splat.p);
        gl.uniform1i(P.splat.u.uU, S.read.u.bind(0));
        gl.uniform1i(P.splat.u.uB, S.read.f.bind(1));
        gl.uniform2f(P.splat.u.point, s.x, s.y);
        gl.uniform2f(P.splat.u.delta, s.dx, s.dy);
        gl.uniform1f(P.splat.u.aspect, aspect);
        gl.uniform1f(P.splat.u.radius, s.radius);
        gl.uniform1f(P.splat.u.ink, s.ink);
        gl.uniform1f(P.splat.u.smallr, SMALLR);
        drawQuad(S.write); S.swap();
      }
      pending.splats.length = 0;

      if (Math.abs(pending.shear) > 1e-4) {
        const amp = Math.max(-6, Math.min(6, pending.shear));   // entry-point bounded; backstop
        const full = 1.5 * machNow();
        gl.useProgram(P.shear.p);
        gl.uniform1i(P.shear.u.uU, S.read.u.bind(0));
        gl.uniform1i(P.shear.u.uB, S.read.f.bind(1));
        gl.uniform1f(P.shear.u.amp, amp);
        const w = Math.max(3 / grid.h, 0.012);
        gl.uniform1f(P.shear.u.width, w);
        gl.uniform1f(P.shear.u.cycles, Math.round(0.44 * aspect / (6.2831853 * w)));
        gl.uniform1f(P.shear.u.seed, amp * 0.50);
        gl.uniform1f(P.shear.u.phase, state.time * 0.21);
        gl.uniform1f(P.shear.u.smallr, SMALLR);
        const s = Math.min(0.9, Math.abs(amp) / full) * 0.60;
        gl.uniform3f(P.shear.u.band, a[0] * s, a[1] * s, a[2] * s);
        drawQuad(S.write); S.swap();
        pending.shear *= 0.55;
        if (Math.abs(pending.shear) < 1e-3) pending.shear = 0;
      }

      for (const b of pending.blasts) {
        gl.useProgram(P.blast.p);
        gl.uniform1i(P.blast.u.uU, S.read.u.bind(0));
        gl.uniform1i(P.blast.u.uB, S.read.f.bind(1));
        gl.uniform2f(P.blast.u.point, b.x, b.y);
        gl.uniform1f(P.blast.u.aspect, aspect);
        gl.uniform1f(P.blast.u.radius, b.radius);
        gl.uniform1f(P.blast.u.dens, b.dens);
        gl.uniform1f(P.blast.u.kick, b.kick);
        gl.uniform1f(P.blast.u.smallr, SMALLR);
        drawQuad(S.write); S.swap();
      }
      pending.blasts.length = 0;
    }

    function stepGas() {
      const pr = activePreset();
      cflReduce();
      applyForcing(pr);
      // Not during the warm-up: nothing is on screen for it to resolve, and the
      // parabola nearly doubles the cost of a step, which is the whole boot stall.
      const G = (RECON === 'ppm' && agitated <= 0 && !warming) ? P.god3 : P.god;
      gl.useProgram(G.p);
      eos(G);
      dtUniforms(G);
      gl.uniform1i(G.u.uU, S.read.u.bind(0));
      gl.uniform1i(G.u.uB, S.read.f.bind(1));
      gl.uniform1i(G.u.uCmax, cmaxTex().bind(2));
      gl.uniform2i(G.u.size, grid.w, grid.h);
      gl.uniform1f(G.u.slopeType, SLOPE_TYPE);
      gl.uniform1f(G.u.dyeDiss, pr.dyeDiss);
      gl.uniform1f(G.u.dxCell, 1 / grid.h);
      gl.uniform1f(G.u.psiDamp, PSI_DAMP);
      gl.uniform1f(G.u.powell, POWELL);
      gl.uniform1f(G.u.fric, FRIC);
      drawQuad(S.write); S.swap();
      state.time += dtNow();
      state.steps++;
    }

    function writeFlow() {
      gl.useProgram(P.flow.p);
      gl.uniform1f(P.flow.u.smallr, SMALLR);
      gl.uniform1i(P.flow.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.flow.u.uB, S.read.f.bind(1));
      drawQuad(flow);

      gl.useProgram(P.disp.p);
      gl.uniform1f(P.disp.u.smallr, SMALLR);
      gl.uniform1i(P.disp.u.uU, S.read.u.bind(0));
      gl.uniform2i(P.disp.u.size, grid.w, grid.h);
      drawQuad(disp);
      smoothDisp();
    }

    function smoothDisp() {
      gl.useProgram(P.dsm.p);
      gl.uniform1i(P.dsm.u.uSrc, disp.bind(0));
      gl.uniform2i(P.dsm.u.size, grid.w, grid.h);
      gl.uniform1f(P.dsm.u.sigma, SHOCK_BLUR);
      drawQuad(disp2);
    }

    function stepDust(pr, dt) {
      const a = dt / Math.max(pr.tau, 1e-4);
      gl.useProgram(P.part.p);
      gl.uniform1i(P.part.u.uPart, part.read.u.bind(0));
      gl.uniform1i(P.part.u.uPartZ, part.read.f.bind(1));
      gl.uniform1i(P.part.u.uFlow, flow.bind(2));
      gl.uniform1f(P.part.u.dragA, a);
      gl.uniform1f(P.part.u.dt, dt);
      gl.uniform1f(P.part.u.aspect, aspect);
      gl.uniform2f(P.part.u.gridSize, grid.w, grid.h);
      gl.uniform1f(P.part.u.brown, pr.brown * (dt * 60));
      gl.uniform1f(P.part.u.charge, CHARGE);
      gl.uniform1f(P.part.u.stream, pr.stream);
      gl.uniform1ui(P.part.u.uFrame, state.steps >>> 0);
      drawQuad(part.write); part.swap();
    }

    // The population's mean drift, reduced on the GPU and never read back: the
    // vertex shader fetches the tip directly, so the lighting scale is measured
    // every time it is used and nothing has to be predicted from tau and M.
    function driftReduce() {
      let src = part.read.f, first = 1;
      for (const d of pdr) {
        gl.useProgram(P.pdrift.p);
        gl.uniform1i(P.pdrift.u.uSrc, src.bind(0));
        gl.uniform2i(P.pdrift.u.srcSize, src.w, src.h);
        gl.uniform1f(P.pdrift.u.first, first);
        drawQuad(d);
        src = d; first = 0;
      }
    }
    function driftTex() { return pdr[pdr.length - 1]; }

    // -------------------------------------------------------------- painting

    // One pass, before the composite, at half resolution.
    function drawLic() {
      gl.useProgram(P.lic.p);
      gl.uniform1i(P.lic.u.uFlow, flow.bind(0));
      gl.uniform1i(P.lic.u.uNoise, noise.bind(1));
      gl.uniform1f(P.lic.u.aspect, aspect);
      // wall time, not sim time: the packets travel at a rate you can see rather
      // than at the clock the gas happens to be running at
      gl.uniform1f(P.lic.u.time, state.wall);
      gl.uniform1f(P.lic.u.b0, B0);
      // one noise texel = NOISE_PX target pixels, and the walk steps by the same
      gl.uniform1f(P.lic.u.stepUv, STEP_PX / lic.h);
      gl.uniform1f(P.lic.u.noiseFreq, lic.h / (256.0 * NOISE_PX));
      gl.uniform1f(P.lic.u.waveLen, WAVE_LEN);
      gl.uniform1f(P.lic.u.waveSpeed, WAVE_SPEED);
      drawQuad(lic);
    }

    function composite(pr) {
      gl.useProgram(P.comp.p);
      gl.uniform1i(P.comp.u.uLic, lic.bind(3));
      gl.uniform1f(P.comp.u.uFieldGain, fieldVis ? FIELD_GAIN : 0.0);
      gl.uniform1i(P.comp.u.uDisp, disp2.bind(0));
      gl.uniform2f(P.comp.u.gridSize, grid.w, grid.h);
      gl.uniform2f(P.comp.u.uRes, canvas.width, canvas.height);
      gl.uniform3f(P.comp.u.uBg, 0.014, 0.015, 0.021);
      gl.uniform3f(P.comp.u.uTint, pr.tint[0], pr.tint[1], pr.tint[2]);
      gl.uniform3f(P.comp.u.uAccent, pr.accent[0], pr.accent[1], pr.accent[2]);
      gl.uniform1f(P.comp.u.uTime, state.wall);
      gl.uniform1f(P.comp.u.uGrain, pr.grain);
      gl.uniform1f(P.comp.u.uDyeGain, pr.dyeGain);
      gl.uniform1f(P.comp.u.uVignette, pr.vignette);
      gl.uniform1f(P.comp.u.uShockGain, SHOCK_GAIN / Math.max(0.22 * machNow(), 0.08));
      gl.uniform1f(P.comp.u.uShockLift, SHOCK_LIFT);
      gl.uniform1f(P.comp.u.uShockPow, SHOCK_POW);
      drawQuad(null);
    }

    function drawDust(pr) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, canvas.width, canvas.height);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE);
      gl.useProgram(P.pdraw.p);
      gl.uniform1i(P.pdraw.u.uPart, part.read.u.bind(0));
      gl.uniform1i(P.pdraw.u.uPartZ, part.read.f.bind(1));
      gl.uniform1i(P.pdraw.u.uDrift, driftTex().bind(2));
      gl.uniform1i(P.pdraw.u.uPW, pSide);
      gl.uniform1f(P.pdraw.u.uPointSize, pr.pointSize * dpr * 3.9);
      gl.uniform1f(P.pdraw.u.uDriftRef, pr.driftRef);
      gl.uniform1f(P.pdraw.u.uInvN, 1 / nPart);
      gl.uniform1f(P.pdraw.u.uAlpha, pr.dustGain * 0.30);
      gl.uniform1f(P.pdraw.u.uVignette, pr.vignette);
      gl.uniform1f(P.pdraw.u.uAspect, aspect);
      gl.uniform3f(P.pdraw.u.uColor, pr.accent[0], pr.accent[1], pr.accent[2]);
      gl.bindVertexArray(emptyVAO);
      gl.drawArrays(gl.POINTS, 0, nDraw);
      gl.disable(gl.BLEND);
    }

    function paint(pr) { if (fieldVis) drawLic(); composite(pr); drawDust(pr); }

    // ----------------------------------------------------------- diagnostics

    const buf = new Float32Array(4 * 64);
    let measureCountdown = 2;

    // The reduction chain and the single readback both metrics share: two maxima
    // in .x/.w and two sums in .y/.z, whatever the metric pass put there.
    function reduceAndRead() {
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
      catch (e) { return null; }
      let mx = 0, sy = 0, sz = 0, mw = 0;
      for (let i = 0; i < n; i++) {
        mx = Math.max(mx, buf[i * 4]);
        sy += buf[i * 4 + 1];
        sz += buf[i * 4 + 2];
        mw = Math.max(mw, buf[i * 4 + 3]);
      }
      return { mx, sy, sz, mw };
    }

    function measure() {
      gl.useProgram(P.metric.p);
      eos(P.metric);
      gl.uniform1i(P.metric.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.metric.u.uB, S.read.f.bind(1));
      drawQuad(met);
      const r = reduceAndRead();
      if (!r) return;
      const cells = grid.w * grid.h;
      if (isFinite(r.mw)) state.machMax = r.mw / CS;
      if (isFinite(r.sy)) state.machRms = Math.sqrt(Math.max(0, r.sy / cells)) / CS;
      if (isFinite(r.sz) && r.sz > 0) state.dens = r.sz / cells;
      if (isFinite(r.mx) && r.mx > 0) state.ctot = r.mx;
      state.maxSig = state.machMax * CS + 2 * CS;

      if (!isFinite(state.machRms) || state.machMax > 400) { state.amp = 2.0; reset(); }
    }

    // How well the cleaning is working. Off the render path: a readback drains a
    // queue that is a couple of hundred steps a second deep, so this is called by
    // hand, not per frame.
    function divbNow() {
      gl.useProgram(P.divb.p);
      gl.uniform1i(P.divb.u.uB, S.read.f.bind(0));
      gl.uniform2i(P.divb.u.size, grid.w, grid.h);
      drawQuad(met);
      const r = reduceAndRead();
      if (!r) return null;
      const cells = grid.w * grid.h;
      state.divb = r.mx;                              // max |div B| dx / |B|
      state.divbRms = Math.sqrt(Math.max(0, r.sy / cells));
      state.bmean = r.sz / cells;
      state.jmax = r.mw;
      return { maxRel: state.divb, rms: state.divbRms, bmean: state.bmean, jmax: state.jmax };
    }

    // Hold the target Mach. Proportional-derivative, and only ever called once per
    // measurement -- per step it is a runaway. The derivative is per measurement
    // interval, which is the natural unit here: it is the rate the error is
    // changing at, and subtracting it means the servo eases off while the flow is
    // still answering rather than pushing until the error changes sign.
    let servoPrev = -1;
    function servo(gain) {
      const pr = activePreset();
      const want = pr.mach * pr.stirGain;
      const s = Math.max(want, 0.4);
      const err = (want - state.machRms) / s;
      const d = servoPrev < 0 ? 0 : (state.machRms - servoPrev) / s;
      servoPrev = state.machRms;
      const u = err - SERVO_D * d;
      const g = gain * Math.max(-SERVO_CLAMP, Math.min(SERVO_CLAMP, u));
      state.amp = Math.max(AMP_MIN, Math.min(AMP_MAX, state.amp * (1 + g)));
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

    const IMP_CAP = 1.0, IMP_LEAK = 0.7;
    // Wall seconds of PLM left to run. Any interaction refreshes it.
    let agitated = 0;
    function agitate(){ agitated = PLM_WINDOW; }
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

    let credit = 0, driftCountdown = 0;
    function tick(pr) {
      const n = TIERS[tier][1];
      credit += n;
      let stepped = false;
      while (credit >= 1) { stepGas(); credit -= 1; stepped = true; }
      if (stepped) writeFlow();
      stepDust(pr, Math.max(n * dtNow(), 1e-9));
      // the lighting scale changes on the timescale the whole population changes
      // on, so it does not need refreshing every frame
      if (--driftCountdown <= 0) { driftCountdown = 8; driftReduce(); }
      paint(pr);
    }

    function bench(n) {
      const pr = activePreset();
      const px = new Float32Array(4);
      const sync = function () {
        gl.bindFramebuffer(gl.FRAMEBUFFER, S.read.fbo);
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
      impulse = Math.max(0, impulse - IMP_LEAK * machNow() * DRIVE_GAIN * dtWall);
      agitated = Math.max(0, agitated - dtWall);
      pimp = Math.max(0, pimp - SPLAT_LEAK * machNow() * dtWall);
      const pr = activePreset();

      if (!frozen) {
        tick(pr);
        if (--measureCountdown <= 0) { measureCountdown = 40; measure(); servo(SERVO_GAIN); }
      } else {
        paint(pr);
      }
      if (reduced && state.wall > 2.5) frozen = true;
    }

    canvas.addEventListener('webglcontextlost', function (e) {
      e.preventDefault(); lost = true; owned = []; part = null; S = null;
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
      // The interaction surface site.js speaks, unchanged from the live site:
      // every impulse is a multiple of the flow's own rms Mach, drawn against a
      // budget that refills at a fixed rate, so a held scroll cannot keep
      // refilling a ceiling.
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
      shear(amount) {
        const M = machNow();
        const want = Math.max(-1.5, Math.min(1.5, amount / 70)) * 0.22 * M * DRIVE_GAIN;
        const room = Math.max(0, IMP_CAP * M * DRIVE_GAIN - impulse) / 2.22;
        const s = Math.sign(want) * Math.min(Math.abs(want), room);
        if (Math.abs(s) < 1e-5) return;
        impulse += Math.abs(s) * 2.22;
        pending.shear += s;
        frozen = false;
        agitate();
      },
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
      // The field visualisation, on a switch. Off costs one texture fetch in the
      // composite and skips the LIC pass entirely.
      setFieldVis(on) { fieldVis = !!on; return fieldVis; },
      fieldVis() { return fieldVis; },
      __bench: bench,
      __divb: divbNow,
      // A column of primitives -- rho, u, v, Bx, By, psi -- for measuring a
      // profile off the clock.
      __col(ix) {
        const pu = new Float32Array(grid.h * 4), pf = new Float32Array(grid.h * 4);
        gl.bindFramebuffer(gl.FRAMEBUFFER, S.read.u.fbo);
        gl.readPixels(ix | 0, 0, 1, grid.h, gl.RGBA, gl.FLOAT, pu);
        gl.bindFramebuffer(gl.FRAMEBUFFER, S.read.f.fbo);
        gl.readPixels(ix | 0, 0, 1, grid.h, gl.RGBA, gl.FLOAT, pf);
        const out = [];
        for (let j = 0; j < grid.h; j++) {
          const r = Math.max(pu[j * 4], 1e-6);
          out.push([r, pu[j * 4 + 1] / r, pu[j * 4 + 2] / r,
            pf[j * 4], pf[j * 4 + 1], pf[j * 4 + 2]]);
        }
        return out;
      },
      // A sample of the grain population: position, velocity (three components),
      // and the drift it was lit by.
      __grains(n) {
        const side = Math.max(1, Math.min(pSide, Math.round(Math.sqrt(n || 64))));
        const a = new Float32Array(side * side * 4), b = new Float32Array(side * side * 4);
        gl.bindFramebuffer(gl.FRAMEBUFFER, part.read.u.fbo);
        gl.readPixels(0, 0, side, side, gl.RGBA, gl.FLOAT, a);
        gl.bindFramebuffer(gl.FRAMEBUFFER, part.read.f.fbo);
        gl.readPixels(0, 0, side, side, gl.RGBA, gl.FLOAT, b);
        const out = [];
        for (let i = 0; i < side * side; i++) {
          out.push([a[i * 4], a[i * 4 + 1], a[i * 4 + 2], a[i * 4 + 3], b[i * 4], b[i * 4 + 1]]);
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
          solver: 'HLLD · ' +
            (RECON === 'ppm' ? 'PPM (CW84) · MonCen slopes' : 'MonCen slopes (slope_type 2)') +
            ' · Dedner GLM · unsplit',
          mgLevels: null,
          time: state.time, sub: TIERS[tier][1], steps: state.steps, amp: state.amp,
          recon: (RECON === 'ppm' && agitated <= 0) ? 'ppm' : 'plm',
          timeScale: 1, targetFps: TARGET_FPS, cfl: CFL,
          drag: 'backward Euler after Cayley Lorentz', driveGain: DRIVE_GAIN,
          // MHD
          beta: BETA, b0: B0, charge: CHARGE, ch: state.ctot,
          psiDamp: PSI_DAMP, powell: POWELL, fric: FRIC, fieldVis: fieldVis,
          divb: state.divb, divbRms: state.divbRms, bmean: state.bmean, jmax: state.jmax,
          alfvenMach: state.machRms * CS / Math.max(state.bmean / Math.sqrt(Math.max(state.dens, 1e-6)), 1e-6)
        };
      }
    };
  }

  global.Field = {
    create(canvas) {
      try {
        // site.js holds the instance privately, and this page wants a key bound to
        // the field visualisation, so publish it. Test page; one global.
        const api = Field(canvas);
        global.__mhdField = api;
        return api;
      }
      catch (e) {
        if (global.console) console.warn('[field-mhd] ' + e.message);
        return null;
      }
    }
  };
})(window);

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
  ivec2 o = ivec2(gl_FragCoord.xy) * 16;
  float m = 0.0;
  for (int j = 0; j < 16; j++) {
    for (int i = 0; i < 16; i++) {
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
  ivec2 o = ivec2(gl_FragCoord.xy) * 16;
  float m = 0.0;
  for (int j = 0; j < 16; j++) {
    for (int i = 0; i < 16; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      m = max(m, texelFetch(uSrc, p, 0).x);
    }
  }
  outColor = vec4(m, 0.0, 0.0, 1.0);
}`;

  const DTDX = `
uniform sampler2D uCmax;
uniform float cfl, smallc;
float stepDtDx(){
  float ctot = max(texelFetch(uCmax, ivec2(0, 0), 0).x, 1e-20);
  return min(cfl / smallc, cfl / ctot);
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
uniform float amp, radius, zeta, aspect, dxCell, smallr;

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
  // the stirrers leave a little dye behind them, as they did before
  float w = exp(-dot((vUv - c0) * vec2(aspect, 1.0), (vUv - c0) * vec2(aspect, 1.0)) / radius)
          + exp(-dot((vUv - c1) * vec2(aspect, 1.0), (vUv - c1) * vec2(aspect, 1.0)) / radius);
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
  // transverse seed: the canonical Kelvin-Helmholtz setup, so the roll-up is
  // real. Applied to momentum, because momentum is what is conserved here.
  const F_SHEAR = HEAD + `
uniform sampler2D uU;
uniform vec3  band;
uniform float amp, width, seed, phase, smallr;
void main(){
  vec4  U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  float r = max(U.x, smallr);
  float y = vUv.y - 0.5;
  float env = exp(-pow(y / (2.6 * width), 2.0));
  vec2  du = vec2(amp * tanh(y / width) * env,
                  seed * sin(6.2831853 * (vUv.x * 5.0 + phase)) * env);
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
  float Y = clamp(0.45 + 0.4 * sin(TAU * (uv.x + 0.3)) * cos(TAU * uv.y), 0.0, 1.0);
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

  const F_PART = HEAD + SPECTRUM + `
uniform sampler2D uPart, uVel;
uniform vec2  stream;
uniform float dragA, dt, aspect, brown;
uniform uint  uFrame;

void main(){
  vec4  P   = texture(uPart, vUv);
  vec2  pos = P.xy;
  vec2  vel = P.zw;
  uvec2 id  = uvec2(gl_FragCoord.xy);

  vec2 ug = texture(uVel, pos).xy;
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
uniform float uPointSize, uDriftNorm, uAlpha, uVignette, uAspect;
out float vBright;
out float vDrift;
` + SPECTRUM + `
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
  vec2 ug = texture(uVel, pos).xy;
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

  const F_COMPOSITE = HEAD + EOS + `
uniform sampler2D uU;
uniform ivec2 size;
uniform vec2  uRes;
uniform vec3  uBg, uTint, uAccent;
uniform float uTime, uGrain, uDyeGain, uVignette, uShockGain;

uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float h1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}

vec4 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0);
}
// Bilinear by hand, for the same reason F_VEL exists: RGBA32F does not filter.
vec4 at(vec2 uv){
  vec2  g  = uv * vec2(size) - 0.5;
  ivec2 i0 = ivec2(floor(g));
  vec2  t  = g - vec2(i0);
  return mix(mix(cell(i0),               cell(i0 + ivec2(1, 0)), t.x),
             mix(cell(i0 + ivec2(0, 1)), cell(i0 + ivec2(1, 1)), t.x), t.y);
}

void main(){
  vec4  q  = toPrim(at(vUv));
  vec2  tx = 1.0 / vec2(size);
  vec4  qe = toPrim(at(vUv + vec2(tx.x, 0.0))), qw = toPrim(at(vUv - vec2(tx.x, 0.0)));
  vec4  qn = toPrim(at(vUv + vec2(0.0, tx.y))), qs = toPrim(at(vUv - vec2(0.0, tx.y)));

  // Dye, through a squared density response: faint gas stays dark and only
  // genuinely dense wisps register, so the gas sits behind the dust.
  float dens = clamp(q.w * uDyeGain, 0.0, 1.0);
  vec3  dye  = uAccent * dens * dens;

  // faint cold gradient so the page never reads as flat black
  float g = smoothstep(1.15, -0.15, vUv.y + vUv.x * 0.22);
  vec3  c = uBg + uTint * g * 0.5;
  c += dye;

  // Convergence, which is the thing an incompressible background cannot show at
  // all: div u < 0 is a shock, and here it is a real discontinuity the Riemann
  // solver captured rather than a gradient in a projected field.
  float div = 0.5 * ((qe.y - qw.y) + (qn.z - qs.z));
  float sh  = clamp(-div * uShockGain, 0.0, 1.0);
  // A highlight on the fronts, not a second background. Cubed, so only genuinely
  // strong convergence registers -- squared covered half the frame.
  c += vec3(1.0, 0.74, 0.46) * sh * sh * sh * 0.42;

  float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
  c *= 1.0 / (1.0 + 0.62 * l);

  float d = length((vUv - 0.5) * vec2(uRes.x / uRes.y, 1.0));
  c *= 1.0 - uVignette * smoothstep(0.30, 1.05, d);

  c += (h1(uvec2(gl_FragCoord.xy), uint(uTime * 60.0)) - 0.5) * uGrain;
  outColor = vec4(max(c, 0.0), 1.0);
}`;

  // Diagnostics: max ctot for the readout, plus sums for rms Mach and mean
  // density. Separate from the CFL chain, which runs every step and is max-only.
  const F_METRICS = HEAD + EOS + `
uniform sampler2D uU;
void main(){
  vec4  q = toPrim(texelFetch(uU, ivec2(gl_FragCoord.xy), 0));
  outColor = vec4(length(q.yz), dot(q.yz, q.yz), q.x, 1.0);
}`;

  const F_REDUCE = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 4;
  float mx = 0.0;
  vec3  s  = vec3(0.0);
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      vec4 v = texelFetch(uSrc, p, 0);
      mx = max(mx, v.x);
      s += v.yzw;
    }
  }
  outColor = vec4(mx, s);
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
    mach: 2.2,          // target rms Mach, scaled by stirGain
    dustGain: 1.00,
    dyeGain: 0.75,
    dyeDiss: 0.55,
    pointSize: 1.5,
    driftRef: 0,        // 0 = derive from tau and the measured Mach
    brown: 0.0,
    bhat: [1.0, 0.0],
    stream: 0.0,        // grain bulk streaming, in units of cs
    stirGain: 1.00,
    zeta: 0.30,         // compressive fraction of the stirring
    grain: 0.030,
    vignette: 0.62
  };

  const PRESETS = {
    hero:      { accent: [0.50, 0.64, 1.00], tau: 0.14, dustGain: 1.48, stirGain: 1.25 },
    bio:       { accent: [0.96, 0.80, 0.50], tau: 0.10, dustGain: 1.05 },
    picdust:   { accent: [1.00, 0.58, 0.18], tau: 0.030, dustGain: 1.45, pointSize: 1.2, stirGain: 1.15 },
    dfmm:      { accent: [0.26, 0.95, 0.72], tau: 0.85,  dustGain: 2.24, pointSize: 1.75, brown: 0.02 },
    mhd:       { accent: [0.60, 0.46, 1.00], tau: 0.12,  dustGain: 1.81, bhat: [0.94, 0.34], zeta: 0.18 },
    cosmicray: { accent: [0.22, 0.86, 1.00], tau: 0.020, dustGain: 2.40, stream: 4.8, bhat: [0.94, 0.34], pointSize: 1.15, driftRef: 6.5 },
    phrike:    { accent: [1.00, 0.24, 0.48], tau: 0.06,  dustGain: 1.38, stirGain: 1.35, zeta: 0.45 },
    papers:    { accent: [0.94, 0.82, 0.52], tau: 0.09,  dustGain: 0.95, stirGain: 0.85 },
    blog:      { accent: [0.56, 0.72, 1.00], tau: 0.11,  dustGain: 0.95, stirGain: 0.8 },
    iron:      { accent: [1.00, 0.13, 0.10], tau: 0.55,  dustGain: 1.00, pointSize: 1.6, stirGain: 0.65, vignette: 0.74, zeta: 0.5 },
    cactus:    { accent: [0.54, 0.84, 0.32], tau: 0.95,  dustGain: 2.06, pointSize: 1.7, brown: 0.04, stirGain: 0.42, grain: 0.040 },
    skate:     { accent: [1.00, 0.48, 0.06], tau: 0.20,  dustGain: 1.19, pointSize: 1.3, stirGain: 1.1 },
    contact:   { accent: [0.92, 0.86, 0.58], tau: 0.10,  dustGain: 1.05, stirGain: 0.9 }
  };

  function resolvePreset(key) {
    const p = Object.assign({}, BASE, PRESETS[key] || PRESETS.hero);
    if (!p.driftRef) p.driftRef = 1.15 + 2.4 * p.tau;
    return p;
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
    const CBLOCK = 16;
    // 1e-4 was reachable: at Mach 5 the rarefactions span five decades and pin
    // the floor. Three decades below the mean is still far outside the physics.
    const SMALLR = 1e-3;
    const SMALLC = 1e-3;
    const CFL = 0.8;                 // RAMSES courant_factor
    const CS = 1.0;

    // [gridH, substeps, grainFrac, dprCap].
    //
    // These grids are much coarser than the incompressible build's 192, and that
    // is the central finding of this experiment rather than a shortcut. An
    // explicit compressible scheme is CFL-limited by the *acoustic* speed, so it
    // needs about ctot/(C cs) * gridH ~ 20 * gridH steps per box crossing time. A
    // projection method needs none of them: its pressure is elliptic, so it can
    // step on the advective timescale alone. At gridH = 176 this background
    // managed t = 0.33 crossing times in twenty-five seconds of wall clock, which
    // is not a background, it is a screensaver of a still image. Coarsening the
    // gas buys dt linearly and pass cost quadratically, and the dust is drawn at
    // native resolution regardless, so the visible loss is fatter shocks.
    const TIERS = [
      [112, 5, 1.00, 1.35],
      [ 96, 4, 0.85, 1.20],
      [ 80, 3, 0.65, 1.05],
      [ 64, 2, 0.45, 1.00]
    ];
    let tier = mobile ? 3 : 1;
    let ceiling = 0;

    const P = {
      god:    program(gl, VERT, F_GODUNOV),
      cmax0:  program(gl, VERT, F_CMAX_STATE),
      cmaxN:  program(gl, VERT, F_CMAX_DOWN),
      stir:   program(gl, VERT, F_STIR),
      splat:  program(gl, VERT, F_SPLAT),
      shear:  program(gl, VERT, F_SHEAR),
      blast:  program(gl, VERT, F_BLAST),
      init:   program(gl, VERT, F_INIT),
      vel:    program(gl, VERT, F_VEL),
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

    let U = null, vel = null, part = null, cmax = [], met = null, red = [];
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
      maxSig: 1, machRms: 0, machMax: 0, dens: 1, div: 0
    };
    const pending = { splats: [], shear: 0, blasts: [] };

    function activePreset() { return state.live || state.preset; }

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
    function dtUniforms(pr) {
      gl.uniform1f(pr.u.cfl, CFL);
      gl.uniform1f(pr.u.smallc, SMALLC);
    }

    let allocW = -1, allocH = -1, allocTier = -1, warmed = false;

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
      if (!warmed) { warmed = true; for (let i = 0; i < 900; i++) stepGas(); }
      writeVel();
      paint(activePreset());
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
    function dtdxNow() { return Math.min(CFL / SMALLC, CFL / Math.max(state.maxSig, 1e-20)); }
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
      gl.uniform1f(P.stir.u.zeta, pr.zeta);
      gl.uniform1f(P.stir.u.aspect, aspect);
      gl.uniform1f(P.stir.u.dxCell, dx);
      gl.uniform1f(P.stir.u.smallr, SMALLR);
      // enough ink that the concentration reaches order unity in the plumes
      gl.uniform3f(P.stir.u.dye, 2.6, 0.0, 0.0);
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
        const amp = Math.max(-6, Math.min(6, pending.shear));
        gl.useProgram(P.shear.p);
        gl.uniform1i(P.shear.u.uU, U.read.bind(0));
        gl.uniform1f(P.shear.u.amp, amp);
        gl.uniform1f(P.shear.u.width, 0.055);
        gl.uniform1f(P.shear.u.seed, amp * 0.30);
        gl.uniform1f(P.shear.u.phase, state.time * 0.21);
        gl.uniform1f(P.shear.u.smallr, SMALLR);
        const s = Math.min(0.9, Math.abs(amp) / 6) * 0.25;
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
      gl.useProgram(P.god.p);
      eos(P.god);
      dtUniforms(P.god);
      gl.uniform1i(P.god.u.uU, U.read.bind(0));
      gl.uniform1i(P.god.u.uCmax, cmaxTex().bind(1));
      gl.uniform2i(P.god.u.size, grid.w, grid.h);
      gl.uniform1f(P.god.u.llf, 0);
      gl.uniform1f(P.god.u.dyeDiss, pr.dyeDiss);
      gl.uniform1f(P.god.u.dxCell, 1 / grid.h);
      drawQuad(U.write); U.swap();
      state.time += dtNow();
      state.steps++;
    }

    function writeVel() {
      gl.useProgram(P.vel.p);
      eos(P.vel);
      gl.uniform1i(P.vel.u.uU, U.read.bind(0));
      drawQuad(vel);
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
      gl.uniform1f(P.part.u.brown, pr.brown * (dt * 60));
      gl.uniform1ui(P.part.u.uFrame, state.steps >>> 0);
      gl.uniform2f(P.part.u.stream, pr.stream * b[0] / bn, pr.stream * b[1] / bn);
      drawQuad(part.write); part.swap();
    }

    // -------------------------------------------------------------- painting

    function composite(pr) {
      gl.useProgram(P.comp.p);
      eos(P.comp);
      gl.uniform1i(P.comp.u.uU, U.read.bind(0));
      gl.uniform2i(P.comp.u.size, grid.w, grid.h);
      gl.uniform2f(P.comp.u.uRes, canvas.width, canvas.height);
      gl.uniform3f(P.comp.u.uBg, 0.014, 0.015, 0.021);
      gl.uniform3f(P.comp.u.uTint, pr.tint[0], pr.tint[1], pr.tint[2]);
      gl.uniform3f(P.comp.u.uAccent, pr.accent[0], pr.accent[1], pr.accent[2]);
      gl.uniform1f(P.comp.u.uTime, state.wall);
      gl.uniform1f(P.comp.u.uGrain, pr.grain);
      gl.uniform1f(P.comp.u.uDyeGain, pr.dyeGain);
      gl.uniform1f(P.comp.u.uVignette, pr.vignette);
      // div is du per cell; a shock is a jump of order cs across a cell
      gl.uniform1f(P.comp.u.uShockGain, 1.0 / (0.9 * CS));
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
      gl.uniform1f(P.pdraw.u.uPointSize, pr.pointSize * dpr * 3.9);
      // the drift scale rides on the measured Mach number, so the lighting keeps
      // its meaning as the driving changes
      const ref = pr.driftRef * Math.max(state.machRms, 0.6);
      gl.uniform1f(P.pdraw.u.uDriftNorm, 1 / Math.max(ref, 0.2));
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
      let mx = 0, sq = 0, sm = 0;
      for (let i = 0; i < n; i++) {
        mx = Math.max(mx, buf[i * 4]);
        sq += buf[i * 4 + 1];
        sm += buf[i * 4 + 2];
      }
      const cells = grid.w * grid.h;
      if (isFinite(mx)) state.machMax = mx / CS;
      if (isFinite(sq)) state.machRms = Math.sqrt(Math.max(0, sq / cells)) / CS;
      if (isFinite(sm) && sm > 0) state.dens = sm / cells;
      state.maxSig = state.machMax * CS + 2 * CS;

      // hold the target Mach: once per measurement, rate-limited, exactly as the
      // compressible sub-page does. Running it per step is a runaway.
      const pr = activePreset();
      const want = pr.mach * pr.stirGain;
      const rel = (want - state.machRms) / Math.max(want, 0.5);
      state.amp = Math.max(0.05, Math.min(400, state.amp * (1 + 0.45 * Math.max(-0.35, Math.min(0.35, rel)))));
      if (!isFinite(state.machRms) || state.machMax > 400) { state.amp = 2.0; reset(); }
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
      try { allocate(true); } catch (e) { frozen = true; dead = true; }
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
      const pr = activePreset();

      if (!frozen) {
        const sub = TIERS[tier][1];
        const t0 = state.time;
        for (let i = 0; i < sub; i++) stepGas();
        writeVel();
        // one dust update per frame with the gas time actually elapsed, which is
        // safe because backward Euler is L-stable for any step
        stepDust(pr, Math.max(state.time - t0, 1e-8));
        if (--measureCountdown <= 0) { measureCountdown = 20; measure(); }
      }
      paint(pr);
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
      // Momentum in units of cs, ink in dye units. The incompressible build took
      // hundreds of cells per second here; this one takes Mach numbers.
      splat(x, y, dx, dy, radius) {
        const sp = Math.hypot(dx, dy);
        const cap = 4.0;
        const k = sp > cap ? cap / sp : 1;
        const m = Math.min(1, sp / cap);
        pending.splats.push({
          x, y, dx: dx * k, dy: dy * k,
          ink: 0.06 + m * 0.30,
          radius: radius || 0.0030
        });
        if (pending.splats.length > 6) pending.splats.splice(0, pending.splats.length - 6);
        frozen = false;
      },
      shear(amount) { pending.shear += amount; frozen = false; },
      blast(x, y, amp, radius) {
        pending.blasts.push({ x, y, dens: 2.2, kick: Math.min(6, (amp || 300) / 60), radius: radius || 0.010 });
        frozen = false;
      },
      setPreset(key) {
        const next = resolvePreset(key);
        state.preset = Object.assign({}, activePreset());
        state.target = next;
        state.blend = 0;
        frozen = false;
      },
      accentOf(key) { return resolvePreset(key).accent; },
      stats() {
        const pr = activePreset();
        return {
          gridW: grid.w, gridH: grid.h, nPart, nDraw, fps: state.fps, tier,
          tau: pr.tau, dpr,
          rms: state.machRms, max: state.machMax, div: state.div,
          mach: state.machRms, dens: state.dens,
          solver: 'HLL · piecewise constant · unsplit',
          mgLevels: null,
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

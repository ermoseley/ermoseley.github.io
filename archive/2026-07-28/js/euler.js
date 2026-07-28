/* ============================================================================
   euler.js — compressible, for once.

   The background on the main page is incompressible: a projection method, an
   elliptic pressure, no shocks anywhere by construction. This is the other half
   of the subject. It solves the isothermal Euler equations

       d(rho)/dt   + div(rho u)                 = 0
       d(rho u)/dt + div(rho u u + rho cs^2 I)  = rho a

   conservatively, with a Godunov scheme, on a periodic square grid:

     reconstruct : piecewise linear, minmod-limited, on the conserved variables
     predict     : MUSCL-Hancock half step, F(U^L) - F(U^R), so the scheme is
                   second order in time without storing a second stage
     flux        : HLL, Davis wave-speed estimate SL = min(uL,uR) - cs,
                   SR = max(uL,uR) + cs -- positivity-preserving for isothermal
     update      : conservative differencing of face fluxes
     directions  : dimensionally split, Strang-alternated each step, so the
                   splitting error stays second order

   Isothermal on purpose. It is the standard closure for the cold neutral and
   molecular ISM, where line cooling pins the temperature over the dynamical
   time; it also makes the sound speed a constant, so the CFL condition is
   honest and cheap, and the state fits in three channels instead of four.

   The timestep is set by the measured maximum signal speed, not guessed:
   dt = CFL * dx / max(|u| + cs), remeasured on the GPU and reduced to a handful
   of texels before it comes back to the CPU.

   Everything you can see is a consequence of that: shocks are real
   discontinuities that the Riemann solver captures, the density PDF is
   log-normal because isothermal supersonic turbulence makes it so, and the
   width of that log-normal follows sigma^2 = ln(1 + b^2 M^2) without anyone
   putting it there.

   Written from scratch. No libraries.
   ========================================================================= */

(function (global) {
  'use strict';

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

  // --------------------------------------------------------------- the solver

  // One dimensionally-split MUSCL-Hancock sweep with an HLL flux. `dir` picks
  // the sweep axis, so the same shader serves x and y. The stencil is i-2..i+2:
  // updating cell i needs the half-step boundary states of i-1 and i+1, which
  // need their slopes, which need their neighbours.
  const F_SWEEP = HEAD + `
uniform sampler2D uU;
uniform ivec2 size;
uniform ivec2 dir;          // (1,0) for an x sweep, (0,1) for a y sweep
uniform float dtdx, cs2, rhoMin;

vec3 fetch(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0).xyz;
}

// Conserved -> flux along the sweep axis. The state is carried as
// (rho, m_par, m_perp) so one expression covers both directions.
vec3 flux(vec3 U){
  float rho = max(U.x, rhoMin);
  float up  = U.y / rho;
  return vec3(U.y, U.y * up + rho * cs2, U.z * up);
}

float minmod(float a, float b){
  return (a * b <= 0.0) ? 0.0 : (abs(a) < abs(b) ? a : b);
}
vec3 slope(vec3 l, vec3 c, vec3 r){
  vec3 a = c - l, b = r - c;
  return vec3(minmod(a.x, b.x), minmod(a.y, b.y), minmod(a.z, b.z));
}

// Rotate so component .y is always along the sweep axis.
vec3 toLocal(vec3 U){ return dir.x == 1 ? U : vec3(U.x, U.z, U.y); }
vec3 toWorld(vec3 U){ return dir.x == 1 ? U : vec3(U.x, U.z, U.y); }

// HLL. With a constant sound speed the Davis estimate is exact enough and
// keeps rho positive for any admissible pair of states.
vec3 hll(vec3 UL, vec3 UR){
  float rL = max(UL.x, rhoMin), rR = max(UR.x, rhoMin);
  float uL = UL.y / rL,         uR = UR.y / rR;
  float cs = sqrt(cs2);
  float SL = min(uL, uR) - cs;
  float SR = max(uL, uR) + cs;
  if (SL >= 0.0) return flux(UL);
  if (SR <= 0.0) return flux(UR);
  vec3 FL = flux(UL), FR = flux(UR);
  return (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL);
}

// Boundary-extrapolated states of one cell, already advanced half a step.
void hancock(ivec2 c, out vec3 sL, out vec3 sR){
  vec3 Um = toLocal(fetch(c - dir));
  vec3 U0 = toLocal(fetch(c));
  vec3 Up = toLocal(fetch(c + dir));
  vec3 d  = slope(Um, U0, Up);
  vec3 L  = U0 - 0.5 * d;
  vec3 R  = U0 + 0.5 * d;
  // MUSCL-Hancock: evolve both interface states by dt/2 using the *same*
  // flux difference, which is what makes this second order in time with no
  // intermediate storage.
  vec3 h  = 0.5 * dtdx * (flux(L) - flux(R));
  sL = L + h;
  sR = R + h;
}

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);

  vec3 aL, aR, bL, bR, cL, cR;
  hancock(c - dir, aL, aR);          // left neighbour
  hancock(c,       bL, bR);          // this cell
  hancock(c + dir, cL, cR);          // right neighbour

  vec3 Fm = hll(aR, bL);             // face at i-1/2
  vec3 Fp = hll(bR, cL);             // face at i+1/2

  vec3 U = toLocal(fetch(c)) - dtdx * (Fp - Fm);
  U.x = max(U.x, rhoMin);
  outColor = vec4(toWorld(U), 1.0);
}`;

  // Large-scale driving, as an Ornstein-Uhlenbeck process.
  //
  // The first version of this used a handful of modes with steadily rotating
  // phases, which is coherent forcing: it compresses the same places over and
  // over, and the density contrast it produced was far larger than the physics
  // warrants. Measured b came out near 4 in the purely compressive limit against
  // an expected 1. Real drivers decorrelate, so this one does too.
  //
  // Each mode carries a complex amplitude driven by
  //     da = -a dt/tc + sqrt(2 dt/tc) dW
  // with tc the box eddy-turnover time, which gives unit variance in the steady
  // state and a finite correlation time. The amplitudes live on the CPU (eight
  // modes is sixteen numbers) and arrive as a uniform array.
  //
  // Each mode is split into a solenoidal part -- the curl of a stream function,
  // divergence-free -- and a compressive part -- the gradient of a potential,
  // curl-free. zeta mixes them: 0 purely solenoidal, 1 purely compressive. That
  // is the knob that moves b in sigma^2 = ln(1 + b^2 M^2).
  const NMODE = 8;
  const F_DRIVE = HEAD + `
uniform sampler2D uU;
uniform float amp, zeta, dt;
uniform vec2  kvec[8];   // integer wavenumbers, so exactly periodic
uniform vec2  amps[8];   // complex OU amplitude per mode

const float TAU = 6.28318530718;

void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 sol = vec2(0.0), comp = vec2(0.0);
  for (int i = 0; i < 8; i++) {
    vec2  k = kvec[i];
    float kk = length(k);
    if (kk < 1e-6) continue;
    float th = TAU * dot(k, vUv);
    // S is d/d(theta) of Re(a e^{i theta}), so grad of that field is -TAU k S
    float S  = amps[i].x * sin(th) + amps[i].y * cos(th);
    vec2  kh = k / kk;
    sol  += S * vec2(-kh.y, kh.x);
    comp += S * -kh;
  }
  float norm = 1.0 / sqrt(float(8));
  // acceleration, so the force on a cell is rho * a: heavy gas is harder to
  // push, which is part of why compressive driving widens the density PDF
  vec2 a = mix(sol, comp, zeta) * amp * norm;
  outColor = vec4(U.x, U.yz + U.x * a * dt, U.w);
}`;

  // A blast: a local overdensity plus the outward momentum to match. In an
  // isothermal gas pressure is rho cs^2, so piling up density *is* piling up
  // pressure, and this relaxes into a real outward-running shock rather than a
  // painted-on ring.
  const F_BLAST = HEAD + `
uniform sampler2D uU;
uniform vec2  point;
uniform float aspect, radius, dens, kick;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float r = length(p) + 1e-5;
  float w = exp(-dot(p, p) / radius);
  float rho = U.x + dens * w;
  vec2  m   = U.yz + kick * w * (p / r) * U.x;
  outColor = vec4(rho, m, U.w);
}`;

  // A stir: momentum along a stroke, no density change.
  const F_PUSH = HEAD + `
uniform sampler2D uU;
uniform vec2  point, delta;
uniform float aspect, radius;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float w = exp(-dot(p, p) / radius);
  outColor = vec4(U.x, U.yz + delta * w * U.x, U.w);
}`;

  const F_INIT = HEAD + `
uniform float seed;
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
  uvec2 id = uvec2(gl_FragCoord.xy);
  // uniform gas at rest, with a whisper of noise so the driving has something
  // to work on rather than starting from an exactly symmetric state
  float n = h1(id, uint(seed)) - 0.5;
  outColor = vec4(1.0 + 0.02 * n, 0.0, 0.0, 1.0);
}`;

  // Per-cell diagnostics, gathered so one readback answers everything.
  //   R: |u| + cs   reduced by MAX -> sets the CFL-limited timestep
  //   G: |u|^2      reduced by sum -> rms Mach number
  //   B: rho        reduced by sum -> TOTAL MASS, which a conservative scheme
  //                 must hold to round-off; this is the check that the flux
  //                 differencing is doing what it claims
  //   A: ln rho     reduced by sum -> mean of ln rho
  // sigma is estimated from the decimated sample instead, which has 9216 points
  // and is plenty for a second moment; the divergence comes from there too.
  const F_METRICS = HEAD + `
uniform sampler2D uU;
uniform float cs2, rhoMin;
void main(){
  vec3  U   = texelFetch(uU, ivec2(gl_FragCoord.xy), 0).xyz;
  float rho = max(U.x, rhoMin);
  vec2  u   = U.yz / rho;
  float q   = dot(u, u);
  float l   = log(rho);
  outColor = vec4(sqrt(q) + sqrt(cs2), q, rho, l);
}`;

  // 4x4 gather. Max on R, sum on the rest.
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

  // Decimation, not averaging: a box mean would change the very distribution
  // this exists to measure. Point samples preserve the marginal PDF of ln rho.
  // Carries the divergence too, so one readback yields both the PDF and the
  // fraction of the volume that is actively shocking.
  const F_DECIMATE = HEAD + `
uniform sampler2D uU;
uniform ivec2 stride, size;
uniform float rhoMin, cs2;
vec3 at(ivec2 c){
  ivec2 q = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, q, 0).xyz;
}
void main(){
  ivec2 p = ivec2(gl_FragCoord.xy) * stride;
  vec3  U = at(p);
  float rho = max(U.x, rhoMin);
  vec3  R = at(p + ivec2(1, 0)), L = at(p - ivec2(1, 0));
  vec3  T = at(p + ivec2(0, 1)), B = at(p - ivec2(0, 1));
  // du per cell. A shock is a jump of order cs across about one cell, so this
  // is already in the natural units for deciding what counts as one.
  float div = 0.5 * ((R.y / max(R.x, rhoMin) - L.y / max(L.x, rhoMin))
                   + (T.z / max(T.x, rhoMin) - B.z / max(B.x, rhoMin)));
  outColor = vec4(log(rho), rho, div, length(U.yz / rho) / sqrt(cs2));
}`;

  // Display. Density on a log stretch is the base image, because that is what
  // the physics is log-normal in; convergence (-div u) is overlaid, because
  // that is where the shocks are.
  const F_SHOW = HEAD + `
uniform sampler2D uU;
uniform vec2  uRes;
uniform vec3  uC0, uC1, uC2, uC3, uShock;
uniform float cs2, rhoMin, logLo, logHi, shockGain, shockMin, mode, vignette;

// Bilinear by hand. The grid is coarse next to the canvas and point-sampling it
// turns every shock into a staircase -- but RGBA32F is not a filterable format
// unless OES_texture_float_linear is present, and setting LINEAR on a texture
// whose format cannot filter makes it INCOMPLETE, after which every read of it
// returns zero, texelFetch included. So the interpolation is written out.
uniform ivec2 size;
vec3 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0).xyz;
}
vec3 at(vec2 uv){
  vec2  g  = uv * vec2(size) - 0.5;
  ivec2 i0 = ivec2(floor(g));
  vec2  t  = g - vec2(i0);
  return mix(mix(cell(i0),              cell(i0 + ivec2(1, 0)), t.x),
             mix(cell(i0 + ivec2(0, 1)), cell(i0 + ivec2(1, 1)), t.x), t.y);
}
vec3 ramp(float t){
  t = clamp(t, 0.0, 1.0);
  if (t < 0.3333) return mix(uC0, uC1, t / 0.3333);
  if (t < 0.6666) return mix(uC1, uC2, (t - 0.3333) / 0.3333);
  return mix(uC2, uC3, (t - 0.6666) / 0.3334);
}
void main(){
  vec3  U   = at(vUv);
  float rho = max(U.x, rhoMin);
  vec2  u   = U.yz / rho;

  vec2  tx = 1.0 / vec2(size);
  vec3  R = at(vUv + vec2(tx.x, 0.0)), L = at(vUv - vec2(tx.x, 0.0));
  vec3  T = at(vUv + vec2(0.0, tx.y)), B = at(vUv - vec2(0.0, tx.y));
  float div = 0.5 * ((R.y / max(R.x, rhoMin) - L.y / max(L.x, rhoMin))
                   + (T.z / max(T.x, rhoMin) - B.z / max(B.x, rhoMin)));

  float ln = log(rho);
  float dn = (ln - logLo) / max(logHi - logLo, 1e-4);
  float mach = length(u) / sqrt(cs2);

  vec3 col;
  if (mode < 0.5) {
    col = ramp(dn);
  } else if (mode < 1.5) {
    // convergence only: the shock network on its own
    col = uShock * clamp((-div - shockMin) * shockGain, 0.0, 1.4);
  } else if (mode < 2.5) {
    col = ramp(clamp(mach / 6.0, 0.0, 1.0));
  } else {
    // divergence, signed: compression warm, expansion cool
    float s = clamp(div * shockGain * 0.5, -1.0, 1.0);
    col = mix(uC1 * 0.6, uShock, 0.5 + 0.5 * (-s));
    col *= abs(s) * 0.9 + 0.05;
  }

  if (mode < 0.5) {
    // shocks laid over the density field, keyed on convergence
    float sh = clamp((-div - shockMin) * shockGain, 0.0, 1.0);
    col += uShock * sh * sh * 0.85;
  }

  float d = length((vUv - 0.5) * vec2(uRes.x / uRes.y, 1.0));
  col *= 1.0 - vignette * smoothstep(0.32, 1.05, d);

  float l = dot(col, vec3(0.2126, 0.7152, 0.0722));
  col *= 1.0 / (1.0 + 0.35 * l);
  outColor = vec4(max(col, 0.0), 1.0);
}`;

  // ------------------------------------------------------------- gl plumbing

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
      const info = gl.getActiveUniform(p, i);
      const nm = info.name.replace(/\[0\]$/, '');
      u[nm] = gl.getUniformLocation(p, nm);
    }
    return { p, u };
  }

  function makeFBO(gl, w, h, filter) {
    const tex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, w, h, 0, gl.RGBA, gl.FLOAT, null);
    const fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) {
      throw new Error('incomplete framebuffer');
    }
    return {
      tex, fbo, w, h,
      bind(unit) { gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, tex); return unit; }
    };
  }

  // ------------------------------------------------------------------ solver

  function Euler(canvas, opts) {
    opts = opts || {};
    const gl = canvas.getContext('webgl2', {
      alpha: false, depth: false, stencil: false, antialias: false,
      premultipliedAlpha: false, preserveDrawingBuffer: false
    });
    if (!gl) throw new Error('no webgl2');
    if (!gl.getExtension('EXT_color_buffer_float')) throw new Error('no float render targets');

    const P = {
      sweep:  program(gl, VERT, F_SWEEP),
      drive:  program(gl, VERT, F_DRIVE),
      blast:  program(gl, VERT, F_BLAST),
      push:   program(gl, VERT, F_PUSH),
      init:   program(gl, VERT, F_INIT),
      metric: program(gl, VERT, F_METRICS),
      reduce: program(gl, VERT, F_REDUCE),
      decim:  program(gl, VERT, F_DECIMATE),
      show:   program(gl, VERT, F_SHOW)
    };

    const vao = gl.createVertexArray();
    gl.bindVertexArray(vao);
    const vb = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vb);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);

    const RHO_MIN = 1e-4;
    const PDF_W = 128, PDF_H = 72;

    const state = {
      // cs = 1 sets the unit of speed: every velocity in here is a Mach number.
      cs: 1.0,
      mach: opts.mach == null ? 4.0 : opts.mach,   // target rms Mach
      zeta: opts.zeta == null ? 0.35 : opts.zeta,  // 0 solenoidal, 1 compressive
      cfl: opts.cfl == null ? 0.4 : opts.cfl,
      mode: 0,
      amp: 6.0,
      running: true,
      // measured
      t: 0, steps: 0, dt: 0, maxSig: 1, machRms: 0, sigma: 0, meanLn: 0,
      mass: 1, mass0: null, massErr: 0, sigmaBar: 0, machBar: 0,
      compress: 0, rhoMax: 1, fps: 60,
      pdf: null, pdfLo: -3, pdfHi: 3
    };

    let U, met, red = [], pdfTex, grid = { w: 0, h: 0 }, dpr = 1;

    // Eight modes at the box scale. k_x is scaled by round(aspect) so the
    // wavenumbers stay integer -- hence exactly periodic -- while being close to
    // isotropic in physical units, since the box is aspect-ratio wide and 1 tall.
    const KBASE = [[1, 0], [0, 1], [1, 1], [1, -1], [2, 0], [0, 2], [2, 1], [1, 2]];
    const kArr = new Float32Array(NMODE * 2);
    const aArr = new Float32Array(NMODE * 2);
    for (let i = 0; i < NMODE * 2; i++) aArr[i] = 0;

    function buildModes() {
      const kx = Math.max(1, Math.round(grid.w / grid.h));
      for (let i = 0; i < NMODE; i++) {
        kArr[i * 2 + 0] = KBASE[i][0] * kx;
        kArr[i * 2 + 1] = KBASE[i][1];
      }
    }

    // Box-Muller, one cached spare
    let spare = null;
    function gauss() {
      if (spare !== null) { const v = spare; spare = null; return v; }
      let u = 0, v = 0;
      while (u === 0) u = Math.random();
      v = Math.random();
      const r = Math.sqrt(-2 * Math.log(u));
      spare = r * Math.sin(2 * Math.PI * v);
      return r * Math.cos(2 * Math.PI * v);
    }

    // da = -a dt/tc + sqrt(2 dt/tc) dW, so <|a|^2> -> 1 with correlation time tc.
    // tc is the box eddy-turnover time: L / (M cs), with L = 1.
    function stirOU(dt) {
      const tc = 1 / Math.max(state.machRms, 0.5);
      const f = Math.min(1, dt / tc);
      const g = Math.sqrt(2 * f);
      for (let i = 0; i < NMODE * 2; i++) aArr[i] = aArr[i] * (1 - f) + g * gauss();
    }
    let owned = [];

    function drawTo(f) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, f ? f.fbo : null);
      gl.viewport(0, 0, f ? f.w : canvas.width, f ? f.h : canvas.height);
      gl.bindVertexArray(vao);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function release() {
      for (const o of owned) { gl.deleteFramebuffer(o.fbo); gl.deleteTexture(o.tex); }
      owned = [];
    }

    function mk(w, h, filter) {
      const o = makeFBO(gl, w, h, filter || gl.NEAREST);
      owned.push(o);
      return o;
    }

    function allocate(n) {
      release();
      dpr = Math.min(global.devicePixelRatio || 1, 1.5);
      canvas.width = Math.max(1, Math.round(canvas.clientWidth * dpr));
      canvas.height = Math.max(1, Math.round(canvas.clientHeight * dpr));
      const aspect = canvas.width / Math.max(1, canvas.height);
      grid.h = n;
      grid.w = Math.max(32, Math.round(n * aspect / 4) * 4);

      U = { a: mk(grid.w, grid.h), b: mk(grid.w, grid.h) };
      met = mk(grid.w, grid.h);
      red = [];
      let w = grid.w, h = grid.h;
      while (w > 4 || h > 4) {
        w = Math.max(1, Math.ceil(w / 4));
        h = Math.max(1, Math.ceil(h / 4));
        red.push(mk(w, h));
      }
      pdfTex = mk(PDF_W, PDF_H);
      buildModes();
      reset();
    }

    function reset() {
      gl.useProgram(P.init.p);
      gl.uniform1f(P.init.u.seed, Math.floor(Math.random() * 65535));
      drawTo(U.a);
      drawTo(U.b);
      state.t = 0; state.steps = 0; state.mass0 = null; state.massErr = 0;
      state.sigmaBar = 0; state.machBar = 0;
      measure();

      // Sanity: the initial condition is a gas of density ~1. If the state
      // reads back as zero, the render target is unusable -- an unfilterable
      // format asked to filter, a driver quirk, anything -- and the honest
      // response is to fail visibly rather than present a still, empty box as
      // though it were physics.
      const probe = new Float32Array(4);
      gl.bindFramebuffer(gl.FRAMEBUFFER, U.a.fbo);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, probe);
      if (!(probe[0] > 0.1)) {
        throw new Error('state texture reads back as ' + probe[0] + ' after init');
      }
    }

    function swap() { const s = U.a; U.a = U.b; U.b = s; }

    // ------------------------------------------------------------- one step

    function sweep(dx, dy, dt) {
      gl.useProgram(P.sweep.p);
      gl.uniform1i(P.sweep.u.uU, U.a.bind(0));
      gl.uniform2i(P.sweep.u.size, grid.w, grid.h);
      gl.uniform2i(P.sweep.u.dir, dx, dy);
      // dx here is the cell size along the sweep: the grid is square-celled, so
      // the shorter dimension sets it and the domain is 1 unit tall.
      gl.uniform1f(P.sweep.u.dtdx, dt * grid.h);
      gl.uniform1f(P.sweep.u.cs2, state.cs * state.cs);
      gl.uniform1f(P.sweep.u.rhoMin, RHO_MIN);
      drawTo(U.b);
      swap();
    }

    function step() {
      // dt from the measured signal speed, with a margin because the maximum
      // can grow between measurements.
      // The margin covers growth in the signal speed between measurements: a
      // shock forming mid-interval can raise max(|u| + cs) sharply, and a stale
      // maximum means a violated CFL condition and a dead solver.
      const dt = state.cfl / (grid.h * Math.max(state.maxSig, state.cs) * 1.6);
      state.dt = dt;

      // Driving amplitude is servoed toward the target Mach number instead of
      // being set open-loop: the dissipation rate of a shock-dominated flow is
      // not something you can predict well enough to hit a Mach number by hand.
      stirOU(dt);
      gl.useProgram(P.drive.p);
      gl.uniform1i(P.drive.u.uU, U.a.bind(0));
      gl.uniform1f(P.drive.u.amp, state.amp);
      gl.uniform1f(P.drive.u.zeta, state.zeta);
      gl.uniform1f(P.drive.u.dt, dt);
      gl.uniform2fv(P.drive.u.kvec, kArr);
      gl.uniform2fv(P.drive.u.amps, aArr);
      drawTo(U.b);
      swap();

      // Strang alternation: xy on even steps, yx on odd, so the splitting
      // error stays second order instead of accumulating a preferred axis.
      if (state.steps & 1) { sweep(0, 1, dt); sweep(1, 0, dt); }
      else                 { sweep(1, 0, dt); sweep(0, 1, dt); }

      state.t += dt;
      state.steps++;
    }

    // ----------------------------------------------------------- diagnostics

    const buf = new Float32Array(4 * 64);
    const pdfBuf = new Float32Array(4 * PDF_W * PDF_H);

    function measure() {
      gl.useProgram(P.metric.p);
      gl.uniform1i(P.metric.u.uU, U.a.bind(0));
      gl.uniform1f(P.metric.u.cs2, state.cs * state.cs);
      gl.uniform1f(P.metric.u.rhoMin, RHO_MIN);
      drawTo(met);

      let src = met;
      for (const d of red) {
        gl.useProgram(P.reduce.p);
        gl.uniform1i(P.reduce.u.uSrc, src.bind(0));
        gl.uniform2i(P.reduce.u.srcSize, src.w, src.h);
        drawTo(d);
        src = d;
      }

      gl.bindFramebuffer(gl.FRAMEBUFFER, src.fbo);
      const n = Math.min(64, src.w * src.h);
      try { gl.readPixels(0, 0, src.w, src.h, gl.RGBA, gl.FLOAT, buf); }
      catch (e) { return; }

      let mx = 0, sq = 0, sm = 0, sl = 0;
      for (let i = 0; i < n; i++) {
        mx = Math.max(mx, buf[i * 4 + 0]);
        sq += buf[i * 4 + 1];
        sm += buf[i * 4 + 2];
        sl += buf[i * 4 + 3];
      }
      const cells = grid.w * grid.h;
      if (isFinite(mx) && mx > 0) state.maxSig = mx;
      if (isFinite(sq)) state.machRms = Math.sqrt(Math.max(0, sq / cells)) / state.cs;
      if (isFinite(sl)) state.meanLn = sl / cells;
      state.machMax = Math.max(0, (state.maxSig - state.cs) / state.cs);
      // mean density is 1 by construction of the initial condition, so the
      // relative drift of the mean is the relative mass error outright
      if (isFinite(sm) && sm > 0) {
        state.mass = sm / cells;
        if (state.mass0 == null) state.mass0 = state.mass;
        state.massErr = Math.abs(state.mass / state.mass0 - 1);
      }

      // Drive amplitude servo. Runs once per measurement -- about ten times a
      // second -- with a fractional error and a hard slew limit, so it settles
      // in a couple of seconds and cannot wind up between measurements.
      const rel = (state.mach - state.machRms) / Math.max(state.mach, 0.5);
      const slew = Math.max(-0.35, Math.min(0.35, rel));
      state.amp = Math.max(0.02, Math.min(500, state.amp * (1 + 0.45 * slew)));

      // Last resort. If the flow has run away or gone non-finite there is
      // nothing to show and nothing to learn from leaving it up; start over
      // rather than present a dead box as a simulation.
      if (!isFinite(state.machRms) || state.maxSig > 400 * state.cs) {
        state.amp = 1.0;
        reset();
      }
    }

    // rms Mach needs a sum of |u|^2, which the max-channel cannot give. Take it
    // from the decimated sample, which is already coming back for the PDF.
    function samplePDF() {
      gl.useProgram(P.decim.p);
      gl.uniform1i(P.decim.u.uU, U.a.bind(0));
      gl.uniform2i(P.decim.u.stride, Math.max(1, Math.floor(grid.w / PDF_W)),
                                     Math.max(1, Math.floor(grid.h / PDF_H)));
      gl.uniform2i(P.decim.u.size, grid.w, grid.h);
      gl.uniform1f(P.decim.u.rhoMin, RHO_MIN);
      gl.uniform1f(P.decim.u.cs2, state.cs * state.cs);
      drawTo(pdfTex);
      gl.bindFramebuffer(gl.FRAMEBUFFER, pdfTex.fbo);
      try { gl.readPixels(0, 0, PDF_W, PDF_H, gl.RGBA, gl.FLOAT, pdfBuf); }
      catch (e) { return; }

      const N = PDF_W * PDF_H;
      const NB = 48;
      // window from the previous call's sigma; it changes far more slowly than
      // the interval between calls
      const lo = state.meanLn - 4.2 * Math.max(state.sigma, 0.05);
      const hi = state.meanLn + 4.2 * Math.max(state.sigma, 0.05);
      const bins = new Float32Array(NB);
      let conv = 0, used = 0, sl = 0, sl2 = 0;
      // a cell counts as shocked when it converges faster than a fifth of the
      // sound speed per cell: strong enough that only a real compression
      // qualifies, loose enough to catch the weak ones
      const thresh = -0.2 * state.cs;
      for (let i = 0; i < N; i++) {
        const l = pdfBuf[i * 4];
        if (!isFinite(l)) continue;
        used++;
        sl += l; sl2 += l * l;
        if (pdfBuf[i * 4 + 2] < thresh) conv++;
        let k = Math.floor((l - lo) / (hi - lo) * NB);
        if (k < 0) k = 0; else if (k >= NB) k = NB - 1;
        bins[k] += 1;
      }
      state.compress = used ? conv / used : 0;
      if (used > 64) {
        const m = sl / used;
        state.sigma = Math.sqrt(Math.max(0, sl2 / used - m * m));
        // A single snapshot of a 2-D box scatters by tens of per cent: sigma is
        // a second moment of a heavy-tailed field, and three box crossings is
        // nowhere near enough samples. Carry an exponential running mean with a
        // ~5 s time constant so the printed number means something.
        const A = 0.02;
        state.sigmaBar = state.sigmaBar > 0 ? state.sigmaBar + (state.sigma - state.sigmaBar) * A : state.sigma;
        state.machBar = state.machBar > 0 ? state.machBar + (state.machRms - state.machBar) * A : state.machRms;
      }
      const width = (hi - lo) / NB;
      for (let k = 0; k < NB; k++) bins[k] /= N * width;   // normalised density
      state.pdf = bins;
      state.pdfLo = lo;
      state.pdfHi = hi;
    }

    // ---------------------------------------------------------------- render

    const PAL = {
      dye: [[0.02, 0.03, 0.09], [0.10, 0.16, 0.42], [0.42, 0.55, 0.92], [0.90, 0.95, 1.00]],
      shock: [1.00, 0.62, 0.24]
    };

    function render() {
      gl.useProgram(P.show.p);
      gl.uniform1i(P.show.u.uU, U.a.bind(0));
      gl.uniform2i(P.show.u.size, grid.w, grid.h);
      gl.uniform2f(P.show.u.uRes, canvas.width, canvas.height);
      const c = PAL.dye;
      gl.uniform3fv(P.show.u.uC0, c[0]);
      gl.uniform3fv(P.show.u.uC1, c[1]);
      gl.uniform3fv(P.show.u.uC2, c[2]);
      gl.uniform3fv(P.show.u.uC3, c[3]);
      gl.uniform3fv(P.show.u.uShock, PAL.shock);
      gl.uniform1f(P.show.u.cs2, state.cs * state.cs);
      gl.uniform1f(P.show.u.rhoMin, RHO_MIN);
      // stretch follows the measured distribution, so the image stays readable
      // as the Mach number and therefore the contrast change
      const s = Math.max(state.sigma, 0.12);
      gl.uniform1f(P.show.u.logLo, state.meanLn - 2.4 * s);
      gl.uniform1f(P.show.u.logHi, state.meanLn + 2.6 * s);
      // div is du per cell and a shock is a jump of order cs across a cell, so
      // cs is the natural unit for both the threshold and the gain.
      gl.uniform1f(P.show.u.shockGain, 1.0 / (0.45 * state.cs));
      gl.uniform1f(P.show.u.shockMin, 0.12 * state.cs);
      gl.uniform1f(P.show.u.mode, state.mode);
      gl.uniform1f(P.show.u.vignette, 0.55);
      drawTo(null);
    }

    // ------------------------------------------------------------------ loop

    let raf = 0, last = 0, frames = 0, fpsT = 0, mCount = 0;
    let STEPS_PER_FRAME = opts.substeps || 3;

    function frame(now) {
      raf = requestAnimationFrame(frame);
      if (!last) { last = now; return; }
      const el = Math.min((now - last) / 1000, 0.1);
      last = now;
      frames++; fpsT += el;
      if (fpsT > 0.5) { state.fps = frames / fpsT; frames = 0; fpsT = 0; }

      if (state.running) {
        for (let i = 0; i < STEPS_PER_FRAME; i++) step();
        // readPixels syncs the pipeline, so diagnostics are periodic
        if (--mCount <= 0) {
          mCount = 4;
          measure();
          samplePDF();
          if (state.onstats) state.onstats(stats());
        }
      }
      render();
    }

    function stats() {
      // sigma is the headline number: for isothermal supersonic turbulence the
      // density PDF is log-normal with sigma^2 = ln(1 + b^2 M^2), and b runs
      // from about 1/3 for purely solenoidal driving to 1 for purely
      // compressive. Solving that for b is a live test of the solver.
      // b from the running means, not from one frame, or the number jitters by
      // a factor and means nothing.
      const M = Math.max(state.machBar || state.machRms, 1e-3);
      const sg = state.sigmaBar || state.sigma;
      const bFit = Math.sqrt(Math.max(0, Math.exp(sg * sg) - 1)) / M;
      return {
        grid: grid.w + ' × ' + grid.h, gridW: grid.w, gridH: grid.h,
        t: state.t, steps: state.steps, dt: state.dt, cfl: state.cfl,
        cs: state.cs, machTarget: state.mach, machMax: state.machMax || 0,
        machRms: state.machRms, zeta: state.zeta, sigma: state.sigma,
        sigmaBar: state.sigmaBar, machBar: state.machBar,
        bFit: bFit, bPred: 1 / 3 + (1 - 1 / 3) * state.zeta,
        sigmaPred: Math.sqrt(Math.log(1 + Math.pow((1 / 3 + (1 - 1 / 3) * state.zeta) * M, 2))),
        compress: state.compress, fps: state.fps, amp: state.amp,
        mass: state.mass, massErr: state.massErr,
        pdf: state.pdf, pdfLo: state.pdfLo, pdfHi: state.pdfHi,
        meanLn: state.meanLn, mode: state.mode
      };
    }

    allocate(opts.n || 288);

    let resizeT = null;
    global.addEventListener('resize', function () {
      clearTimeout(resizeT);
      resizeT = setTimeout(function () { allocate(grid.h); }, 250);
    });

    raf = requestAnimationFrame(frame);

    return {
      kind: 'webgl2',
      stats,
      onStats(fn) { state.onstats = fn; },
      setMach(m) { state.mach = Math.max(0.2, Math.min(20, m)); },
      setZeta(z) { state.zeta = Math.max(0, Math.min(1, z)); },
      setMode(m) { state.mode = m; },
      setResolution(n) { allocate(Math.max(64, Math.min(512, n))); },
      pause() { state.running = false; },
      play() { state.running = true; },
      toggle() { state.running = !state.running; return state.running; },
      setSpeed(n) { STEPS_PER_FRAME = Math.max(1, Math.min(8, n | 0)); },
      reset,
      blast(x, y, dens, kick, radius) {
        gl.useProgram(P.blast.p);
        gl.uniform1i(P.blast.u.uU, U.a.bind(0));
        gl.uniform2f(P.blast.u.point, x, y);
        gl.uniform1f(P.blast.u.aspect, grid.w / grid.h);
        gl.uniform1f(P.blast.u.radius, radius == null ? 0.0009 : radius);
        gl.uniform1f(P.blast.u.dens, dens == null ? 6.0 : dens);
        gl.uniform1f(P.blast.u.kick, kick == null ? 4.0 : kick);
        drawTo(U.b);
        swap();
      },
      push(x, y, dx, dy, radius) {
        gl.useProgram(P.push.p);
        gl.uniform1i(P.push.u.uU, U.a.bind(0));
        gl.uniform2f(P.push.u.point, x, y);
        gl.uniform2f(P.push.u.delta, dx, dy);
        gl.uniform1f(P.push.u.aspect, grid.w / grid.h);
        gl.uniform1f(P.push.u.radius, radius == null ? 0.0016 : radius);
        drawTo(U.b);
        swap();
      },
      destroy() { cancelAnimationFrame(raf); release(); }
    };
  }

  global.Euler = {
    create(canvas, opts) {
      try { return Euler(canvas, opts); }
      catch (e) {
        if (global.console) console.warn('[euler] ' + e.message);
        return null;
      }
    }
  };
})(window);

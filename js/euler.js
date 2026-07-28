/* ============================================================================
   euler.js — compressible, for once.

   The background on the main page is incompressible: a projection method, an
   elliptic pressure, no shocks anywhere by construction. This is the other half
   of the subject.

   The scheme is the one from mini-ramses-ism's GPU hydro path (gpu_hydro.cuf),
   cut down to two dimensions, to an isothermal equation of state, and to the
   cheapest configuration that still captures a shock:

     variables   : conserved (rho, rho u, rho v). Isothermal closure P = rho cs^2
                   with cs constant, so there is no energy equation and no sqrt
                   anywhere in the solver -- the sound speed is a uniform.
     reconstruct : none. Piecewise constant, so the interface states are the cell
                   averages: first-order Godunov, no slopes, no limiter, no
                   predictor half-step.
     flux        : HLL, with RAMSES's wave-speed estimate
                     SL = min(min(uL,uR) - cs, 0)
                     SR = max(max(uL,uR) + cs, 0)
                   Clamping both speeds through zero is what makes this cheap:
                   SL <= 0 <= SR always, so the central formula is always the
                   right one, the supersonic branches disappear, and SR - SL >=
                   2 cs means the division needs no guard. LLF (Rusanov) is
                   selectable and cheaper still -- one speed, used symmetrically.
     update      : unsplit and conservative,
                     U += dt/dx [ (Fx_i - Fx_{i+1}) + (Fy_j - Fy_{j+1}) ]
     timestep    : RAMSES's, dt = C dx / (|u| + |v| + ndim cs), the directional
                   sum rather than max(|u| + cs) because the update is unsplit,
                   with the maximum measured on the GPU rather than assumed.
     floors      : smallr on density, where the reference applies it. The
                   pressure floor the reference needs is automatic here, since
                   P = rho cs^2 cannot go negative once rho cannot.

   Isothermal is the right closure for this and not a shortcut. It is standard for
   the cold neutral and molecular ISM, where line cooling pins the temperature
   over a dynamical time; it keeps a driven box supersonic indefinitely, where an
   adiabatic one would shock-heat until there was nothing left to look at; and it
   makes the sound speed a constant, which removes the only sqrt in the hot path.

   Precision is single throughout, which is less a choice than the only option:
   WebGL2's highp float is IEEE binary32 and there is no double path at all. The
   state textures are RGBA32F to match, with three of the four channels used.

   GLSL ES has no -ffast-math switch, so the fast-math intent is written into the
   code instead: reciprocal multiplies rather than divides (as the reference does
   via rinv), no sqrt, no pow, no exp in the hot path, and a branch-free Riemann
   solver thanks to the clamped wave speeds.

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

  // Shared equation-of-state helpers, mirroring the reference's
  // conserved_2_primitive / primitive_2_conserved. Primitives are carried as
  // (rho, u_normal, u_transverse) so one flux expression serves both directions --
  // the same trick the reference uses by always putting the normal velocity in
  // velocity_x and rotating around it.
  const EOS = `
uniform float cs, cs2, smallr;

vec3 toPrim(vec3 U){
  float r  = max(U.x, smallr);
  float ri = 1.0 / r;                       // one divide, then multiplies
  return vec3(r, U.y * ri, U.z * ri);
}
vec3 toCons(vec3 q){ return vec3(q.x, q.x * q.y, q.x * q.z); }
`;

  // The whole solver: one pass. Piecewise-constant states, four Riemann problems,
  // one unsplit conservative update.
  //
  // Each face is solved twice, once for each of the cells that share it. The
  // reference computes a face once and stores it, which is right on a CPU; here a
  // stored flux would cost a second full-screen pass, and on this GPU a pass
  // costs far more than the arithmetic it saves. Four Riemann solves in one pass
  // beat two solves in two passes.
  const F_GODUNOV = HEAD + EOS + `
uniform sampler2D uU;
uniform ivec2 size;
uniform float dtdx, llf;

vec3 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0).xyz;
}

// HLL flux for a pair of primitive states whose .y component is the velocity
// along the face normal. Returns the conserved flux in the same convention.
// No sqrt: the sound speed is a constant in an isothermal gas.
vec3 hll(vec3 qL, vec3 qR){
  float SL, SR;
  if (llf > 0.5) {
    // LLF / Rusanov: one speed, used both ways. Two fewer min/max than HLL, and
    // more diffusive by about the same margin.
    float sm = max(abs(qL.y), abs(qR.y)) + cs;
    SL = -sm; SR = sm;
  } else {
    // Clamped through zero, as in the reference. SR - SL >= 2 cs > 0, so the
    // division below can never be by zero and needs no guard.
    SL = min(min(qL.y, qR.y) - cs, 0.0);
    SR = max(max(qL.y, qR.y) + cs, 0.0);
  }

  vec3 UL = toCons(qL), UR = toCons(qR);
  vec3 FL = vec3(UL.y, qL.y * UL.y + qL.x * cs2, qL.y * UL.z);
  vec3 FR = vec3(UR.y, qR.y * UR.y + qR.x * cs2, qR.y * UR.z);

  return (SR * FL - SL * FR + SR * SL * (UR - UL)) / (SR - SL);
}

// swap the velocity components, so a y-face reuses the x-face solver
vec3 swapv(vec3 a){ return vec3(a.x, a.z, a.y); }

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);

  vec3 U0 = cell(c);
  vec3 q0 = toPrim(U0);
  vec3 qw = toPrim(cell(c - ivec2(1, 0)));
  vec3 qe = toPrim(cell(c + ivec2(1, 0)));
  vec3 qs = toPrim(cell(c - ivec2(0, 1)));
  vec3 qn = toPrim(cell(c + ivec2(0, 1)));

  vec3 Fw = hll(qw, q0);                          // face at i-1/2
  vec3 Fe = hll(q0, qe);                          // face at i+1/2
  vec3 Fs = swapv(hll(swapv(qs), swapv(q0)));     // face at j-1/2
  vec3 Fn = swapv(hll(swapv(q0), swapv(qn)));     // face at j+1/2

  vec3 U = U0 + ((Fw - Fe) + (Fs - Fn)) * dtdx;
  U.x = max(U.x, smallr);
  outColor = vec4(U, 1.0);
}`;

  // Large-scale driving, as an Ornstein-Uhlenbeck process.
  //
  // Each mode carries a complex amplitude driven by
  //     da = -a dt/tc + sqrt(2 dt/tc) dW
  // with tc the box eddy-turnover time, which gives unit variance in the steady
  // state and a finite correlation time. An earlier version used steadily
  // rotating phases, which is coherent forcing: it compresses the same places
  // over and over and produced a density contrast far larger than the physics
  // warrants. The amplitudes live on the CPU -- eight modes is sixteen numbers --
  // and arrive as a uniform array.
  //
  // Each mode splits into a solenoidal part (the curl of a stream function,
  // divergence-free) and a compressive part (the gradient of a potential,
  // curl-free). zeta mixes them: 0 purely solenoidal, 1 purely compressive. That
  // is the knob that moves b in sigma^2 = ln(1 + b^2 M^2).
  const NMODE = 8;
  const F_DRIVE = HEAD + `
uniform sampler2D uU;
uniform float amp, zeta, dt, smallr;
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
    float S  = amps[i].x * sin(th) + amps[i].y * cos(th);
    vec2  kh = k / kk;
    sol  += S * vec2(-kh.y, kh.x);
    comp += S * -kh;
  }
  // acceleration, so the force on a cell is rho * a: heavy gas is harder to
  // push, which is part of why compressive driving widens the density PDF
  vec2 a = mix(sol, comp, zeta) * amp * (1.0 / sqrt(8.0));
  outColor = vec4(U.x, U.yz + max(U.x, smallr) * a * dt, 1.0);
}`;

  // A blast. In an isothermal gas the pressure is rho cs^2, so piling up density
  // *is* piling up pressure: this needs no energy injection and no painted-on
  // velocity ring. The Riemann solver turns the overpressure into an outward
  // shock by itself, which is the point of having one.
  const F_BLAST = HEAD + `
uniform sampler2D uU;
uniform vec2  point;
uniform float aspect, radius, dens, smallr;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float w = exp(-dot(p, p) / radius);
  float r0 = max(U.x, smallr);
  float r1 = U.x + dens * w;
  // the new mass arrives with the local velocity, so the blast is a pressure
  // perturbation and not also a momentum one
  outColor = vec4(r1, U.yz * (r1 / r0), 1.0);
}`;

  // A stir: momentum along a stroke.
  const F_PUSH = HEAD + `
uniform sampler2D uU;
uniform vec2  point, delta;
uniform float aspect, radius, smallr;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float w = exp(-dot(p, p) / radius);
  outColor = vec4(U.x, U.yz + delta * w * max(U.x, smallr), 1.0);
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
void main(){
  // uniform gas at rest, with a whisper of density noise so the driving has
  // something to act on rather than an exactly symmetric state
  outColor = vec4(1.0 + 0.02 * (h1(uvec2(gl_FragCoord.xy), 17u) - 0.5), 0.0, 0.0, 1.0);
}`;

  // Per-cell diagnostics, gathered so one readback answers everything.
  //   R: |u| + |v| + 2 cs   reduced by MAX -> the RAMSES unsplit CFL bound
  //   G: |u|^2 / cs^2       reduced by sum -> rms Mach number
  //   B: rho                reduced by sum -> total mass, which a conservative
  //                         scheme must hold to round-off
  //   A: ln rho             reduced by sum -> mean of ln rho
  const F_METRICS = HEAD + EOS + `
uniform sampler2D uU;
void main(){
  vec3  q = toPrim(texelFetch(uU, ivec2(gl_FragCoord.xy), 0).xyz);
  outColor = vec4(abs(q.y) + abs(q.z) + 2.0 * cs,
                  dot(q.yz, q.yz) / cs2, q.x, log(q.x));
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
  // Carries divergence and Mach too, so one readback yields the PDF, the shocked
  // fraction and the peak Mach number.
  const F_DECIMATE = HEAD + EOS + `
uniform sampler2D uU;
uniform ivec2 stride, size;
vec3 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0).xyz;
}
void main(){
  ivec2 p = ivec2(gl_FragCoord.xy) * stride;
  vec3 q  = toPrim(cell(p));
  vec3 qe = toPrim(cell(p + ivec2(1, 0))), qw = toPrim(cell(p - ivec2(1, 0)));
  vec3 qn = toPrim(cell(p + ivec2(0, 1))), qs = toPrim(cell(p - ivec2(0, 1)));
  // du per cell. A shock is a jump of order cs across about one cell, so this is
  // already in the natural units for deciding what counts as one.
  float div = 0.5 * ((qe.y - qw.y) + (qn.z - qs.z));
  outColor = vec4(log(q.x), q.x, div, length(q.yz) / cs);
}`;

  // Display. Density on a log stretch is the base image, because that is what
  // the physics is log-normal in; convergence is overlaid, because that is where
  // the shocks are.
  const F_SHOW = HEAD + EOS + `
uniform sampler2D uU;
uniform ivec2 size;
uniform vec2  uRes;
uniform vec3  uC0, uC1, uC2, uC3, uShock;
uniform float logLo, logHi, shockGain, shockMin, mode, vignette;

vec3 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0).xyz;
}
// Bilinear by hand. The grid is coarse next to the canvas and point-sampling it
// turns every shock into a staircase -- but RGBA32F is not a filterable format
// unless OES_texture_float_linear is present, and setting LINEAR on a texture
// whose format cannot filter makes it INCOMPLETE, after which every read of it
// returns zero, texelFetch included.
vec3 at(vec2 uv){
  vec2  g  = uv * vec2(size) - 0.5;
  ivec2 i0 = ivec2(floor(g));
  vec2  t  = g - vec2(i0);
  return mix(mix(cell(i0),               cell(i0 + ivec2(1, 0)), t.x),
             mix(cell(i0 + ivec2(0, 1)), cell(i0 + ivec2(1, 1)), t.x), t.y);
}
vec3 ramp(float t){
  t = clamp(t, 0.0, 1.0);
  if (t < 0.3333) return mix(uC0, uC1, t / 0.3333);
  if (t < 0.6666) return mix(uC1, uC2, (t - 0.3333) / 0.3333);
  return mix(uC2, uC3, (t - 0.6666) / 0.3334);
}
void main(){
  vec3  q  = toPrim(at(vUv));
  vec2  tx = 1.0 / vec2(size);
  vec3  qe = toPrim(at(vUv + vec2(tx.x, 0.0))), qw = toPrim(at(vUv - vec2(tx.x, 0.0)));
  vec3  qn = toPrim(at(vUv + vec2(0.0, tx.y))), qs = toPrim(at(vUv - vec2(0.0, tx.y)));
  float div = 0.5 * ((qe.y - qw.y) + (qn.z - qs.z));

  float ln   = log(q.x);
  float dn   = (ln - logLo) / max(logHi - logLo, 1e-4);
  float mach = length(q.yz) / cs;

  vec3 col;
  if (mode < 0.5) {
    col = ramp(dn);
    float sh = clamp((-div - shockMin) * shockGain, 0.0, 1.0);
    col += uShock * sh * sh * 0.85;
  } else if (mode < 1.5) {
    col = uShock * clamp((-div - shockMin) * shockGain, 0.0, 1.4);
  } else if (mode < 2.5) {
    col = ramp(clamp(mach / 6.0, 0.0, 1.0));
  } else {
    // divergence, signed: compression warm, expansion cool
    float sg = clamp(div * shockGain * 0.5, -1.0, 1.0);
    col = mix(uC1 * 0.6, uShock, 0.5 + 0.5 * (-sg));
    col *= abs(sg) * 0.9 + 0.05;
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
      const nm = gl.getActiveUniform(p, i).name.replace(/\[0\]$/, '');
      u[nm] = gl.getUniformLocation(p, nm);
    }
    return { p, u };
  }

  function makeFBO(gl, w, h) {
    const tex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
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
      god:    program(gl, VERT, F_GODUNOV),
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

    const SMALLR = 1e-4;
    const PDF_W = 128, PDF_H = 72;

    const state = {
      // cs = 1 sets the unit of speed, so every velocity in here is a Mach number
      cs: 1.0,
      mach: opts.mach == null ? 4.0 : opts.mach,   // target rms Mach
      zeta: opts.zeta == null ? 0.35 : opts.zeta,  // 0 solenoidal, 1 compressive
      cfl: opts.cfl == null ? 0.5 : opts.cfl,      // RAMSES courant_factor
      llf: 0,
      mode: 0,
      amp: 6.0,
      running: true,
      t: 0, steps: 0, dt: 0, maxSig: 1, machRms: 0, sigma: 0, meanLn: 0,
      machMax: 0, mass: 1, mass0: null, massErr: 0,
      sigmaBar: 0, machBar: 0, compress: 0, fps: 60,
      pdf: null, pdfLo: -3, pdfHi: 3
    };

    let U, met, red = [], pdfTex, grid = { w: 0, h: 0 }, dpr = 1;
    let owned = [];

    // Eight modes at the box scale. k_x is scaled by round(aspect) so the
    // wavenumbers stay integer -- hence exactly periodic -- while being close to
    // isotropic in physical units, since the box is aspect-ratio wide and 1 tall.
    const KBASE = [[1, 0], [0, 1], [1, 1], [1, -1], [2, 0], [0, 2], [2, 1], [1, 2]];
    const kArr = new Float32Array(NMODE * 2);
    const aArr = new Float32Array(NMODE * 2);

    function buildModes() {
      const kx = Math.max(1, Math.round(grid.w / grid.h));
      for (let i = 0; i < NMODE; i++) {
        kArr[i * 2 + 0] = KBASE[i][0] * kx;
        kArr[i * 2 + 1] = KBASE[i][1];
      }
    }

    let spare = null;
    function gauss() {
      if (spare !== null) { const v = spare; spare = null; return v; }
      let u = 0;
      while (u === 0) u = Math.random();
      const r = Math.sqrt(-2 * Math.log(u)), th = 2 * Math.PI * Math.random();
      spare = r * Math.sin(th);
      return r * Math.cos(th);
    }

    // da = -a dt/tc + sqrt(2 dt/tc) dW, so <|a|^2> -> 1 with correlation time tc.
    function stirOU(dt) {
      const tc = 1 / Math.max(state.machRms, 0.5);
      const f = Math.min(1, dt / tc), g = Math.sqrt(2 * f);
      for (let i = 0; i < NMODE * 2; i++) aArr[i] = aArr[i] * (1 - f) + g * gauss();
    }

    function drawTo(f) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, f ? f.fbo : null);
      gl.viewport(0, 0, f ? f.w : canvas.width, f ? f.h : canvas.height);
      gl.bindVertexArray(vao);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function mk(w, h) {
      const o = makeFBO(gl, w, h);
      owned.push(o);
      return o;
    }

    function release() {
      for (const o of owned) { gl.deleteFramebuffer(o.fbo); gl.deleteTexture(o.tex); }
      owned = [];
    }

    // Every program that touches the state needs the same equation of state.
    function eos(pr) {
      gl.uniform1f(pr.u.cs, state.cs);
      gl.uniform1f(pr.u.cs2, state.cs * state.cs);
      gl.uniform1f(pr.u.smallr, SMALLR);
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

    function swap() { const s = U.a; U.a = U.b; U.b = s; }

    function reset() {
      gl.useProgram(P.init.p);
      drawTo(U.a);
      drawTo(U.b);
      state.t = 0; state.steps = 0; state.mass0 = null; state.massErr = 0;
      state.sigmaBar = 0; state.machBar = 0;
      measure();

      // Sanity: the initial condition is a gas of density ~1. If the state reads
      // back as zero the render target is unusable, and the honest response is to
      // fail visibly rather than present a still, empty box as though it were
      // physics.
      const probe = new Float32Array(4);
      gl.bindFramebuffer(gl.FRAMEBUFFER, U.a.fbo);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, probe);
      if (!(probe[0] > 0.1)) {
        throw new Error('state texture reads back as ' + probe[0] + ' after init');
      }
    }

    // ------------------------------------------------------------- one step

    function step() {
      // RAMSES's timestep: ctot sums the directional signal speeds because the
      // update is unsplit, which is stricter than max(|u| + c) and is why the
      // courant factor can sit at 0.5. The margin covers growth in ctot between
      // measurements -- a shock forming mid-interval raises it sharply, and a
      // stale maximum means a violated CFL condition and a dead solver.
      const dt = state.cfl / (grid.h * Math.max(state.maxSig, state.cs) * 1.6);
      state.dt = dt;

      stirOU(dt);
      gl.useProgram(P.drive.p);
      gl.uniform1i(P.drive.u.uU, U.a.bind(0));
      gl.uniform1f(P.drive.u.amp, state.amp);
      gl.uniform1f(P.drive.u.zeta, state.zeta);
      gl.uniform1f(P.drive.u.dt, dt);
      gl.uniform1f(P.drive.u.smallr, SMALLR);
      gl.uniform2fv(P.drive.u.kvec, kArr);
      gl.uniform2fv(P.drive.u.amps, aArr);
      drawTo(U.b);
      swap();

      // one unsplit Godunov pass
      gl.useProgram(P.god.p);
      eos(P.god);
      gl.uniform1i(P.god.u.uU, U.a.bind(0));
      gl.uniform2i(P.god.u.size, grid.w, grid.h);
      gl.uniform1f(P.god.u.dtdx, dt * grid.h);   // square cells, box 1 unit tall
      gl.uniform1f(P.god.u.llf, state.llf);
      drawTo(U.b);
      swap();

      state.t += dt;
      state.steps++;
    }

    // ----------------------------------------------------------- diagnostics

    const buf = new Float32Array(4 * 64);
    const pdfBuf = new Float32Array(4 * PDF_W * PDF_H);

    function measure() {
      gl.useProgram(P.metric.p);
      eos(P.metric);
      gl.uniform1i(P.metric.u.uU, U.a.bind(0));
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
      if (isFinite(sq)) state.machRms = Math.sqrt(Math.max(0, sq / cells));
      if (isFinite(sl)) state.meanLn = sl / cells;
      if (isFinite(sm) && sm > 0) {
        state.mass = sm / cells;
        if (state.mass0 == null) state.mass0 = state.mass;
        state.massErr = Math.abs(state.mass / state.mass0 - 1);
      }

      // Drive amplitude servo. Once per measurement -- about ten times a second
      // -- with a fractional error and a hard slew limit, so it settles in a
      // couple of seconds and cannot wind up between measurements. Running it
      // per step is a runaway: at ~500 steps/s a 2% gain compounds to x148 per
      // second and breaks the CFL condition before the next measurement lands.
      const rel = (state.mach - state.machRms) / Math.max(state.mach, 0.5);
      const slew = Math.max(-0.35, Math.min(0.35, rel));
      state.amp = Math.max(0.02, Math.min(500, state.amp * (1 + 0.45 * slew)));

      if (!isFinite(state.machRms) || state.maxSig > 400 * state.cs) {
        state.amp = 1.0;
        reset();
      }
    }

    function samplePDF() {
      gl.useProgram(P.decim.p);
      eos(P.decim);
      gl.uniform1i(P.decim.u.uU, U.a.bind(0));
      gl.uniform2i(P.decim.u.stride, Math.max(1, Math.floor(grid.w / PDF_W)),
                                     Math.max(1, Math.floor(grid.h / PDF_H)));
      gl.uniform2i(P.decim.u.size, grid.w, grid.h);
      drawTo(pdfTex);
      gl.bindFramebuffer(gl.FRAMEBUFFER, pdfTex.fbo);
      try { gl.readPixels(0, 0, PDF_W, PDF_H, gl.RGBA, gl.FLOAT, pdfBuf); }
      catch (e) { return; }

      const N = PDF_W * PDF_H, NB = 48;
      const lo = state.meanLn - 4.2 * Math.max(state.sigma, 0.05);
      const hi = state.meanLn + 4.2 * Math.max(state.sigma, 0.05);
      const bins = new Float32Array(NB);
      let conv = 0, used = 0, sl = 0, sl2 = 0, mmax = 0;
      // a cell counts as shocked when it converges faster than a fifth of the
      // sound speed per cell
      const thresh = -0.2 * state.cs;
      for (let i = 0; i < N; i++) {
        const l = pdfBuf[i * 4];
        if (!isFinite(l)) continue;
        used++;
        sl += l; sl2 += l * l;
        mmax = Math.max(mmax, pdfBuf[i * 4 + 3]);
        if (pdfBuf[i * 4 + 2] < thresh) conv++;
        let k = Math.floor((l - lo) / (hi - lo) * NB);
        if (k < 0) k = 0; else if (k >= NB) k = NB - 1;
        bins[k] += 1;
      }
      const width = (hi - lo) / NB;
      for (let k = 0; k < NB; k++) bins[k] /= N * width;
      state.pdf = bins;
      state.pdfLo = lo;
      state.pdfHi = hi;
      state.compress = used ? conv / used : 0;
      if (used > 64) {
        const m = sl / used;
        state.sigma = Math.sqrt(Math.max(0, sl2 / used - m * m));
        // from the decimated sample, so it is a lower bound on the true peak
        state.machMax = mmax;
        // A single snapshot of a 2-D box scatters by tens of per cent, so carry
        // running means with a ~5 s time constant for the printed numbers.
        const A = 0.02;
        state.sigmaBar = state.sigmaBar > 0 ? state.sigmaBar + (state.sigma - state.sigmaBar) * A : state.sigma;
        state.machBar = state.machBar > 0 ? state.machBar + (state.machRms - state.machBar) * A : state.machRms;
      }
    }

    // ---------------------------------------------------------------- render

    const PAL = {
      dye: [[0.02, 0.03, 0.09], [0.10, 0.16, 0.42], [0.42, 0.55, 0.92], [0.90, 0.95, 1.00]],
      shock: [1.00, 0.62, 0.24]
    };

    function render() {
      gl.useProgram(P.show.p);
      eos(P.show);
      gl.uniform1i(P.show.u.uU, U.a.bind(0));
      gl.uniform2i(P.show.u.size, grid.w, grid.h);
      gl.uniform2f(P.show.u.uRes, canvas.width, canvas.height);
      const c = PAL.dye;
      gl.uniform3fv(P.show.u.uC0, c[0]);
      gl.uniform3fv(P.show.u.uC1, c[1]);
      gl.uniform3fv(P.show.u.uC2, c[2]);
      gl.uniform3fv(P.show.u.uC3, c[3]);
      gl.uniform3fv(P.show.u.uShock, PAL.shock);
      gl.uniform1f(P.show.u.cs, state.cs);
      // stretch follows the measured distribution, so the image stays readable
      // as the Mach number and therefore the contrast change
      const s = Math.max(state.sigma, 0.12);
      gl.uniform1f(P.show.u.logLo, state.meanLn - 2.4 * s);
      gl.uniform1f(P.show.u.logHi, state.meanLn + 2.6 * s);
      gl.uniform1f(P.show.u.shockGain, 1.0 / (0.45 * state.cs));
      gl.uniform1f(P.show.u.shockMin, 0.12 * state.cs);
      gl.uniform1f(P.show.u.mode, state.mode);
      gl.uniform1f(P.show.u.vignette, 0.55);
      drawTo(null);
    }

    // ------------------------------------------------------------------ loop

    let raf = 0, last = 0, frames = 0, fpsT = 0, mCount = 0;
    let STEPS_PER_FRAME = opts.substeps || 4;

    function frame(now) {
      raf = requestAnimationFrame(frame);
      if (!last) { last = now; return; }
      const el = Math.min((now - last) / 1000, 0.1);
      last = now;
      frames++; fpsT += el;
      if (fpsT > 0.5) { state.fps = frames / fpsT; frames = 0; fpsT = 0; }

      if (state.running) {
        for (let i = 0; i < STEPS_PER_FRAME; i++) step();
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
      // sigma^2 = ln(1 + b^2 M^2) with b from the driving mixture. See the page
      // for how far that should be trusted -- the short answer is not far, in a
      // 2-D box over a few crossing times.
      const M = Math.max(state.machBar || state.machRms, 1e-3);
      const sg = state.sigmaBar || state.sigma;
      const bPred = 1 / 3 + (1 - 1 / 3) * state.zeta;
      return {
        grid: grid.w + ' × ' + grid.h, gridW: grid.w, gridH: grid.h,
        t: state.t, steps: state.steps, dt: state.dt, cfl: state.cfl,
        cs: state.cs, llf: state.llf,
        machTarget: state.mach, machMax: state.machMax, machRms: state.machRms,
        machBar: state.machBar, zeta: state.zeta,
        sigma: state.sigma, sigmaBar: state.sigmaBar,
        bFit: Math.sqrt(Math.max(0, Math.exp(sg * sg) - 1)) / M,
        bPred: bPred,
        sigmaPred: Math.sqrt(Math.log(1 + Math.pow(bPred * M, 2))),
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
      setRiemann(kind) { state.llf = kind === 'llf' ? 1 : 0; },
      setMode(m) { state.mode = m; },
      setResolution(n) { allocate(Math.max(64, Math.min(512, n))); },
      setSpeed(n) { STEPS_PER_FRAME = Math.max(1, Math.min(12, n | 0)); },
      pause() { state.running = false; },
      play() { state.running = true; },
      toggle() { state.running = !state.running; return state.running; },
      reset,
      blast(x, y, dens, radius) {
        gl.useProgram(P.blast.p);
        gl.uniform1i(P.blast.u.uU, U.a.bind(0));
        gl.uniform2f(P.blast.u.point, x, y);
        gl.uniform1f(P.blast.u.aspect, grid.w / grid.h);
        gl.uniform1f(P.blast.u.radius, radius == null ? 0.0009 : radius);
        gl.uniform1f(P.blast.u.dens, dens == null ? 7.0 : dens);
        gl.uniform1f(P.blast.u.smallr, SMALLR);
        drawTo(U.b);
        swap();
        // A blast injects mass on purpose, so the running mass-conservation
        // readout has to re-baseline or it would report the injection as solver
        // error -- which is exactly the opposite of what it is there to measure.
        state.mass0 = null;
      },
      push(x, y, dx, dy, radius) {
        gl.useProgram(P.push.p);
        gl.uniform1i(P.push.u.uU, U.a.bind(0));
        gl.uniform2f(P.push.u.point, x, y);
        gl.uniform2f(P.push.u.delta, dx, dy);
        gl.uniform1f(P.push.u.aspect, grid.w / grid.h);
        gl.uniform1f(P.push.u.radius, radius == null ? 0.004 : radius);
        gl.uniform1f(P.push.u.smallr, SMALLR);
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

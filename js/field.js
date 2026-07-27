/* ============================================================================
   field.js — the live background.

   This is not a decorative particle toy. It is a 2-D incompressible
   Navier-Stokes solver coupled to a population of Lagrangian dust
   superparticles, run every frame on the GPU:

     gas    : semi-Lagrangian advection -> vorticity confinement ->
              divergence -> Jacobi pressure projection -> gradient subtract
              on a periodic, square-celled grid.
     dust   : x' = v, v' = (u_gas - v)/tau, integrated with the exact
              exponential (analytic) drag integrator so that stiff, tightly
              coupled grains stay stable at any timestep. Grains are
              rendered into a fading trail buffer, so what you see are real
              particle paths through a real velocity field.

   Because the drag law is honest, so is the behaviour: small Stokes number
   grains trace the gas, large Stokes number grains decouple, lag the flow,
   and concentrate preferentially in regions of high strain. That clumping
   is not painted on.

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

  const F_HEAD = `#version 300 es
precision highp float;
precision highp sampler2D;
in vec2 vUv;
out vec4 outColor;
`;

  // semi-Lagrangian backtrace. velocity is stored in cells/second; texel is
  // 1/grid, so vel*dt*texel is the displacement in uv.
  const F_ADVECT = F_HEAD + `
uniform sampler2D uVel, uSrc;
uniform vec2  texel;
uniform float dt, dissipation;
void main(){
  vec2 vel   = texture(uVel, vUv).xy;
  vec2 coord = vUv - dt * vel * texel;
  outColor   = texture(uSrc, coord) / (1.0 + dissipation * dt);
}`;

  const F_DIVERGENCE = F_HEAD + `
uniform sampler2D uVel;
uniform vec2 texel;
void main(){
  float L = texture(uVel, vUv - vec2(texel.x, 0.0)).x;
  float R = texture(uVel, vUv + vec2(texel.x, 0.0)).x;
  float B = texture(uVel, vUv - vec2(0.0, texel.y)).y;
  float T = texture(uVel, vUv + vec2(0.0, texel.y)).y;
  outColor = vec4(0.5 * ((R - L) + (T - B)), 0.0, 0.0, 1.0);
}`;

  const F_JACOBI = F_HEAD + `
uniform sampler2D uP, uDiv;
uniform vec2 texel;
void main(){
  float L = texture(uP, vUv - vec2(texel.x, 0.0)).x;
  float R = texture(uP, vUv + vec2(texel.x, 0.0)).x;
  float B = texture(uP, vUv - vec2(0.0, texel.y)).x;
  float T = texture(uP, vUv + vec2(0.0, texel.y)).x;
  float d = texture(uDiv, vUv).x;
  outColor = vec4((L + R + B + T - d) * 0.25, 0.0, 0.0, 1.0);
}`;

  const F_GRADSUB = F_HEAD + `
uniform sampler2D uP, uVel;
uniform vec2 texel;
void main(){
  float L = texture(uP, vUv - vec2(texel.x, 0.0)).x;
  float R = texture(uP, vUv + vec2(texel.x, 0.0)).x;
  float B = texture(uP, vUv - vec2(0.0, texel.y)).x;
  float T = texture(uP, vUv + vec2(0.0, texel.y)).x;
  vec2  v = texture(uVel, vUv).xy - 0.5 * vec2(R - L, T - B);
  outColor = vec4(v, 0.0, 1.0);
}`;

  const F_CURL = F_HEAD + `
uniform sampler2D uVel;
uniform vec2 texel;
void main(){
  float L = texture(uVel, vUv - vec2(texel.x, 0.0)).y;
  float R = texture(uVel, vUv + vec2(texel.x, 0.0)).y;
  float B = texture(uVel, vUv - vec2(0.0, texel.y)).x;
  float T = texture(uVel, vUv + vec2(0.0, texel.y)).x;
  outColor = vec4(0.5 * ((R - L) - (T - B)), 0.0, 0.0, 1.0);
}`;

  // Fedkiw-style vorticity confinement: f = eps * omega * (N.y, -N.x),
  // N = grad|omega| / |grad|omega||. Restores small-scale swirl lost to the
  // diffusive first-order advection at this resolution.
  const F_VORTICITY = F_HEAD + `
uniform sampler2D uVel, uCurl;
uniform vec2  texel;
uniform float curlAmt, dt;
void main(){
  float L = abs(texture(uCurl, vUv - vec2(texel.x, 0.0)).x);
  float R = abs(texture(uCurl, vUv + vec2(texel.x, 0.0)).x);
  float B = abs(texture(uCurl, vUv - vec2(0.0, texel.y)).x);
  float T = abs(texture(uCurl, vUv + vec2(0.0, texel.y)).x);
  float C = texture(uCurl, vUv).x;
  vec2  g = 0.5 * vec2(R - L, T - B);
  vec2  N = g / (length(g) + 1e-5);
  vec2  f = curlAmt * C * vec2(N.y, -N.x);
  outColor = vec4(texture(uVel, vUv).xy + f * dt, 0.0, 1.0);
}`;

  const F_SPLAT = F_HEAD + `
uniform sampler2D uTarget;
uniform float aspect, radius;
uniform vec2  point;
uniform vec3  color;
void main(){
  vec2 p = vUv - point;
  p.x *= aspect;
  vec3 s = exp(-dot(p, p) / radius) * color;
  outColor = vec4(texture(uTarget, vUv).xyz + s, 1.0);
}`;

  // A hyperbolic-tangent shear layer, localised near mid-screen: the
  // canonical Kelvin-Helmholtz setup, with a sinusoidal transverse seed so
  // the instability actually rolls up instead of sitting there.
  const F_SHEAR = F_HEAD + `
uniform sampler2D uVel;
uniform float amp, width, seed, phase;
void main(){
  vec2  v = texture(uVel, vUv).xy;
  float y = vUv.y - 0.5;
  float env = exp(-pow(y / (2.6 * width), 2.0));
  v.x += amp * tanh(y / width) * env;
  v.y += seed * sin(6.2831853 * (vUv.x * 3.0 + phase)) * env;
  outColor = vec4(v, 0.0, 1.0);
}`;

  const F_BAND = F_HEAD + `
uniform sampler2D uTarget;
uniform float width, strength;
uniform vec3  color;
void main(){
  float y = vUv.y - 0.5;
  float env = exp(-pow(y / (1.5 * width), 2.0));
  outColor = vec4(texture(uTarget, vUv).xyz + color * strength * env, 1.0);
}`;

  // radial blast: a crude but honest impulsive outflow, used as a
  // shock-tube-flavoured punctuation when entering some sections.
  const F_BLAST = F_HEAD + `
uniform sampler2D uVel;
uniform vec2  point;
uniform float aspect, amp, radius;
void main(){
  vec2 p = vUv - point;
  p.x *= aspect;
  float r = length(p) + 1e-4;
  float w = exp(-r * r / radius);
  outColor = vec4(texture(uVel, vUv).xy + amp * w * (p / r), 0.0, 1.0);
}`;

  // Illustrative guide field: relax the velocity toward its component along
  // B-hat. Stands in for magnetic tension suppressing cross-field motion,
  // which is what produces field-aligned striations.
  const F_GUIDEFIELD = F_HEAD + `
uniform sampler2D uVel;
uniform vec2  bhat;
uniform float rate, dt;
void main(){
  vec2 v = texture(uVel, vUv).xy;
  vec2 par = dot(v, bhat) * bhat;
  outColor = vec4(mix(v, par, clamp(rate * dt, 0.0, 1.0)), 0.0, 1.0);
}`;

  const F_FADE = F_HEAD + `
uniform sampler2D uSrc;
uniform float amount;
void main(){ outColor = texture(uSrc, vUv) * amount; }`;

  // Per-cell metrics, then a 4x4 reduction, so the on-screen diagnostics are
  // measured from the actual field rather than invented.
  const F_METRICS = F_HEAD + `
uniform sampler2D uVel;
uniform vec2 texel;
void main(){
  vec2  v = texture(uVel, vUv).xy;
  float L = texture(uVel, vUv - vec2(texel.x, 0.0)).x;
  float R = texture(uVel, vUv + vec2(texel.x, 0.0)).x;
  float B = texture(uVel, vUv - vec2(0.0, texel.y)).y;
  float T = texture(uVel, vUv + vec2(0.0, texel.y)).y;
  float div = 0.5 * ((R - L) + (T - B));
  float sp  = length(v);
  outColor = vec4(sp * sp, sp, abs(div), 1.0);
}`;

  const F_REDUCE = F_HEAD + `
uniform sampler2D uSrc;
uniform vec2 srcTexel;
void main(){
  vec4 acc = vec4(0.0);
  float mx = 0.0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      vec2 o = (vec2(float(i), float(j)) - 1.5) * srcTexel;
      vec4 s = texture(uSrc, vUv + o);
      acc += s;
      mx = max(mx, s.y);
    }
  }
  acc /= 16.0;
  outColor = vec4(acc.x, mx, acc.z, 1.0);
}`;

  const V_PARTICLE = `#version 300 es
precision highp float;
uniform sampler2D uPart;
uniform int   uPW;
uniform float uPointSize, uSpeedNorm;
out float vBright;
void main(){
  int id = gl_VertexID;
  ivec2 tc = ivec2(id % uPW, id / uPW);
  vec4 P = texelFetch(uPart, tc, 0);
  gl_Position = vec4(P.xy * 2.0 - 1.0, 0.0, 1.0);
  float sp = length(P.zw) * uSpeedNorm;
  vBright  = clamp(0.075 + sp * 0.62, 0.0, 1.05);
  gl_PointSize = uPointSize;
}`;

  const F_PARTICLE = `#version 300 es
precision highp float;
in float vBright;
out vec4 outColor;
uniform vec3 uColor;
void main(){
  vec2  d  = gl_PointCoord - 0.5;
  float a  = exp(-dot(d, d) * 11.0);
  outColor = vec4(uColor * vBright * a, a);
}`;

  const F_PARTICLE_UPDATE = F_HEAD + `
uniform sampler2D uPart, uVel;
uniform vec2  texel, stream;
uniform float dt, tau, seed, reseed, brown;

float hash(vec2 p){
  p = fract(p * vec2(123.34, 456.21));
  p += dot(p, p + 45.32);
  return fract(p.x * p.y);
}

void main(){
  vec4 P   = texture(uPart, vUv);
  vec2 pos = P.xy;
  vec2 vel = P.zw;

  vec2  ug = texture(uVel, pos).xy;
  float a  = exp(-dt / max(tau, 1.0e-4));   // exact drag integrator
  vel = ug + (vel - ug) * a;

  if (brown > 0.0) {
    vel += brown * (vec2(hash(vUv + seed), hash(vUv * 1.7 + seed + 3.1)) - 0.5);
  }

  pos = fract(pos + dt * (vel + stream) * texel);

  if (hash(vUv * 91.7 + floor(seed * 3.0)) < reseed) {
    pos = vec2(hash(vUv * 5.3 + seed), hash(vUv * 8.9 + seed * 2.0));
    vel = vec2(0.0);
  }
  outColor = vec4(pos, vel);
}`;

  const F_COMPOSITE = F_HEAD + `
uniform sampler2D uDye, uTrail;
uniform vec2  uRes;
uniform vec3  uBg, uTint;
uniform float uTime, uGrain, uDustGain, uDyeGain, uVignette;

float hash(vec2 p){
  p = fract(p * vec2(123.34, 456.21));
  p += dot(p, p + 45.32);
  return fract(p.x * p.y);
}

void main(){
  vec3 dye  = texture(uDye,   vUv).rgb * uDyeGain;
  vec3 dust = texture(uTrail, vUv).rgb * uDustGain;

  // faint cold gradient so the page never reads as flat black
  float g = smoothstep(1.15, -0.15, vUv.y + vUv.x * 0.22);
  vec3  c = uBg + uTint * g * 0.5;

  c += dye;
  c += dust;

  // Roll off on luminance, not per channel: a per-channel Reinhard drives
  // every bright pixel toward white and throws the section accent away.
  float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
  c *= 1.0 / (1.0 + 0.62 * l);

  float d = length((vUv - 0.5) * vec2(uRes.x / uRes.y, 1.0));
  c *= 1.0 - uVignette * smoothstep(0.30, 1.05, d);

  c += (hash(vUv * uRes + fract(uTime) * 71.3) - 0.5) * uGrain;
  outColor = vec4(max(c, 0.0), 1.0);
}`;

  // ------------------------------------------------------------- gl helpers

  function compile(gl, type, src) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
      throw new Error('shader: ' + gl.getShaderInfoLog(s) + '\n' + src.slice(0, 400));
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
      const name = gl.getActiveUniform(p, i).name.replace(/\[0\]$/, '');
      u[name] = gl.getUniformLocation(p, name);
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
    gl.viewport(0, 0, w, h);
    gl.clearColor(0, 0, 0, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);

    const ok = gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE;
    if (!ok) throw new Error('incomplete framebuffer ' + internal);

    return {
      tex, fbo, w, h, texel: [1 / w, 1 / h],
      bind(unit) { gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, tex); return unit; }
    };
  }

  function makeDouble(gl, w, h, internal, format, type, filter, wrap) {
    let a = makeFBO(gl, w, h, internal, format, type, filter, wrap);
    let b = makeFBO(gl, w, h, internal, format, type, filter, wrap);
    return {
      w, h, texel: a.texel,
      get read() { return a; },
      get write() { return b; },
      swap() { const t = a; a = b; b = t; }
    };
  }

  // -------------------------------------------------------------- presets

  // tau is the grain stopping time in seconds. The gas turns over on a
  // timescale of order 1 s here, so tau ~ 0.03 is a tightly coupled tracer
  // and tau ~ 1 is a grain that goes ballistic between eddies.
  const BASE = {
    accent: [0.85, 0.88, 0.95],
    tint: [0.020, 0.024, 0.040],
    tau: 0.10,
    dustGain: 1.00,
    dyeGain: 0.80,
    forceGain: 1.00,
    curl: 12.0,
    velDiss: 0.22,
    dyeDiss: 0.55,
    pointSize: 1.35,
    trailFade: 0.90,
    brown: 0.0,
    guide: 0.0,
    bhat: [1.0, 0.0],
    stream: 0.0,
    stirGain: 1.00,
    grain: 0.030,
    vignette: 0.62
  };

  const PRESETS = {
    hero:      { accent: [0.78, 0.84, 1.00], tau: 0.14, dustGain: 1.25, curl: 16, stirGain: 1.25, trailFade: 0.915 },
    bio:       { accent: [0.93, 0.89, 0.80], tau: 0.10, dustGain: 1.05, curl: 12 },
    picdust:   { accent: [1.00, 0.66, 0.26], tau: 0.030, dustGain: 1.45, curl: 15, pointSize: 1.2, trailFade: 0.90, stirGain: 1.15 },
    dfmm:      { accent: [0.42, 0.95, 0.78], tau: 0.85,  dustGain: 1.20, curl: 9,  pointSize: 1.5, trailFade: 0.945, brown: 0.05 },
    mhd:       { accent: [0.68, 0.60, 1.00], tau: 0.12,  dustGain: 1.15, curl: 8,  guide: 3.4, bhat: [0.94, 0.34], trailFade: 0.935 },
    cosmicray: { accent: [0.40, 0.90, 1.00], tau: 0.020, dustGain: 1.05, curl: 10, stream: 190, bhat: [0.94, 0.34], pointSize: 1.0, trailFade: 0.955, guide: 2.2 },
    phrike:    { accent: [1.00, 0.36, 0.55], tau: 0.06,  dustGain: 1.30, curl: 22, velDiss: 0.10, trailFade: 0.905, stirGain: 1.35 },
    papers:    { accent: [0.90, 0.87, 0.79], tau: 0.09,  dustGain: 0.95, curl: 11, stirGain: 0.85 },
    blog:      { accent: [0.80, 0.85, 0.94], tau: 0.11,  dustGain: 0.95, curl: 10, stirGain: 0.8 },
    iron:      { accent: [1.00, 0.20, 0.16], tau: 0.55,  dustGain: 1.30, curl: 5, velDiss: 0.60, pointSize: 1.9, trailFade: 0.90, stirGain: 0.65, vignette: 0.74 },
    cactus:    { accent: [0.62, 0.82, 0.44], tau: 0.95,  dustGain: 0.90, curl: 3, velDiss: 0.42, pointSize: 1.5, trailFade: 0.955, brown: 0.09, stirGain: 0.42, grain: 0.040 },
    skate:     { accent: [1.00, 0.56, 0.14], tau: 0.20,  dustGain: 1.25, curl: 17, velDiss: 0.34, pointSize: 1.3, trailFade: 0.895, stirGain: 1.1 },
    contact:   { accent: [0.90, 0.90, 0.86], tau: 0.10,  dustGain: 1.05, curl: 12, stirGain: 0.9 }
  };

  function resolvePreset(key) {
    return Object.assign({}, BASE, PRESETS[key] || PRESETS.hero);
  }

  // ------------------------------------------------------------- the field

  function Field(canvas) {
    const gl = canvas.getContext('webgl2', {
      alpha: false, depth: false, stencil: false, antialias: false,
      premultipliedAlpha: false, preserveDrawingBuffer: false,
      powerPreference: 'high-performance'
    });
    if (!gl) throw new Error('no webgl2');
    if (!gl.getExtension('EXT_color_buffer_float')) throw new Error('no float render targets');
    gl.getExtension('OES_texture_float_linear');

    const reduced = global.matchMedia &&
      global.matchMedia('(prefers-reduced-motion: reduce)').matches;
    const mobile = Math.min(global.innerWidth, global.innerHeight) < 720;

    // quality tiers: [gridH, dyeScale, particleSide, trailScale, dprCap]
    const TIERS = [
      [208, 2.0, 232, 0.80, 1.75],
      [176, 1.8, 192, 0.70, 1.50],
      [144, 1.6, 152, 0.58, 1.25],
      [112, 1.4, 112, 0.50, 1.00]
    ];
    let tier = mobile ? 2 : 0;

    const quad = gl.createVertexArray();
    gl.bindVertexArray(quad);
    const vb = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vb);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    gl.bindVertexArray(null);

    const emptyVAO = gl.createVertexArray();

    const P = {
      advect:   program(gl, VERT, F_ADVECT),
      diverge:  program(gl, VERT, F_DIVERGENCE),
      jacobi:   program(gl, VERT, F_JACOBI),
      gradsub:  program(gl, VERT, F_GRADSUB),
      curl:     program(gl, VERT, F_CURL),
      vort:     program(gl, VERT, F_VORTICITY),
      splat:    program(gl, VERT, F_SPLAT),
      shear:    program(gl, VERT, F_SHEAR),
      band:     program(gl, VERT, F_BAND),
      blast:    program(gl, VERT, F_BLAST),
      guide:    program(gl, VERT, F_GUIDEFIELD),
      fade:     program(gl, VERT, F_FADE),
      metrics:  program(gl, VERT, F_METRICS),
      reduce:   program(gl, VERT, F_REDUCE),
      comp:     program(gl, VERT, F_COMPOSITE),
      pupdate:  program(gl, VERT, F_PARTICLE_UPDATE),
      pdraw:    program(gl, V_PARTICLE, F_PARTICLE)
    };

    let vel, dye, prs, div, crl, part, trail, met, red1, red2, red3;
    let grid = { w: 0, h: 0 }, dyeRes = { w: 0, h: 0 }, pSide = 0, nPart = 0;
    let dpr = 1;

    const state = {
      preset: resolvePreset('hero'),
      target: resolvePreset('hero'),
      blend: 1,
      time: 0,
      running: true,
      fps: 60,
      diag: { rms: 0, max: 0, div: 0 }
    };

    const pending = { splats: [], shear: 0, blasts: [] };
    const RGBA16 = [gl.RGBA16F, gl.RGBA, gl.HALF_FLOAT];
    const RG16 = [gl.RG16F, gl.RG, gl.HALF_FLOAT];
    const R16 = [gl.R16F, gl.RED, gl.HALF_FLOAT];

    function drawQuad(fbo) {
      if (fbo) { gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.fbo); gl.viewport(0, 0, fbo.w, fbo.h); }
      else { gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.viewport(0, 0, canvas.width, canvas.height); }
      gl.bindVertexArray(quad);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function allocate() {
      const t = TIERS[tier];
      dpr = Math.min(global.devicePixelRatio || 1, t[4]);
      const cw = Math.max(1, Math.round(canvas.clientWidth * dpr));
      const ch = Math.max(1, Math.round(canvas.clientHeight * dpr));
      canvas.width = cw; canvas.height = ch;

      const aspect = cw / ch;
      grid.h = t[0];
      grid.w = Math.max(48, Math.round(t[0] * aspect));
      dyeRes.h = Math.round(grid.h * t[1]);
      dyeRes.w = Math.round(grid.w * t[1]);
      pSide = t[2];
      nPart = pSide * pSide;

      const R = gl.REPEAT, C = gl.CLAMP_TO_EDGE, L = gl.LINEAR, N = gl.NEAREST;
      vel = makeDouble(gl, grid.w, grid.h, RG16[0], RG16[1], RG16[2], L, R);
      dye = makeDouble(gl, dyeRes.w, dyeRes.h, RGBA16[0], RGBA16[1], RGBA16[2], L, R);
      prs = makeDouble(gl, grid.w, grid.h, R16[0], R16[1], R16[2], L, R);
      div = makeFBO(gl, grid.w, grid.h, R16[0], R16[1], R16[2], L, R);
      crl = makeFBO(gl, grid.w, grid.h, R16[0], R16[1], R16[2], L, R);
      // 32-bit: |u|^2 overflows half-float once speeds pass ~256 cells/s
      met = makeFBO(gl, grid.w, grid.h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);

      const r1w = Math.max(1, grid.w >> 2), r1h = Math.max(1, grid.h >> 2);
      const r2w = Math.max(1, r1w >> 2), r2h = Math.max(1, r1h >> 2);
      const r3w = Math.max(1, r2w >> 2), r3h = Math.max(1, r2h >> 2);
      red1 = makeFBO(gl, r1w, r1h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);
      red2 = makeFBO(gl, r2w, r2h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);
      red3 = makeFBO(gl, r3w, r3h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);

      const tw = Math.max(1, Math.round(cw * t[3])), th = Math.max(1, Math.round(ch * t[3]));
      trail = makeDouble(gl, tw, th, RGBA16[0], RGBA16[1], RGBA16[2], L, C);

      part = makeDouble(gl, pSide, pSide, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);
      seedParticles();
      seedField();
    }

    function seedParticles() {
      const data = new Float32Array(nPart * 4);
      for (let i = 0; i < nPart; i++) {
        data[i * 4 + 0] = Math.random();
        data[i * 4 + 1] = Math.random();
        data[i * 4 + 2] = 0;
        data[i * 4 + 3] = 0;
      }
      for (const f of [part.read, part.write]) {
        gl.bindTexture(gl.TEXTURE_2D, f.tex);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, pSide, pSide, 0, gl.RGBA, gl.FLOAT, data);
      }
    }

    // Spin the field up before the first frame is ever shown, so the page
    // opens on developed turbulence rather than on a blank box relaxing.
    function seedField() {
      const a = state.target.accent;
      for (let i = 0; i < 11; i++) {
        const ph = i * 2.399;
        pending.splats.push({
          x: 0.5 + 0.38 * Math.cos(ph), y: 0.5 + 0.34 * Math.sin(ph * 1.31),
          dx: 90 * Math.sin(ph * 1.7), dy: 90 * Math.cos(ph * 2.3),
          color: [a[0] * 0.5, a[1] * 0.5, a[2] * 0.5], radius: 0.013
        });
      }
      const dt = 1 / 60;
      for (let i = 0; i < 90; i++) {
        state.time += dt;
        step(dt);
      }
      state.time = 0;
    }

    // ------------------------------------------------------------ stepping

    function lerpPreset(dt) {
      if (state.blend >= 1) return;
      state.blend = Math.min(1, state.blend + dt * 1.35);
      const k = state.blend;
      const e = k * k * (3 - 2 * k);
      const cur = state.preset, tgt = state.target;
      const out = {};
      for (const key in tgt) {
        const a = cur[key], b = tgt[key];
        if (Array.isArray(b)) out[key] = b.map((v, i) => a[i] + (v - a[i]) * e);
        else if (typeof b === 'number') out[key] = a + (b - a) * e;
        else out[key] = b;
      }
      state.live = out;
    }

    function activePreset() { return state.live || state.preset; }

    function applySplats(pr) {
      if (!pending.splats.length) return;
      const aspect = grid.w / grid.h;
      for (const s of pending.splats) {
        gl.useProgram(P.splat.p);
        gl.uniform1i(P.splat.u.uTarget, vel.read.bind(0));
        gl.uniform1f(P.splat.u.aspect, aspect);
        gl.uniform2f(P.splat.u.point, s.x, s.y);
        gl.uniform3f(P.splat.u.color, s.dx, s.dy, 0);
        gl.uniform1f(P.splat.u.radius, s.radius);
        drawQuad(vel.write); vel.swap();

        gl.uniform1i(P.splat.u.uTarget, dye.read.bind(0));
        gl.uniform3f(P.splat.u.color, s.color[0], s.color[1], s.color[2]);
        gl.uniform1f(P.splat.u.radius, s.radius * 1.5);
        drawQuad(dye.write); dye.swap();
      }
      pending.splats.length = 0;
    }

    function applyShear(pr) {
      if (Math.abs(pending.shear) < 1e-4) return;
      const amp = Math.max(-170, Math.min(170, pending.shear));
      gl.useProgram(P.shear.p);
      gl.uniform1i(P.shear.u.uVel, vel.read.bind(0));
      gl.uniform1f(P.shear.u.amp, amp);
      gl.uniform1f(P.shear.u.width, 0.055);
      gl.uniform1f(P.shear.u.seed, amp * 0.30);
      gl.uniform1f(P.shear.u.phase, state.time * 0.21);
      drawQuad(vel.write); vel.swap();

      const a = pr.accent, s = Math.min(0.95, Math.abs(amp) / 170 * 0.95);
      gl.useProgram(P.band.p);
      gl.uniform1i(P.band.u.uTarget, dye.read.bind(0));
      gl.uniform1f(P.band.u.width, 0.055);
      gl.uniform1f(P.band.u.strength, s);
      gl.uniform3f(P.band.u.color, a[0], a[1], a[2]);
      drawQuad(dye.write); dye.swap();

      pending.shear *= 0.55;
      if (Math.abs(pending.shear) < 1e-3) pending.shear = 0;
    }

    function applyBlasts() {
      if (!pending.blasts.length) return;
      const aspect = grid.w / grid.h;
      for (const b of pending.blasts) {
        gl.useProgram(P.blast.p);
        gl.uniform1i(P.blast.u.uVel, vel.read.bind(0));
        gl.uniform2f(P.blast.u.point, b.x, b.y);
        gl.uniform1f(P.blast.u.aspect, aspect);
        gl.uniform1f(P.blast.u.amp, b.amp);
        gl.uniform1f(P.blast.u.radius, b.radius);
        drawQuad(vel.write); vel.swap();
      }
      pending.blasts.length = 0;
    }

    // slow, solenoidal-ish large-scale driving: two counter-rotating
    // wandering stirrers. Without it the field decays and the page dies.
    function stir(dt, pr) {
      const t = state.time;
      // tuned so that |u|_max settles near 250-350 cells/s, i.e. a Courant
      // displacement of a few cells per step rather than tens
      const g = pr.stirGain * 17;
      for (let i = 0; i < 2; i++) {
        const ph = t * (0.11 + i * 0.052) + i * 2.4;
        const x = 0.5 + 0.31 * Math.cos(ph * 1.7 + i);
        const y = 0.5 + 0.27 * Math.sin(ph * 1.31 + i * 2.0);
        const dir = i === 0 ? 1 : -1;
        const a = pr.accent;
        pending.splats.push({
          x, y,
          dx: dir * g * Math.cos(ph * 2.3),
          dy: dir * g * Math.sin(ph * 2.9),
          color: [a[0] * 0.048, a[1] * 0.048, a[2] * 0.048],
          radius: 0.024
        });
      }
    }

    function step(dt) {
      const pr = activePreset();

      stir(dt, pr);
      applySplats(pr);
      applyShear(pr);
      applyBlasts();

      // vorticity confinement
      if (pr.curl > 0.01) {
        gl.useProgram(P.curl.p);
        gl.uniform1i(P.curl.u.uVel, vel.read.bind(0));
        gl.uniform2f(P.curl.u.texel, vel.texel[0], vel.texel[1]);
        drawQuad(crl);

        gl.useProgram(P.vort.p);
        gl.uniform1i(P.vort.u.uVel, vel.read.bind(0));
        gl.uniform1i(P.vort.u.uCurl, crl.bind(1));
        gl.uniform2f(P.vort.u.texel, vel.texel[0], vel.texel[1]);
        gl.uniform1f(P.vort.u.curlAmt, pr.curl);
        gl.uniform1f(P.vort.u.dt, dt);
        drawQuad(vel.write); vel.swap();
      }

      if (pr.guide > 0.01) {
        const b = pr.bhat, n = Math.hypot(b[0], b[1]) || 1;
        gl.useProgram(P.guide.p);
        gl.uniform1i(P.guide.u.uVel, vel.read.bind(0));
        gl.uniform2f(P.guide.u.bhat, b[0] / n, b[1] / n);
        gl.uniform1f(P.guide.u.rate, pr.guide);
        gl.uniform1f(P.guide.u.dt, dt);
        drawQuad(vel.write); vel.swap();
      }

      // projection
      gl.useProgram(P.diverge.p);
      gl.uniform1i(P.diverge.u.uVel, vel.read.bind(0));
      gl.uniform2f(P.diverge.u.texel, vel.texel[0], vel.texel[1]);
      drawQuad(div);

      gl.useProgram(P.jacobi.p);
      gl.uniform2f(P.jacobi.u.texel, vel.texel[0], vel.texel[1]);
      gl.uniform1i(P.jacobi.u.uDiv, div.bind(1));
      const iters = tier <= 1 ? 30 : 20;
      for (let i = 0; i < iters; i++) {
        gl.uniform1i(P.jacobi.u.uP, prs.read.bind(0));
        drawQuad(prs.write); prs.swap();
      }

      gl.useProgram(P.gradsub.p);
      gl.uniform2f(P.gradsub.u.texel, vel.texel[0], vel.texel[1]);
      gl.uniform1i(P.gradsub.u.uP, prs.read.bind(0));
      gl.uniform1i(P.gradsub.u.uVel, vel.read.bind(1));
      drawQuad(vel.write); vel.swap();

      // advect velocity, then dye
      gl.useProgram(P.advect.p);
      gl.uniform2f(P.advect.u.texel, vel.texel[0], vel.texel[1]);
      gl.uniform1f(P.advect.u.dt, dt);
      gl.uniform1f(P.advect.u.dissipation, pr.velDiss);
      gl.uniform1i(P.advect.u.uVel, vel.read.bind(0));
      gl.uniform1i(P.advect.u.uSrc, vel.read.bind(0));
      drawQuad(vel.write); vel.swap();

      gl.uniform1f(P.advect.u.dissipation, pr.dyeDiss);
      gl.uniform1i(P.advect.u.uVel, vel.read.bind(0));
      gl.uniform1i(P.advect.u.uSrc, dye.read.bind(1));
      drawQuad(dye.write); dye.swap();

      // dust
      gl.useProgram(P.pupdate.p);
      gl.uniform1i(P.pupdate.u.uPart, part.read.bind(0));
      gl.uniform1i(P.pupdate.u.uVel, vel.read.bind(1));
      gl.uniform2f(P.pupdate.u.texel, vel.texel[0], vel.texel[1]);
      gl.uniform1f(P.pupdate.u.dt, dt);
      gl.uniform1f(P.pupdate.u.tau, pr.tau);
      gl.uniform1f(P.pupdate.u.seed, state.time);
      gl.uniform1f(P.pupdate.u.reseed, 0.0018);
      gl.uniform1f(P.pupdate.u.brown, pr.brown * 240);
      const b = pr.bhat, bn = Math.hypot(b[0], b[1]) || 1;
      gl.uniform2f(P.pupdate.u.stream, pr.stream * b[0] / bn, pr.stream * b[1] / bn);
      drawQuad(part.write); part.swap();
    }

    function renderDust(pr) {
      // fade the trail buffer, then splat this frame's grains into it
      gl.useProgram(P.fade.p);
      gl.uniform1i(P.fade.u.uSrc, trail.read.bind(0));
      gl.uniform1f(P.fade.u.amount, pr.trailFade);
      drawQuad(trail.write);

      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE);
      gl.useProgram(P.pdraw.p);
      gl.uniform1i(P.pdraw.u.uPart, part.read.bind(0));
      gl.uniform1i(P.pdraw.u.uPW, pSide);
      gl.uniform1f(P.pdraw.u.uPointSize, pr.pointSize * dpr * TIERS[tier][3] * 1.4);
      gl.uniform1f(P.pdraw.u.uSpeedNorm, 1 / 130);
      const a = pr.accent;
      gl.uniform3f(P.pdraw.u.uColor, a[0], a[1], a[2]);
      gl.bindVertexArray(emptyVAO);
      gl.drawArrays(gl.POINTS, 0, nPart);
      gl.disable(gl.BLEND);
      trail.swap();
    }

    function composite(pr) {
      gl.useProgram(P.comp.p);
      gl.uniform1i(P.comp.u.uDye, dye.read.bind(0));
      gl.uniform1i(P.comp.u.uTrail, trail.read.bind(1));
      gl.uniform2f(P.comp.u.uRes, canvas.width, canvas.height);
      gl.uniform3f(P.comp.u.uBg, 0.014, 0.015, 0.021);
      gl.uniform3f(P.comp.u.uTint, pr.tint[0], pr.tint[1], pr.tint[2]);
      gl.uniform1f(P.comp.u.uTime, state.time);
      gl.uniform1f(P.comp.u.uGrain, pr.grain);
      gl.uniform1f(P.comp.u.uDustGain, pr.dustGain);
      gl.uniform1f(P.comp.u.uDyeGain, pr.dyeGain);
      gl.uniform1f(P.comp.u.uVignette, pr.vignette);
      drawQuad(null);
    }

    // measured, not invented
    const readBuf = new Float32Array(4 * 64);
    let measureCountdown = 2;
    function measure() {
      gl.useProgram(P.metrics.p);
      gl.uniform1i(P.metrics.u.uVel, vel.read.bind(0));
      gl.uniform2f(P.metrics.u.texel, vel.texel[0], vel.texel[1]);
      drawQuad(met);

      let src = met, dst = red1;
      for (const d of [red1, red2, red3]) {
        gl.useProgram(P.reduce.p);
        gl.uniform1i(P.reduce.u.uSrc, src.bind(0));
        gl.uniform2f(P.reduce.u.srcTexel, 1 / src.w, 1 / src.h);
        drawQuad(d);
        src = d;
      }
      gl.bindFramebuffer(gl.FRAMEBUFFER, src.fbo);
      const n = Math.min(64, src.w * src.h);
      try {
        gl.readPixels(0, 0, src.w, src.h, gl.RGBA, gl.FLOAT, readBuf);
      } catch (e) { return; }
      let s2 = 0, mx = 0, dv = 0;
      for (let i = 0; i < n; i++) {
        s2 += readBuf[i * 4 + 0];
        mx = Math.max(mx, readBuf[i * 4 + 1]);
        dv += readBuf[i * 4 + 2];
      }
      const rms = Math.sqrt(Math.max(0, s2 / n));
      state.diag.rms = isFinite(rms) ? rms : state.diag.rms;
      state.diag.max = isFinite(mx) ? mx : state.diag.max;
      state.diag.div = isFinite(dv) ? dv / n : state.diag.div;
    }

    // ---------------------------------------------------------------- loop

    let last = 0, acc = 0, frames = 0, fpsT = 0, slow = 0;
    let frozen = false;

    function frame(now) {
      if (!state.running) { last = now; requestAnimationFrame(frame); return; }
      if (!last) last = now;
      let dt = (now - last) / 1000;
      last = now;
      dt = Math.min(dt, 1 / 45);
      state.time += dt;

      frames++; fpsT += dt;
      if (fpsT > 0.5) {
        state.fps = frames / fpsT;
        frames = 0; fpsT = 0;
        if (state.fps < 34 && tier < TIERS.length - 1) {
          if (++slow >= 3) { slow = 0; tier++; allocate(); }
        } else if (state.fps > 52) slow = 0;
      }

      lerpPreset(dt);
      const pr = activePreset();

      if (!frozen) {
        step(dt);
        renderDust(pr);
        // readPixels forces a pipeline sync, so keep it rare
        if (--measureCountdown <= 0) { measureCountdown = 42; measure(); }
      }
      composite(pr);

      if (reduced && state.time > 3.0) frozen = true;
      requestAnimationFrame(frame);
    }

    // ------------------------------------------------------------- public

    let resizeTimer = null;
    function onResize() {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(() => { allocate(); frozen = false; }, 180);
    }

    allocate();
    requestAnimationFrame(frame);
    global.addEventListener('resize', onResize);
    document.addEventListener('visibilitychange', () => { state.running = !document.hidden; last = 0; });

    return {
      kind: 'webgl2',
      splat(x, y, dx, dy, radius) {
        const pr = activePreset(), a = pr.accent;
        const sp = Math.hypot(dx, dy);
        if (sp > 320) { dx *= 320 / sp; dy *= 320 / sp; }
        const m = Math.min(1, sp / 320);
        pending.splats.push({
          x, y, dx, dy,
          color: [a[0] * (0.30 + m * 1.1), a[1] * (0.30 + m * 1.1), a[2] * (0.30 + m * 1.1)],
          radius: radius || 0.0085
        });
        if (pending.splats.length > 90) pending.splats.splice(0, pending.splats.length - 90);
        frozen = false;
      },
      shear(amount) { pending.shear += amount; frozen = false; },
      blast(x, y, amp, radius) {
        pending.blasts.push({ x, y, amp: amp || 700, radius: radius || 0.010 });
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
        return {
          gridW: grid.w, gridH: grid.h, nPart, fps: state.fps, tier,
          tau: activePreset().tau, dpr,
          rms: state.diag.rms, max: state.diag.max, div: state.diag.div,
          jacobi: tier <= 1 ? 30 : 20
        };
      }
    };
  }

  // ----------------------------------------------------- canvas2d fallback

  function Fallback(canvas) {
    const ctx = canvas.getContext('2d');
    const reduced = global.matchMedia &&
      global.matchMedia('(prefers-reduced-motion: reduce)').matches;
    let w = 0, h = 0, dpr = Math.min(global.devicePixelRatio || 1, 1.5);
    const N = Math.min(2600, Math.round((global.innerWidth * global.innerHeight) / 900));
    const px = new Float32Array(N), py = new Float32Array(N);
    const vx = new Float32Array(N), vy = new Float32Array(N);
    let accent = [0.78, 0.84, 1.0], tau = 0.12, t = 0, shear = 0;

    function size() {
      w = canvas.width = Math.round(canvas.clientWidth * dpr);
      h = canvas.height = Math.round(canvas.clientHeight * dpr);
      ctx.fillStyle = '#04050a';
      ctx.fillRect(0, 0, w, h);
    }
    for (let i = 0; i < N; i++) { px[i] = Math.random(); py[i] = Math.random(); }
    size();
    global.addEventListener('resize', size);

    // divergence-free by construction: u = curl of a scalar stream function
    function flow(x, y, tt) {
      const s = 2.7;
      const psx = Math.sin(s * y + tt * 0.31) * Math.cos(s * 0.7 * x - tt * 0.23);
      const psy = Math.cos(s * x - tt * 0.19) * Math.sin(s * 0.6 * y + tt * 0.27);
      return [0.11 * (psx + 0.5 * psy), 0.11 * (psy - 0.5 * psx) + shear * (y - 0.5) * 0.9];
    }

    function frame() {
      t += 1 / 60;
      shear *= 0.94;
      ctx.fillStyle = 'rgba(4,5,10,0.11)';
      ctx.fillRect(0, 0, w, h);
      const a = 1 - Math.exp(-(1 / 60) / tau);
      const col = `rgba(${(accent[0] * 255) | 0},${(accent[1] * 255) | 0},${(accent[2] * 255) | 0},0.55)`;
      ctx.fillStyle = col;
      for (let i = 0; i < N; i++) {
        const g = flow(px[i], py[i], t);
        vx[i] += (g[0] - vx[i]) * a;
        vy[i] += (g[1] - vy[i]) * a;
        px[i] += vx[i] / 60; py[i] += vy[i] / 60;
        if (px[i] < 0) px[i] += 1; if (px[i] > 1) px[i] -= 1;
        if (py[i] < 0) py[i] += 1; if (py[i] > 1) py[i] -= 1;
        ctx.fillRect(px[i] * w, py[i] * h, dpr, dpr);
      }
      if (!reduced) requestAnimationFrame(frame);
    }
    requestAnimationFrame(frame);

    return {
      kind: 'canvas2d',
      splat(x, y, dx, dy) { shear += dy / 40000; },
      shear(amount) { shear += amount / 900; },
      blast() {},
      setPreset(key) { const p = resolvePreset(key); accent = p.accent; tau = p.tau; },
      accentOf(key) { return resolvePreset(key).accent; },
      stats() { return { gridW: 0, gridH: 0, nPart: N, fps: 60, tier: 9, tau, dpr, rms: 0, max: 0, div: 0, jacobi: 0 }; }
    };
  }

  global.Field = {
    create(canvas) {
      try { return Field(canvas); }
      catch (e) {
        if (global.console) console.warn('[field] falling back to canvas2d:', e.message);
        try { return Fallback(canvas); } catch (e2) { return null; }
      }
    }
  };
})(window);

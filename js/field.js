/* ============================================================================
   field.js — the live background.

   This is not a decorative particle toy. It is a 2-D incompressible
   Navier-Stokes solver coupled to a population of Lagrangian dust
   superparticles, run every frame on the GPU:

     gas    : semi-Lagrangian advection -> vorticity confinement ->
              divergence -> Jacobi pressure projection -> gradient subtract
              (with a speed governor) on a periodic, square-celled grid.
     dust   : x' = v, v' = (u_gas - v)/tau, advanced by a first-order implicit
              (backward Euler) update. The population is polydisperse: every
              grain draws a stopping time from a log-uniform spectrum 1.4
              decades wide, reproduced from its own texel by a hash rather
              than stored, so the update stays one divide, one multiply-add
              and one multiply per grain -- the cheapest form that is still
              L-stable, which matters because the tightly coupled end of the
              spectrum has tau far below the step.

   Grains are drawn as one soft sprite each, straight to the screen, and what
   sets a grain's brightness and size is its *drift*, |v - u_gas|. That is the
   only quantity in the drag law that does work on a grain: it vanishes for
   perfectly coupled dust and is largest exactly where decoupled grains pile
   up, so the picture lights up its own physics.

   Because the drag law is honest, so is the behaviour: small Stokes number
   grains trace the gas, large Stokes number grains decouple, lag the flow,
   and concentrate preferentially in regions of high strain. That clumping
   is not painted on.

   The simulation clock deliberately runs at one tenth of wall time, and the
   frame rate is capped at 30, so the field drifts rather than churns. Every
   rate-like term is scaled by the step, so neither choice changes the physics.

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

  // Damped Jacobi. omega = 4/5 is the optimal smoothing factor for the 2-D
  // five-point Laplacian, which is what a multigrid smoother wants; omega = 1
  // recovers plain Jacobi. Used at every level of the V-cycle.
  const F_JACOBI = F_HEAD + `
uniform sampler2D uP, uDiv;
uniform vec2  texel;
uniform float omega;
void main(){
  float L = texture(uP, vUv - vec2(texel.x, 0.0)).x;
  float R = texture(uP, vUv + vec2(texel.x, 0.0)).x;
  float B = texture(uP, vUv - vec2(0.0, texel.y)).x;
  float T = texture(uP, vUv + vec2(0.0, texel.y)).x;
  float d = texture(uDiv, vUv).x;
  float p = texture(uP, vUv).x;
  outColor = vec4(p + omega * ((L + R + B + T - d) * 0.25 - p), 0.0, 0.0, 1.0);
}`;

  // Fused residual + restriction. One coarse cell gathers the residual
  // r = d - A p of the four fine cells it covers. The 1/4 from averaging them
  // and the factor 4 from the coarser cell spacing (h_c^2/h_f^2, because this
  // operator has h^2 absorbed into the right-hand side) cancel exactly, so the
  // coarse right-hand side is simply the sum. Fusing the two passes also means
  // no residual buffer has to exist at any level.
  const F_MG_RESTRICT = F_HEAD + `
uniform sampler2D uP, uDiv;
uniform ivec2 fineSize;
float pf(ivec2 c){
  return texelFetch(uP, ivec2((c.x + fineSize.x) % fineSize.x,
                              (c.y + fineSize.y) % fineSize.y), 0).x;
}
void main(){
  ivec2 f0 = ivec2(gl_FragCoord.xy) * 2;
  float sum = 0.0;
  for (int b = 0; b < 2; b++) {
    for (int a = 0; a < 2; a++) {
      ivec2 f = f0 + ivec2(a, b);
      float lap = pf(f + ivec2(-1, 0)) + pf(f + ivec2(1, 0))
                + pf(f + ivec2(0, -1)) + pf(f + ivec2(0, 1)) - 4.0 * pf(f);
      sum += texelFetch(uDiv, f, 0).x - lap;
    }
  }
  outColor = vec4(sum, 0.0, 0.0, 1.0);
}`;

  // Prolongation, added straight onto the fine pressure. A cell-centred fine
  // cell sits at coarse coordinate (i + 0.5)/2 - 0.5, which always lands a
  // quarter or three quarters of the way between two coarse centres -- the
  // classic 3/4-1/4 bilinear stencil. Written out with texelFetch rather than
  // leaning on hardware filtering, because linear filtering of 32-bit float
  // textures needs an extension that is not guaranteed.
  const F_MG_PROLONG = F_HEAD + `
uniform sampler2D uP, uCoarse;
uniform ivec2 coarseSize;
float cf(ivec2 c){
  return texelFetch(uCoarse, ivec2((c.x + coarseSize.x) % coarseSize.x,
                                   (c.y + coarseSize.y) % coarseSize.y), 0).x;
}
void main(){
  ivec2 f  = ivec2(gl_FragCoord.xy);
  vec2  cc = (vec2(f) + 0.5) * 0.5 - 0.5;
  ivec2 i0 = ivec2(floor(cc));
  vec2  t  = cc - vec2(i0);
  float v  = mix(mix(cf(i0),               cf(i0 + ivec2(1, 0)), t.x),
                 mix(cf(i0 + ivec2(0, 1)), cf(i0 + ivec2(1, 1)), t.x), t.y);
  outColor = vec4(texelFetch(uP, f, 0).x + v, 0.0, 0.0, 1.0);
}`;

  // Gradient subtract, plus a smooth speed governor. The clock runs an order of
  // magnitude slower than wall time, which means injected energy also lingers an
  // order of magnitude longer; without a ceiling, a minute of enthusiastic mouse
  // movement would pump the field straight back up to a speed we deliberately
  // left behind. The limiter asymptotes to vmax instead of clipping, so it never
  // shows an edge, and the next projection cleans up the divergence it makes.
  const F_GRADSUB = F_HEAD + `
uniform sampler2D uP, uVel;
uniform vec2  texel;
uniform float vmax;
void main(){
  float L = texture(uP, vUv - vec2(texel.x, 0.0)).x;
  float R = texture(uP, vUv + vec2(texel.x, 0.0)).x;
  float B = texture(uP, vUv - vec2(0.0, texel.y)).x;
  float T = texture(uP, vUv + vec2(0.0, texel.y)).x;
  vec2  v = texture(uVel, vUv).xy - 0.5 * vec2(R - L, T - B);
  float sp = length(v);
  v *= 1.0 / (1.0 + max(0.0, sp / vmax - 1.0));
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
  v.y += seed * sin(6.2831853 * (vUv.x * 5.0 + phase)) * env;
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

  // The grain size spectrum. Real dust populations are polydisperse, and one
  // stopping time for everything is the single least honest thing a picture
  // like this can do: it collapses the coupled and decoupled families into one
  // indistinguishable fluid. The rank is hashed out of the grain's own texel,
  // so the spectrum costs no storage, never has to be kept in step with the
  // particle state, and survives reseeding -- the distribution is a property
  // of the population, not of any individual grain's history.
  // Integer hash, not the usual fract(sin(...)) or fract(p*k) float trick. Those
  // lose their entropy as the input grows: fract() of a float near 5e5 is
  // quantised to steps of ~0.03, so the result stops being uniform on [0,1] and
  // a test like `hash(...) < 1e-4` fires orders of magnitude too often. That is
  // not a subtle bias -- it made the reseed rate 20x what it was set to. This is
  // exact in uint arithmetic at every magnitude, and being integer it is also
  // bit-identical between the vertex and fragment stages, which matters because
  // two different shaders have to agree on which grain is which.
  const S_SPECTRUM = `
uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float hash1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}
// rank in [0,1] -> tau multiplier, log-uniform over +-0.7 dex about the
// preset's tau, so the preset still names the median grain.
float stRank(ivec2 tc){ return hash1(uvec2(tc), 7u); }
float tauMul(float s){ return exp((s * 2.0 - 1.0) * 1.6118); }
`;

  const V_PARTICLE = `#version 300 es
precision highp float;
uniform sampler2D uPart, uVel;
uniform int   uPW;
uniform vec2  uStream;
uniform float uPointSize, uDriftNorm, uAlpha, uVignette, uAspect;
out float vBright;
out float vDrift;
` + S_SPECTRUM + `
void main(){
  int id = gl_VertexID;
  ivec2 tc = ivec2(id % uPW, id / uPW);
  vec4  P   = texelFetch(uPart, tc, 0);
  vec2  pos = P.xy;
  float st  = stRank(tc);

  // Drift, not speed. A grain swept along with the gas is doing nothing
  // interesting however fast it is moving; a grain slipping through the gas is
  // the whole subject. One texture fetch buys the difference. uStream is part
  // of the grain's real velocity, so it belongs inside the difference too.
  vec2  ug = texture(uVel, pos).xy;
  vDrift   = clamp(length(P.zw + uStream - ug) * uDriftNorm, 0.0, 1.0);

  gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);

  // Sprites go straight to the screen, so they carry the vignette themselves.
  vec2  q  = (pos - 0.5) * vec2(uAspect, 1.0);
  float vg = 1.0 - uVignette * smoothstep(0.30, 1.05, length(q));
  vBright  = uAlpha * vg * (0.30 + 0.70 * vDrift) * (0.72 + 0.58 * st);
  gl_PointSize = uPointSize * (0.80 + 0.90 * st + 0.60 * vDrift);
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
  // Compact core, and the (1 - r2) factor takes the sprite to exactly zero at
  // the rim: a bare Gaussian leaves a step there, and with a hundred thousand
  // overlapping grains that step is visible as square tiling.
  float a = exp(-r2 * 2.7) * (1.0 - r2);
  // Drifting grains do not just brighten, they go hot: drag work against the
  // gas is the only thing heating a grain in this picture.
  vec3  c = mix(uColor, mix(uColor, vec3(1.0), 0.42), vDrift * vDrift);
  float b = a * vBright;
  outColor = vec4(c * b, b);
}`;

  // Cheapest defensible drag: first-order implicit (backward Euler). For
  // dv/dt = (u - v)/tau,
  //
  //     v^{n+1} = (v^n + (dt/tau) u) / (1 + dt/tau)
  //
  // tau varies across the population but not in time, so the only frame scalar
  // needed is dragA = dt/tau_median; a grain rescales it by its own place in
  // the size spectrum. The per-grain cost is one divide, one multiply-add and
  // one multiply. No exp(). Backward Euler is only first-order accurate but it
  // is L-stable, so arbitrarily stiff grains (dt >> tau) relax monotonically to
  // v = u instead of ringing or blowing up -- which is what lets one step serve
  // the whole spectrum at once.
  const F_PARTICLE_UPDATE = F_HEAD + `
uniform sampler2D uPart, uVel;
uniform vec2  stream, uStep;
uniform float dragA, reseed, brown;
uniform uint  uFrame;
` + S_SPECTRUM + `
void main(){
  vec4 P   = texture(uPart, vUv);
  vec2 pos = P.xy;
  vec2 vel = P.zw;

  // gl_FragCoord at a texel centre is (i + 0.5, j + 0.5), so this is the same
  // integer texel the draw pass fetches, and both agree on the grain's size.
  uvec2 id = uvec2(gl_FragCoord.xy);

  vec2 ug = texture(uVel, pos).xy;
  float a = dragA / tauMul(stRank(ivec2(id)));
  vel = (vel + a * ug) / (1.0 + a);           // backward Euler, per-grain tau

  if (brown > 0.0) {
    vel += brown * (vec2(hash1(id, uFrame + 11u), hash1(id, uFrame + 523u)) - 0.5);
  }

  pos = fract(pos + (vel + stream) * uStep);

  // Recycling. A grain used to be reborn at rest, which was invisible when
  // brightness came from absolute speed -- but brightness now comes from drift,
  // and a grain at rest in moving gas has the *largest* drift there is. So every
  // recycled grain flashed on at full brightness and full size before settling.
  // Reborn entrained in the local flow instead: zero drift, dimmest possible,
  // and it is what actually happens to a grain swept into frame.
  if (hash1(id, uFrame) < reseed) {
    pos = vec2(hash1(id, uFrame + 3u), hash1(id, uFrame + 977u));
    vel = texture(uVel, pos).xy;
  }
  outColor = vec4(pos, vel);
}`;

  const F_COMPOSITE = F_HEAD + `
uniform sampler2D uDye;
uniform vec2  uRes;
uniform vec3  uBg, uTint;
uniform float uTime, uGrain, uDyeGain, uVignette;

uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float hash1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}

void main(){
  vec3 raw = texture(uDye, vUv).rgb;
  // Squared density response, the way an optically-thin column actually
  // behaves: faint gas stays dark and only genuinely dense wisps register.
  // A linear response turns every mixed region into a broad pale wash that
  // competes with the dust, and the dust is the subject.
  float dens = clamp(max(max(raw.r, raw.g), raw.b), 0.0, 1.0);
  vec3  dye  = raw * dens * uDyeGain;

  // faint cold gradient so the page never reads as flat black
  float g = smoothstep(1.15, -0.15, vUv.y + vUv.x * 0.22);
  vec3  c = uBg + uTint * g * 0.5;

  c += dye;

  // Roll off on luminance, not per channel: a per-channel Reinhard drives
  // every bright pixel toward white and throws the section accent away.
  float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
  c *= 1.0 / (1.0 + 0.62 * l);

  float d = length((vUv - 0.5) * vec2(uRes.x / uRes.y, 1.0));
  c *= 1.0 - uVignette * smoothstep(0.30, 1.05, d);

  c += (hash1(uvec2(gl_FragCoord.xy), uint(uTime * 60.0)) - 0.5) * uGrain;
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

  // -------------------------------------------------------------- presets

  // tau is the grain stopping time in seconds. The gas turns over on a
  // timescale of order 1 s here, so tau ~ 0.03 is a tightly coupled tracer
  // and tau ~ 1 is a grain that goes ballistic between eddies.
  // Ceiling on |u|, in cells per simulation second.
  const VMAX = 420;

  const BASE = {
    accent: [0.62, 0.72, 0.98],
    tint: [0.020, 0.024, 0.040],
    tau: 0.10,
    dustGain: 1.00,
    dyeGain: 0.75,
    forceGain: 1.00,
    curl: 12.0,
    velDiss: 0.22,
    dyeDiss: 0.55,
    pointSize: 1.5,
    // Drift speed, in cells/second, at which a grain reads as fully decoupled.
    // 0 means "derive it from tau" (see resolvePreset). Chapters that give the
    // dust a bulk stream of its own set it explicitly.
    driftRef: 0,
    brown: 0.0,
    guide: 0.0,
    bhat: [1.0, 0.0],
    stream: 0.0,
    stirGain: 1.00,
    grain: 0.030,
    vignette: 0.62
  };

  // dustGain now scales instantaneous grain brightness rather than the
  // steady state of a trail buffer, where the level went as gain x persistence.
  // The gains below carry each chapter's old persistence folded in, so the
  // relative weight of the dust from chapter to chapter is unchanged.
  const PRESETS = {
    hero:      { accent: [0.50, 0.64, 1.00], tau: 0.14, dustGain: 1.48, curl: 16, stirGain: 1.25 },
    bio:       { accent: [0.96, 0.80, 0.50], tau: 0.10, dustGain: 1.05, curl: 12 },
    picdust:   { accent: [1.00, 0.58, 0.18], tau: 0.030, dustGain: 1.45, curl: 15, pointSize: 1.2, stirGain: 1.15 },
    dfmm:      { accent: [0.26, 0.95, 0.72], tau: 0.85,  dustGain: 2.24, curl: 9,  pointSize: 1.75, brown: 0.05 },
    mhd:       { accent: [0.60, 0.46, 1.00], tau: 0.12,  dustGain: 1.81, curl: 8,  guide: 3.4, bhat: [0.94, 0.34] },
    cosmicray: { accent: [0.22, 0.86, 1.00], tau: 0.020, dustGain: 2.40, curl: 10, stream: 190, bhat: [0.94, 0.34], pointSize: 1.15, driftRef: 260, guide: 2.2 },
    phrike:    { accent: [1.00, 0.24, 0.48], tau: 0.06,  dustGain: 1.38, curl: 22, velDiss: 0.10, stirGain: 1.35 },
    papers:    { accent: [0.94, 0.82, 0.52], tau: 0.09,  dustGain: 0.95, curl: 11, stirGain: 0.85 },
    blog:      { accent: [0.56, 0.72, 1.00], tau: 0.11,  dustGain: 0.95, curl: 10, stirGain: 0.8 },
    iron:      { accent: [1.00, 0.13, 0.10], tau: 0.55,  dustGain: 1.00, curl: 5, velDiss: 0.60, pointSize: 1.6, stirGain: 0.65, vignette: 0.74 },
    cactus:    { accent: [0.54, 0.84, 0.32], tau: 0.95,  dustGain: 2.06, curl: 3, velDiss: 0.42, pointSize: 1.7, brown: 0.09, stirGain: 0.42, grain: 0.040 },
    skate:     { accent: [1.00, 0.48, 0.06], tau: 0.20,  dustGain: 1.19, curl: 17, velDiss: 0.34, pointSize: 1.3, stirGain: 1.1 },
    contact:   { accent: [0.92, 0.86, 0.58], tau: 0.10,  dustGain: 1.05, curl: 12, stirGain: 0.9 }
  };

  function resolvePreset(key) {
    const p = Object.assign({}, BASE, PRESETS[key] || PRESETS.hero);
    // A decoupled grain drifts further from the gas simply because tau is
    // larger, so one fixed saturation scale would peg every heavy-grain chapter
    // at full brightness and throw away exactly the structure the brightness is
    // there to show. Scale the bar with the stopping time instead. Resolved to
    // a concrete number here, never left falsy, so blending between chapters
    // interpolates it smoothly rather than snapping at the end of the fade.
    if (!p.driftRef) p.driftRef = 46 + 95 * p.tau;
    return p;
  }

  // ------------------------------------------------------------- the field

  function Field(canvas) {
    const gl = canvas.getContext('webgl2', {
      alpha: false, depth: false, stencil: false, antialias: false,
      premultipliedAlpha: false, preserveDrawingBuffer: false,
      powerPreference: 'default'
    });
    if (!gl) throw new Error('no webgl2');
    if (!gl.getExtension('EXT_color_buffer_float')) throw new Error('no float render targets');
    gl.getExtension('OES_texture_float_linear');

    const reduced = global.matchMedia &&
      global.matchMedia('(prefers-reduced-motion: reduce)').matches;
    const mobile = Math.min(global.innerWidth, global.innerHeight) < 720;

    // Grains are allocated once, at the top count, and never reallocated: a
    // tier change draws a shorter prefix of the same population instead. Since
    // grain i is scattered independently of i, a prefix is a uniform random
    // subsample -- so quality can move without teleporting a single grain.
    const PART_SIDE = 320;   // 102,400 grains

    // quality tiers: [gridH, dyeScale, unused, grainFrac, dprCap]
    //
    // The dust is the subject, so it is the last thing given up: a tier drops
    // solver resolution and device pixels first and keeps most of the grains.
    // The gas is a smooth field and survives coarsening; the dust is where all
    // the structure lives, and thinning it is immediately visible.
    const TIERS = [
      [192, 1.8, 208, 1.00, 1.35],
      [160, 1.8, 176, 0.85, 1.20],
      [128, 1.5, 144, 0.65, 1.05],
      [ 96, 1.3, 104, 0.45, 1.00]
    ];
    let tier = mobile ? 3 : 1;
    let ceiling = 0;

    // Geometry and programs live in functions, not one-shot consts, because a
    // lost context invalidates every GL object and they all have to be remade.
    let quad = null, vb = null, emptyVAO = null, P = null;

    function buildGeometry() {
      quad = gl.createVertexArray();
      gl.bindVertexArray(quad);
      vb = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, vb);
      gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
      gl.enableVertexAttribArray(0);
      gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
      gl.bindVertexArray(null);
      emptyVAO = gl.createVertexArray();
    }

    function buildPrograms() {
      P = {
        advect:   program(gl, VERT, F_ADVECT),
        diverge:  program(gl, VERT, F_DIVERGENCE),
        jacobi:   program(gl, VERT, F_JACOBI),
        mgrest:   program(gl, VERT, F_MG_RESTRICT),
        mgprol:   program(gl, VERT, F_MG_PROLONG),
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
    }

    // Every render target ever handed out, so the previous set can actually be
    // deleted instead of leaked when the size or the quality tier changes.
    let owned = [];
    let persistent = [];

    function mkFBO(w, h, i, f, t, fil, wr) {
      const o = makeFBO(gl, w, h, i, f, t, fil, wr);
      owned.push(o);
      return o;
    }

    function mkDouble(w, h, i, f, t, fil, wr) {
      let a = mkFBO(w, h, i, f, t, fil, wr);
      let b = mkFBO(w, h, i, f, t, fil, wr);
      return {
        w, h, texel: a.texel,
        get read() { return a; },
        get write() { return b; },
        swap() { const s = a; a = b; b = s; }
      };
    }

    function releaseAll() {
      const all = owned.concat(persistent);
      for (let i = 0; i < all.length; i++) {
        gl.deleteFramebuffer(all[i].fbo);
        gl.deleteTexture(all[i].tex);
      }
      owned = [];
      persistent = [];
      part = null;
      vel = null;
    }

    let vel, dye, prs, div, crl, part, met, red1, red2, red3;
    let mgP = [], mgD = [], mgDim = [];
    let grid = { w: 0, h: 0 }, dyeRes = { w: 0, h: 0 }, pSide = 0, nPart = 0, nDraw = 0;
    let dpr = 1;

    const state = {
      preset: resolvePreset('hero'),
      target: resolvePreset('hero'),
      blend: 1,
      time: 0,
      running: true,
      fps: 30,
      wall: 0,
      diag: { rms: 0, max: 0, div: 0 }
    };

    const pending = { splats: [], shear: 0, blasts: [] };
    const RGBA16 = [gl.RGBA16F, gl.RGBA, gl.HALF_FLOAT];
    const RG16 = [gl.RG16F, gl.RG, gl.HALF_FLOAT];
    const R16 = [gl.R16F, gl.RED, gl.HALF_FLOAT];
    // The pressure hierarchy is 32-bit on purpose. A residual is a difference of
    // terms of size ~4p, so at half-float precision (~3 decimal digits) the
    // cancellation destroys exactly the small quantity multigrid exists to
    // resolve, and the V-cycle stalls on roundoff instead of converging.
    const R32 = [gl.R32F, gl.RED, gl.FLOAT];

    function drawQuad(fbo) {
      if (fbo) { gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.fbo); gl.viewport(0, 0, fbo.w, fbo.h); }
      else { gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.viewport(0, 0, canvas.width, canvas.height); }
      gl.bindVertexArray(quad);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    let allocW = -1, allocH = -1, allocTier = -1, warmed = false;

    function allocate(force) {
      const t = TIERS[tier];
      dpr = Math.min(global.devicePixelRatio || 1, t[4]);
      const cw = Math.max(1, Math.round(canvas.clientWidth * dpr));
      const ch = Math.max(1, Math.round(canvas.clientHeight * dpr));

      // Mobile browsers fire resize on every URL-bar show/hide, which means
      // continuously while scrolling. Rebuilding a full set of float render
      // targets on each one is what exhausted GPU memory and lost the context.
      if (!force && cw === allocW && ch === allocH && tier === allocTier) return;
      allocW = cw; allocH = ch; allocTier = tier;

      // Hold on to the old set until the new one has been filled from it. A
      // tier change used to clear every buffer and teleport every grain, which
      // read exactly like the page crashing and restarting.
      const prev = owned;
      const old = vel ? { vel: vel, dye: dye } : null;
      owned = [];

      canvas.width = cw; canvas.height = ch;

      const aspect = cw / ch;
      grid.h = t[0];
      // The hierarchy halves the grid MG_LEVELS times, so both dimensions have
      // to stay even the whole way down. Every tier's height is already a
      // multiple of 16; rounding the width to one costs at most ~2% of the
      // aspect ratio and buys exact 2:1 grid transfers, which is worth far more
      // than 2% of aspect to the convergence rate.
      const MGQ = 1 << MG_LEVELS;
      grid.w = Math.max(MGQ * 3, Math.round(t[0] * aspect / MGQ) * MGQ);
      dyeRes.h = Math.round(grid.h * t[1]);
      dyeRes.w = Math.round(grid.w * t[1]);

      const R = gl.REPEAT, C = gl.CLAMP_TO_EDGE, L = gl.LINEAR, N = gl.NEAREST;
      vel = mkDouble(grid.w, grid.h, RG16[0], RG16[1], RG16[2], L, R);
      dye = mkDouble(dyeRes.w, dyeRes.h, RGBA16[0], RGBA16[1], RGBA16[2], L, R);
      prs = mkDouble(grid.w, grid.h, R32[0], R32[1], R32[2], L, R);
      div = mkFBO(grid.w, grid.h, R32[0], R32[1], R32[2], L, R);

      // coarse levels of the V-cycle. NEAREST because both transfer operators
      // are written out explicitly; REPEAT because the domain is periodic.
      mgP = []; mgD = [];
      mgDim = [{ w: grid.w, h: grid.h }];
      let mw = grid.w, mh = grid.h;
      for (let l = 1; l <= MG_LEVELS; l++) {
        mw >>= 1; mh >>= 1;
        mgDim.push({ w: mw, h: mh });
        mgP.push(mkDouble(mw, mh, R32[0], R32[1], R32[2], N, R));
        mgD.push(mkFBO(mw, mh, R32[0], R32[1], R32[2], N, R));
      }
      crl = mkFBO(grid.w, grid.h, R16[0], R16[1], R16[2], L, R);
      // 32-bit: |u|^2 overflows half-float once speeds pass ~256 cells/s
      met = mkFBO(grid.w, grid.h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);

      const r1w = Math.max(1, grid.w >> 2), r1h = Math.max(1, grid.h >> 2);
      const r2w = Math.max(1, r1w >> 2), r2h = Math.max(1, r1h >> 2);
      const r3w = Math.max(1, r2w >> 2), r3h = Math.max(1, r2h >> 2);
      red1 = mkFBO(r1w, r1h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);
      red2 = mkFBO(r2w, r2h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);
      red3 = mkFBO(r3w, r3h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);

      // Grains are allocated once for the life of the context. Nothing about
      // them depends on the canvas size, because they are drawn at native
      // resolution straight to the screen and never through a buffer of their
      // own -- which is also what removed two full-canvas float targets and a
      // full-canvas read-modify-write from every frame.
      if (!part) {
        pSide = PART_SIDE;
        nPart = pSide * pSide;
        part = mkDouble(pSide, pSide, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);
        persistent = owned.splice(owned.length - 2, 2);
        seedParticles();
      }
      // How many of them this tier can afford to draw. The update pass always
      // runs the whole population -- 100k fragments is noise next to the grid
      // work -- so a tier change costs nothing but fill.
      nDraw = Math.max(4096, Math.round(nPart * t[3]));

      // Resample the continuous fields into their new resolution rather than
      // starting from black. Velocity is held in cells per second, so a change
      // of grid has to rescale it or the flow speed jumps.
      if (old) {
        copyInto(old.vel.read, vel.read, grid.w / Math.max(1, old.vel.w));
        copyInto(old.dye.read, dye.read, 1);
      }

      for (let i = 0; i < prev.length; i++) {
        gl.deleteFramebuffer(prev[i].fbo);
        gl.deleteTexture(prev[i].tex);
      }

      if (!warmed) { seedField(140); warmed = true; }

      // canvas.width above cleared the drawing buffer; repaint now so the
      // reallocation never shows through as a black frame.
      paint(activePreset());
    }

    // straight blit, filtering handles the change of resolution
    function copyInto(src, dst, scale) {
      gl.useProgram(P.fade.p);
      gl.uniform1i(P.fade.u.uSrc, src.bind(0));
      gl.uniform1f(P.fade.u.amount, scale);
      drawQuad(dst);
    }

    // One tier change is invisible now, but a machine sitting exactly on the
    // boundary could still trade back and forth forever. Give it a budget.
    let tierChanges = 0;
    function retune(delta) {
      if (tierChanges >= 4) return;
      const next = tier + delta;
      if (next < 0 || next > TIERS.length - 1) return;
      tier = next;
      if (delta > 0) ceiling = tier;   // never climb back past a rejected tier
      tierChanges++;
      safeAllocate();
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
    function seedField(steps) {
      const a = state.target.accent;
      for (let i = 0; i < 11; i++) {
        const ph = i * 2.399;
        pending.splats.push({
          x: 0.5 + 0.38 * Math.cos(ph), y: 0.5 + 0.34 * Math.sin(ph * 1.31),
          dx: 90 * Math.sin(ph * 1.7), dy: 90 * Math.cos(ph * 2.3),
          color: [a[0] * 0.030, a[1] * 0.030, a[2] * 0.030], radius: 0.013
        });
      }
      // Long enough to pass a few dissipation times (1/velDiss ~ 4.5 s) so the
      // page opens on a saturated field. Semi-Lagrangian advection is
      // unconditionally stable, so a coarse warm-up step is safe.
      const dt = 1 / 30;
      for (let i = 0; i < steps; i++) {
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
        // Dye gets its own footprint. The ambient stirrers want a broad wash;
        // the pointer wants a tight one, because a wide bloom under the cursor
        // reads as a glow effect rather than as ink going into a fluid.
        gl.uniform1f(P.splat.u.radius, s.dyeRadius || s.radius * 1.5);
        drawQuad(dye.write); dye.swap();
      }
      pending.splats.length = 0;
    }

    function applyShear(pr, dt) {
      if (Math.abs(pending.shear) < 1e-4) return;
      const amp = Math.max(-280, Math.min(280, pending.shear));
      gl.useProgram(P.shear.p);
      gl.uniform1i(P.shear.u.uVel, vel.read.bind(0));
      gl.uniform1f(P.shear.u.amp, amp);
      gl.uniform1f(P.shear.u.width, 0.055);
      gl.uniform1f(P.shear.u.seed, amp * 0.30);
      gl.uniform1f(P.shear.u.phase, state.time * 0.21);
      drawQuad(vel.write); vel.swap();

      const a = pr.accent, s = Math.min(0.95, Math.abs(amp) / 280 * 0.95) * 0.22;
      gl.useProgram(P.band.p);
      gl.uniform1i(P.band.u.uTarget, dye.read.bind(0));
      gl.uniform1f(P.band.u.width, 0.055);
      gl.uniform1f(P.band.u.strength, s);
      gl.uniform3f(P.band.u.color, a[0], a[1], a[2]);
      drawQuad(dye.write); dye.swap();

      // Decay per unit simulation time rather than per frame, so a scroll or a
      // chapter turn drives the same roll-up whatever rate the page is running.
      pending.shear *= Math.pow(0.55, dt / 0.00333);
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
    // The amplitude is an injection *rate*, scaled by the step, so the driving
    // per unit simulation time is independent of frame rate and of TIME_SCALE.
    function stir(sdt, pr) {
      const t = state.time;
      const g = pr.stirGain * 17 * (sdt * 60);
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
          color: [a[0] * 0.0045, a[1] * 0.0045, a[2] * 0.0045],
          radius: 0.024
        });
      }
    }

    // ------------------------------------------------------------ multigrid
    //
    // Jacobi is a smoother, not a solver. Its convergence factor for the
    // Poisson operator is 1 - O((pi/N)^2): it damps error at the grid scale in a
    // few sweeps and error at the domain scale essentially never, so reaching a
    // converged projection needs O(N^2) sweeps -- thousands, at N = 192. Sixteen
    // sweeps were a smoothing pass wearing a solver's name, and the measured
    // divergence wandered with the flow to prove it.
    //
    // A V-cycle instead carries the smooth error down to grids where it *is*
    // short-wavelength, kills it there for almost nothing, and interpolates the
    // correction back. Cost is bounded by the geometric series 1 + 1/4 + 1/16 +
    // ... = 4/3 of the fine-grid smoothing work, so this is both cheaper than 16
    // sweeps and convergent at every wavelength rather than only the short ones.
    // Shape chosen by measurement, not by the textbook. Cost here is very nearly
    // linear in the *number of passes* (~110 us each) and almost independent of
    // the work inside one, because at 336x192 a pass is 64k fragments and this
    // GPU eats that for nothing. So multigrid's classic advantage -- less
    // arithmetic -- buys nothing, and its cost is simply its pass count:
    //
    //     passes = levels * (pre + post + 2) + coarse
    //
    // Measured relative residual ||d - Ap|| / ||d|| from a cold start, one solve:
    //
    //     jacobi x16   16 passes   0.297      <- what this replaced
    //     jacobi x32   32 passes   0.194
    //     V(1,1) x3    12 passes   0.165
    //     V(2,1) x3    16 passes   0.118      <- this
    //     V(2,1) x5    32 passes   0.114
    //
    // Two coarsenings at the same cost as the old sixteen sweeps, for 2.5x less
    // residual; going deeper doubles the pass count to gain almost nothing.
    const MG_LEVELS = 2, MG_PRE = 2, MG_POST = 1, MG_COARSE = 6, MG_OMEGA = 0.8;

    function mgSmooth(Pl, Dl, dim, n, omega) {
      gl.useProgram(P.jacobi.p);
      gl.uniform2f(P.jacobi.u.texel, 1 / dim.w, 1 / dim.h);
      gl.uniform1f(P.jacobi.u.omega, omega);
      gl.uniform1i(P.jacobi.u.uDiv, Dl.bind(1));
      for (let i = 0; i < n; i++) {
        gl.uniform1i(P.jacobi.u.uP, Pl.read.bind(0));
        drawQuad(Pl.write); Pl.swap();
      }
    }

    function mgZero(f) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, f.fbo);
      gl.viewport(0, 0, f.w, f.h);
      gl.clearColor(0, 0, 0, 1);
      gl.clear(gl.COLOR_BUFFER_BIT);
    }

    function vcycle() {
      const Pl = [prs].concat(mgP);
      const Dl = [div].concat(mgD);

      for (let l = 0; l < MG_LEVELS; l++) {
        mgSmooth(Pl[l], Dl[l], mgDim[l], MG_PRE, MG_OMEGA);
        gl.useProgram(P.mgrest.p);
        gl.uniform1i(P.mgrest.u.uP, Pl[l].read.bind(0));
        gl.uniform1i(P.mgrest.u.uDiv, Dl[l].bind(1));
        gl.uniform2i(P.mgrest.u.fineSize, mgDim[l].w, mgDim[l].h);
        drawQuad(Dl[l + 1]);
        // every level below the finest solves for a *correction*, so its
        // initial guess is zero, not last frame's field.
        mgZero(Pl[l + 1].read);
      }

      mgSmooth(Pl[MG_LEVELS], Dl[MG_LEVELS], mgDim[MG_LEVELS], MG_COARSE, 1.0);

      for (let l = MG_LEVELS - 1; l >= 0; l--) {
        gl.useProgram(P.mgprol.p);
        gl.uniform1i(P.mgprol.u.uP, Pl[l].read.bind(0));
        gl.uniform1i(P.mgprol.u.uCoarse, Pl[l + 1].read.bind(1));
        gl.uniform2i(P.mgprol.u.coarseSize, mgDim[l + 1].w, mgDim[l + 1].h);
        drawQuad(Pl[l].write); Pl[l].swap();
        mgSmooth(Pl[l], Dl[l], mgDim[l], MG_POST, MG_OMEGA);
      }
    }

    function step(dt) {
      const pr = activePreset();

      stir(dt, pr);
      applySplats(pr);
      applyShear(pr, dt);
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

      // The finest level keeps last frame's pressure as its initial guess: the
      // flow changes little in 1/60 s, so the warm start is most of the answer.
      vcycle();

      gl.useProgram(P.gradsub.p);
      gl.uniform2f(P.gradsub.u.texel, vel.texel[0], vel.texel[1]);
      gl.uniform1f(P.gradsub.u.vmax, VMAX);
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

      // dust: one backward-Euler drag update per grain. tau varies across the
      // population but not in time, so the median grain's dt/tau is the only
      // scalar the CPU has to supply; each grain rescales it by its own place
      // in the size spectrum.
      gl.useProgram(P.pupdate.p);
      gl.uniform1i(P.pupdate.u.uPart, part.read.bind(0));
      gl.uniform1i(P.pupdate.u.uVel, vel.read.bind(1));
      gl.uniform1f(P.pupdate.u.dragA, dt / Math.max(pr.tau, 1e-4));
      // pre-multiplied drift: velocity (cells/s) -> uv displacement this step
      gl.uniform2f(P.pupdate.u.uStep, dt * vel.texel[0], dt * vel.texel[1]);
      // Recycling rate per unit simulation time. This used to be four times
      // faster, to stop a caustic sharpening indefinitely into a persistent
      // trail buffer. There is no trail buffer any more, only instantaneous
      // grain positions, so the filaments are free to develop for ~17 s of
      // simulation -- many eddy turnovers -- and that is where the structure
      // in the dust comes from. The slow refresh remains only to stop every
      // grain ending up in one attractor and leaving the rest of the frame bare.
      gl.uniform1f(P.pupdate.u.reseed, 0.06 * dt);
      gl.uniform1f(P.pupdate.u.brown, pr.brown * 240 * (dt * 60));
      gl.uniform1ui(P.pupdate.u.uFrame, frameNo);
      const b = pr.bhat, bn = Math.hypot(b[0], b[1]) || 1;
      gl.uniform2f(P.pupdate.u.stream, pr.stream * b[0] / bn, pr.stream * b[1] / bn);
      drawQuad(part.write); part.swap();
    }

    // The gas: one full-screen pass, tone-mapped and vignetted.
    function composite(pr) {
      gl.useProgram(P.comp.p);
      gl.uniform1i(P.comp.u.uDye, dye.read.bind(0));
      gl.uniform2f(P.comp.u.uRes, canvas.width, canvas.height);
      gl.uniform3f(P.comp.u.uBg, 0.014, 0.015, 0.021);
      gl.uniform3f(P.comp.u.uTint, pr.tint[0], pr.tint[1], pr.tint[2]);
      gl.uniform1f(P.comp.u.uTime, state.wall);
      gl.uniform1f(P.comp.u.uGrain, pr.grain);
      gl.uniform1f(P.comp.u.uDyeGain, pr.dyeGain);
      gl.uniform1f(P.comp.u.uVignette, pr.vignette);
      drawQuad(null);
    }

    // The dust: one additive sprite per grain, at native resolution, on top of
    // the gas. Every grain used to be smeared into a persistent trail buffer
    // held at a fraction of the canvas and then blurred back up, which cost two
    // float targets and a read-modify-write per frame and still looked soft.
    // Drawing them where they actually are is both cheaper and sharper.
    function drawDust(pr) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, canvas.width, canvas.height);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE);
      gl.useProgram(P.pdraw.p);
      gl.uniform1i(P.pdraw.u.uPart, part.read.bind(0));
      gl.uniform1i(P.pdraw.u.uVel, vel.read.bind(1));
      gl.uniform1i(P.pdraw.u.uPW, pSide);
      gl.uniform1f(P.pdraw.u.uPointSize, pr.pointSize * dpr * 1.7);
      gl.uniform1f(P.pdraw.u.uDriftNorm, 1 / Math.max(1, pr.driftRef));
      gl.uniform1f(P.pdraw.u.uAlpha, pr.dustGain * 0.30);
      gl.uniform1f(P.pdraw.u.uVignette, pr.vignette);
      gl.uniform1f(P.pdraw.u.uAspect, canvas.width / Math.max(1, canvas.height));
      const s = pr.bhat, sn = Math.hypot(s[0], s[1]) || 1;
      gl.uniform2f(P.pdraw.u.uStream, pr.stream * s[0] / sn, pr.stream * s[1] / sn);
      const a = pr.accent;
      gl.uniform3f(P.pdraw.u.uColor, a[0], a[1], a[2]);
      gl.bindVertexArray(emptyVAO);
      gl.drawArrays(gl.POINTS, 0, nDraw);
      gl.disable(gl.BLEND);
    }

    function paint(pr) { composite(pr); drawDust(pr); }

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

    // The simulation clock runs an order of magnitude slower than wall time.
    // Every rate-like term (driving, Brownian kicks, reseeding) is scaled by the
    // step, so this slows the flow down without changing its character, its
    // Reynolds number, or the balance between driving and dissipation. It also
    // buys a much better Courant number: a few tenths of a cell per step.
    const TIME_SCALE = 0.1;

    // This used to be capped at 30, on the argument that the field is slow
    // enough that a higher rate buys nothing. That was true while a grain was
    // smeared into a trail: the smear covered the gap between frames. A hard
    // little sprite does not, and at 30 fps a hundred thousand of them strobe.
    // So the gate goes to the display rate. Note the physics does not change:
    // the step is scaled by real elapsed time, so the flow evolves at the same
    // rate per second either way.
    const TARGET_FPS = 60;
    const FRAME_MIN = 1 / TARGET_FPS - 0.002;

    // Quality is judged against a floor, not against the target. A machine
    // holding a steady 40 is doing fine and must not be ratcheted down to the
    // cheapest tier merely for missing 60.
    const FPS_FLOOR = 34;

    let last = 0, frames = 0, fpsT = 0, slow = 0, fast = 0, starved = 0, frameNo = 0;
    let frozen = false, lost = false, dead = false;

    function frame(now) {
      if (dead) return;                         // stop asking for frames at all
      requestAnimationFrame(frame);
      if (lost || gl.isContextLost()) return;
      if (!state.running) { last = 0; return; }
      if (!last) { last = now; return; }

      const elapsed = (now - last) / 1000;
      if (elapsed < FRAME_MIN) return;          // hold the frame rate down
      last = now;

      frameNo = (frameNo + 1) >>> 0;
      const dtWall = Math.min(elapsed, 1 / 12);
      const sdt = dtWall * TIME_SCALE;
      state.time += sdt;
      state.wall += dtWall;

      frames++; fpsT += dtWall;
      if (state.wall < 2.5) { frames = 0; fpsT = 0; }
      else if (fpsT > 0.75) {
        state.fps = frames / fpsT;
        frames = 0; fpsT = 0;
        // Measured throughput is the only signal worth trusting, so start in
        // the middle and let it settle both ways rather than opening at full
        // quality and hoping. A tier that has already proved too expensive is
        // never retried, which is what stops the two states oscillating.
        if (state.fps < FPS_FLOOR - 6 && tier < TIERS.length - 1) {
          fast = 0;
          if (++slow >= 3) { slow = 0; retune(+1); }
        } else if (state.fps > FPS_FLOOR + 12 && tier > ceiling) {
          slow = 0; starved = 0;
          if (++fast >= 5) { fast = 0; retune(-1); }
        } else if (state.fps > FPS_FLOOR) {
          slow = 0; starved = 0; fast = 0;
        }
        // Already at the cheapest tier and still crawling: this machine should
        // not be running a fluid solver behind a web page. Stop stepping for
        // good and leave the last field on screen. Compositing continues, both
        // because one full-screen quad is nearly free and because the drawing
        // buffer is not preserved between frames.
        if (tier === TIERS.length - 1 && state.fps < 14 && ++starved >= 3) {
          frozen = true;
        }
      }

      lerpPreset(dtWall);
      const pr = activePreset();

      if (!frozen) {
        step(sdt);
        // readPixels forces a pipeline sync, so keep it rare
        if (--measureCountdown <= 0) { measureCountdown = 150; measure(); }
      }
      // Painted even when frozen: the drawing buffer is not preserved between
      // frames, so skipping this would blank the page rather than hold it.
      paint(pr);

      if (reduced && state.wall > 2.5) frozen = true;
    }

    // ------------------------------------------------------------- public

    // Allocation can fail outright when memory is tight. Step down a tier and
    // try once more; if even that fails, stop stepping rather than throw from
    // inside a resize handler and leave half a solver behind.
    function safeAllocate(force) {
      try {
        allocate(force);
      } catch (err) {
        releaseAll();
        if (tier < TIERS.length - 1) {
          tier++;
          allocW = allocH = allocTier = -1;
          try { allocate(true); return; } catch (e2) { releaseAll(); }
        }
        frozen = true;
        dead = true;
      }
    }

    let resizeTimer = null;
    function onResize() {
      clearTimeout(resizeTimer);
      // allocate() no-ops unless the pixel size really changed, so a scroll-
      // driven storm of resize events costs nothing.
      resizeTimer = setTimeout(function () { safeAllocate(); frozen = false; }, 250);
    }

    // A lost context used to leave a dead canvas. Preventing the default is what
    // makes restoration possible at all; on restore, every GL object is invalid,
    // so programs, geometry and render targets are all rebuilt.
    canvas.addEventListener('webglcontextlost', function (e) {
      e.preventDefault();
      lost = true;
      owned = [];
      persistent = [];
      part = null;
      vel = null;
    }, false);

    canvas.addEventListener('webglcontextrestored', function () {
      try {
        buildGeometry();
        buildPrograms();
        allocW = allocH = allocTier = -1;
        warmed = false;
        allocate(true);
        dead = false;
        lost = false;
        last = 0;
        frozen = false;
      } catch (err) {
        dead = true;
      }
    }, false);

    buildGeometry();
    buildPrograms();
    safeAllocate(true);
    requestAnimationFrame(frame);
    global.addEventListener('resize', onResize);
    document.addEventListener('visibilitychange', function () {
      state.running = !document.hidden;
      last = 0;
    });

    return {
      kind: 'webgl2',
      splat(x, y, dx, dy, radius) {
        const pr = activePreset(), a = pr.accent;
        const sp = Math.hypot(dx, dy);
        if (sp > 380) { dx *= 380 / sp; dy *= 380 / sp; }
        const m = Math.min(1, sp / 380);
        // A narrow splat injects less momentum for the same amplitude, so the
        // dye is concentrated rather than merely dimmer: brightness is raised
        // to keep the stroke reading at the same strength it did when it was
        // spread over four times the area.
        const r = radius || 0.0030;
        pending.splats.push({
          x, y, dx, dy,
          color: [a[0] * (0.14 + m * 0.75), a[1] * (0.14 + m * 0.75), a[2] * (0.14 + m * 0.75)],
          radius: r,
          dyeRadius: r * 1.15
        });
        if (pending.splats.length > 6) pending.splats.splice(0, pending.splats.length - 6);
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
          solver: 'multigrid V(' + MG_PRE + ',' + MG_POST + '), ' + (MG_LEVELS + 1) + ' levels',
          mgLevels: MG_LEVELS + 1,
          timeScale: TIME_SCALE, targetFps: TARGET_FPS, vmax: VMAX,
          drag: 'backward Euler'
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

    // divergence-free by construction: u = curl of a scalar stream function.
    // Amplitudes match the WebGL2 path's slowed clock.
    function flow(x, y, tt) {
      const s = 2.7;
      const psx = Math.sin(s * y + tt * 0.31) * Math.cos(s * 0.7 * x - tt * 0.23);
      const psy = Math.cos(s * x - tt * 0.19) * Math.sin(s * 0.6 * y + tt * 0.27);
      return [0.011 * (psx + 0.5 * psy), 0.011 * (psy - 0.5 * psx) + shear * (y - 0.5) * 0.9];
    }

    const STEP = 1 / 30;
    let lastT = 0;

    function frame(now) {
      requestAnimationFrame(frame);
      if (!lastT) { lastT = now; return; }
      if ((now - lastT) / 1000 < STEP - 0.002) return;
      lastT = now;

      t += STEP * 0.1;                       // same one-tenth clock
      shear *= 0.94;
      ctx.fillStyle = 'rgba(4,5,10,0.038)';  // longer trails at 30 fps
      ctx.fillRect(0, 0, w, h);

      // backward Euler, matching the GPU path
      const a = (STEP * 0.1) / tau;
      const b = 1 / (1 + a);
      const col = `rgba(${(accent[0] * 255) | 0},${(accent[1] * 255) | 0},${(accent[2] * 255) | 0},0.42)`;
      ctx.fillStyle = col;
      for (let i = 0; i < N; i++) {
        const g = flow(px[i], py[i], t);
        vx[i] = (vx[i] + a * g[0]) * b;
        vy[i] = (vy[i] + a * g[1]) * b;
        px[i] += vx[i] * STEP; py[i] += vy[i] * STEP;
        if (px[i] < 0) px[i] += 1; if (px[i] > 1) px[i] -= 1;
        if (py[i] < 0) py[i] += 1; if (py[i] > 1) py[i] -= 1;
        ctx.fillRect(px[i] * w, py[i] * h, dpr, dpr);
      }
      if (reduced && t > 0.4) return;
    }
    requestAnimationFrame(frame);

    return {
      kind: 'canvas2d',
      splat(x, y, dx, dy) { shear += dy / 40000; },
      shear(amount) { shear += amount / 900; },
      blast() {},
      setPreset(key) { const p = resolvePreset(key); accent = p.accent; tau = p.tau; },
      accentOf(key) { return resolvePreset(key).accent; },
      stats() { return { gridW: 0, gridH: 0, nPart: N, fps: 30, tier: 9, tau, dpr, rms: 0, max: 0, div: 0, jacobi: 0, timeScale: 0.1, targetFps: 24, vmax: 0, drag: 'backward Euler' }; }
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

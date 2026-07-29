// gpu-riemann2d.mjs -- run riemann2d-hlld.glsl on the actual GPU and compare it with
// mirror-riemann2d.mjs.
//
// The mirror is token-identical to the GLSL, which says the two texts agree; it does not
// say that GLSL evaluates that text the way JS does. A driver may contract a multiply and
// an add into an fma, reassociate, or hand min/max a different tie-break. So this renders
// the real shader over 4096 corners in one 64x64 pass, reads the EMF back out of an
// RGBA32F attachment and compares it with the mirror's double-precision value. Agreement
// is expected at float32 round-off, not below it: the GPU has 24 bits and the mirror has 53.
//
// Uses the cdp.mjs driver that is already in this directory. Never --virtual-time-budget:
// it hangs on WebGL pages.
//
// node gpu-riemann2d.mjs

import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { riemann2dMirror } from './mirror-riemann2d.mjs';

const HERE = new URL('.', import.meta.url).pathname;
const N = 64;                       // N x N corners per pass
const F = Math.fround;

function mulberry32(a) {
  return function () {
    a |= 0; a = a + 0x6D2B79F5 | 0;
    let t = Math.imul(a ^ a >>> 15, 1 | a);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}
const rng = mulberry32(0xC0FFEE11);
const lg = (lo, hi) => Math.exp(Math.log(lo) + (Math.log(hi) - Math.log(lo)) * rng());
const sgn = () => (rng() < 0.5 ? -1 : 1);

const CS2 = 1.0, GAM = 5 / 3, SMALLR = 1e-3, SMALLC = 1e-3;

// One case per output texel; its four quadrants live at 2*c + (0|1, 0|1) of the inputs,
// with the first index the x side and the second the y side, matching the naming.
function cases(adia) {
  const uU = new Float32Array(2 * N * 2 * N * 4), uB = new Float32Array(2 * N * 2 * N * 4);
  const list = [];
  for (let cy = 0; cy < N; cy++) for (let cx = 0; cx < N; cx++) {
    const q = [];
    for (let k = 0; k < 4; k++) {
      const r = lg(SMALLR * 2, 30);
      const p = adia ? lg(1e-2, 30) : r * CS2;
      const z = rng();
      const A = z < 0.15 ? 0 : sgn() * lg(1e-3, 8);
      const B = z > 0.85 ? 0 : sgn() * lg(1e-3, 8);
      q.push({ r: F(r), p: F(p), A: F(A), B: F(B) });
    }
    // velocities over a range wide enough that every leaf of the selection is visited
    const c = 6.0 * rng();
    for (let k = 0; k < 4; k++) { q[k].u = F(c * (2 * rng() - 1) + 4 * (rng() - 0.5)); q[k].v = F(c * (2 * rng() - 1) + 4 * (rng() - 0.5)); }
    for (let k = 0; k < 4; k++) {
      const jx = k >> 1, jy = k & 1;                     // k = 0 LL, 1 LR, 2 RL, 3 RR
      const i = ((2 * cy + jy) * 2 * N + (2 * cx + jx)) * 4;
      uU[i] = q[k].r; uU[i + 1] = q[k].u; uU[i + 2] = q[k].v; uU[i + 3] = 0.5;
      uB[i] = q[k].A; uB[i + 1] = q[k].B; uB[i + 2] = 0.0; uB[i + 3] = q[k].p;
    }
    list.push(q);
  }
  return { uU, uB, list };
}

const piece = readFileSync(HERE + 'riemann2d-hlld.glsl', 'utf8');
const SHADER = `#version 300 es
precision highp float;
precision highp sampler2D;
out vec4 outColor;
uniform float cs, cs2, smallr, adia, gam, smallc;
uniform sampler2D uU, uB;

struct Prim { vec4 h; vec4 f; };   // h = (rho, u, v, Y),  f = (Bx, By, psi, p)
struct Flux { vec4 h; float bn; float bt; float ps; float en; };
float sound2(float r, float p){ return adia > 0.5 ? gam * p / r : cs2; }

${piece}

Prim q(ivec2 c){ return Prim(texelFetch(uU, c, 0), texelFetch(uB, c, 0)); }
void main(){
  ivec2 b = 2 * ivec2(gl_FragCoord.xy);
  outColor = vec4(riemann2d(q(b + ivec2(0,0)), q(b + ivec2(0,1)),
                            q(b + ivec2(1,0)), q(b + ivec2(1,1)), smallc));
}`;

const iso = cases(false), adi = cases(true);
const b64 = a => Buffer.from(a.buffer, a.byteOffset, a.byteLength).toString('base64');
mkdirSync('/tmp/r2d-gpu', { recursive: true });
writeFileSync('/tmp/r2d-gpu/page.html', '<!doctype html><html><body><canvas id="c"></canvas></body></html>');
writeFileSync('/tmp/r2d-gpu/run.js', `
const N = ${N};
const dec = s => new Float32Array(Uint8Array.from(atob(s), c => c.charCodeAt(0)).buffer);
const DATA = { iso: { uU: dec("${b64(iso.uU)}"), uB: dec("${b64(iso.uB)}") },
               adi: { uU: dec("${b64(adi.uU)}"), uB: dec("${b64(adi.uB)}") } };
const src = ${JSON.stringify(SHADER)};
const gl = document.getElementById('c').getContext('webgl2');
if (!gl) return { error: 'no webgl2' };
if (!gl.getExtension('EXT_color_buffer_float')) return { error: 'no float attachments' };
const sh = (t, s) => { const o = gl.createShader(t); gl.shaderSource(o, s); gl.compileShader(o);
  if (!gl.getShaderParameter(o, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(o)); return o; };
const vs = sh(gl.VERTEX_SHADER, '#version 300 es\\nin vec2 aPos;\\nvoid main(){ gl_Position = vec4(aPos,0.0,1.0); }');
const fs = sh(gl.FRAGMENT_SHADER, src);
const prog = gl.createProgram(); gl.attachShader(prog, vs); gl.attachShader(prog, fs); gl.linkProgram(prog);
if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) return { error: 'link: ' + gl.getProgramInfoLog(prog) };
gl.useProgram(prog);
const buf = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, buf);
gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1,-1, 3,-1, -1,3]), gl.STATIC_DRAW);
const loc = gl.getAttribLocation(prog, 'aPos'); gl.enableVertexAttribArray(loc);
gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
function tex(unit, data) {
  const t = gl.createTexture(); gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, t);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, 2*N, 2*N, 0, gl.RGBA, gl.FLOAT, data);
  return t;
}
const out = gl.createTexture(); gl.bindTexture(gl.TEXTURE_2D, out);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, N, N, 0, gl.RGBA, gl.FLOAT, null);
const fbo = gl.createFramebuffer(); gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, out, 0);
gl.viewport(0, 0, N, N);
const res = {};
for (const mode of ['iso', 'adi']) {
  tex(0, DATA[mode].uU); tex(1, DATA[mode].uB);
  gl.uniform1i(gl.getUniformLocation(prog, 'uU'), 0);
  gl.uniform1i(gl.getUniformLocation(prog, 'uB'), 1);
  gl.uniform1f(gl.getUniformLocation(prog, 'cs'), 1.0);
  gl.uniform1f(gl.getUniformLocation(prog, 'cs2'), ${CS2});
  gl.uniform1f(gl.getUniformLocation(prog, 'smallr'), ${SMALLR});
  gl.uniform1f(gl.getUniformLocation(prog, 'smallc'), ${SMALLC});
  gl.uniform1f(gl.getUniformLocation(prog, 'gam'), ${GAM});
  gl.uniform1f(gl.getUniformLocation(prog, 'adia'), mode === 'adi' ? 1.0 : 0.0);
  gl.drawArrays(gl.TRIANGLES, 0, 3);
  const px = new Float32Array(N * N * 4);
  gl.readPixels(0, 0, N, N, gl.RGBA, gl.FLOAT, px);
  // base64, not an array of numbers: the CDP transport truncates a returned JSON
  // value at 64 KB and 4096 doubles spelled out do not fit.
  const one = new Float32Array(N * N);
  for (let i = 0; i < N * N; i++) one[i] = px[i * 4];
  const b = new Uint8Array(one.buffer);
  let str = ''; for (let i = 0; i < b.length; i++) str += String.fromCharCode(b[i]);
  res[mode] = btoa(str);
}
res.renderer = gl.getParameter(gl.getExtension('WEBGL_debug_renderer_info').UNMASKED_RENDERER_WEBGL);
return res;
`);

const raw = execFileSync('node', [HERE + 'cdp.mjs', '--url', 'file:///tmp/r2d-gpu/page.html',
  '--wait', '1500', '--eval-file', '/tmp/r2d-gpu/run.js'], { encoding: 'utf8', maxBuffer: 1 << 28 });
const got = JSON.parse(raw).result;
if (!got || got.error) { console.log('GPU pass failed: ' + JSON.stringify(got)); process.exit(1); }
for (const m of ['iso', 'adi']) got[m] = new Float32Array(Buffer.from(got[m], 'base64').buffer);

const EPS32 = 1.1920929e-7;
console.log('riemann2d on the GPU vs the mirror');
console.log(got.renderer);
console.log('\n' + 'mode'.padEnd(8) + 'n'.padEnd(8) + 'max |dE|'.padEnd(12) + 'max rel(scale)'.padEnd(16) + 'in float32 eps');
let worstAll = 0;
for (const [mode, list] of [['iso', iso.list], ['adi', adi.list]]) {
  const U = { adia: mode === 'adi' ? 1 : 0, gam: GAM, cs2: CS2, smallr: SMALLR };
  let maxAbs = 0, maxRel = 0, worst = null;
  for (let i = 0; i < list.length; i++) {
    const q = list[i];
    const P = q.map(x => ({ h: { x: x.r, y: x.u, z: x.v, w: 0.5 }, f: { x: x.A, y: x.B, z: 0, w: x.p } }));
    const M = riemann2dMirror(P[0], P[1], P[2], P[3], SMALLC, U);
    // The natural size of the answer, not just of the inputs: in the star region ustar
    // and vstar can exceed every quadrant velocity, so E can be several times the largest
    // input EMF, and normalising by the inputs alone overstates the error.
    const scale = Math.max(...q.map(x => Math.abs(x.u * x.B - x.v * x.A)), Math.abs(M.E), 1e-30);
    const d = Math.abs(got[mode][i] - M.E);
    if (d > maxAbs) maxAbs = d;
    if (d / scale > maxRel) { maxRel = d / scale; worst = { i, gpu: got[mode][i], mirror: M.E, branch: M.branch, scale }; }
  }
  console.log(mode.padEnd(8) + String(list.length).padEnd(8) + maxAbs.toExponential(3).padEnd(12) +
    maxRel.toExponential(3).padEnd(16) + (maxRel / EPS32).toFixed(1));
  if (worst) console.log('   worst: case ' + worst.i + ' leaf ' + worst.branch + '  gpu ' + worst.gpu + '  mirror ' + worst.mirror);
  worstAll = Math.max(worstAll, maxRel);
}
const ok = worstAll < 32 * EPS32;
console.log('\n' + (ok ? 'PASS' : 'FAIL') + '  (threshold 32 x float32 eps = ' + (32 * EPS32).toExponential(2) + ')');
process.exit(ok ? 0 : 1);

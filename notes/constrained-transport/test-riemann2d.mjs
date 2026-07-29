// test-riemann2d.mjs -- mirror (my GLSL, transcribed) vs ref (the Fortran, transcribed).
//
// Both sides get bit-identical inputs, rounded to float32 first because that is what a
// texture can actually hold. Both sides then compute in JS doubles, so every difference
// reported below is an ALGEBRAIC difference between my port and the Fortran -- a guard
// that bound, a reassociated product, a dropped term -- and not a precision artefact.
// That is a sharper test than emulating float32 throughout would be: anything at or
// below ~1e-16 relative is reassociation only, and float32 round-off (eps = 1.2e-7) is
// nine orders of magnitude coarser.
//
// node test-riemann2d.mjs

import { riemann2dRef, findSpeedFastRef } from './ref-riemann2d.mjs';
import { riemann2dMirror } from './mirror-riemann2d.mjs';
import { readFileSync } from 'node:fs';

const F = Math.fround;

// ---------------------------------------------------------------- transcription check
// The mirror is only worth testing if it really is the GLSL. Rather than assert that,
// check it: pull the three function bodies out of both files, normalise away the
// differences that are forced by the language (Math., const/float, the uniform struct
// the shader gets for free, the object the mirror has to return to report its branch)
// and compare what is left character for character. Any real edit to one file and not
// the other shows up here as a mismatch.
function body(src, header) {
  const i = src.indexOf(header);
  if (i < 0) return null;
  let j = src.indexOf('{', i), depth = 0, k = j;
  for (; k < src.length; k++) {
    if (src[k] === '{') depth++;
    else if (src[k] === '}') { depth--; if (depth === 0) break; }
  }
  return src.slice(j + 1, k);
}
function normalise(s) {
  return s
    .replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/[^\n]*/g, '')
    .replace(/return\s*\{\s*E:\s*([\s\S]*?),\s*branch:\s*\d+,\s*diag\s*\};/g, 'return $1;')
    .replace(/(const|let)\s+diag\s*=\s*\{[^}]*\};/g, '')
    .replace(/\bMath\./g, '').replace(/\b(const|float)\s+/g, '')
    .replace(/\bU\./g, '').replace(/,\s*U\)/g, ')').replace(/,\s*U,/g, ',')
    .replace(/\s+/g, '');
}
const glslSrc = readFileSync(new URL('./riemann2d-hlld.glsl', import.meta.url), 'utf8');
const mirrorSrc = readFileSync(new URL('./mirror-riemann2d.mjs', import.meta.url), 'utf8');
const PAIRS = [
  ['float cfast2(', 'function cfast2('],
  ['float irs(', 'function irs('],
  ['float riemann2d(', 'export function riemann2dMirror(']
];
const transcript = [];
for (const [g, m] of PAIRS) {
  const a = normalise(body(glslSrc, g) ?? '#missing-glsl'), b = normalise(body(mirrorSrc, m) ?? '#missing-js');
  let at = -1;
  for (let i = 0; i < Math.max(a.length, b.length); i++) if (a[i] !== b[i]) { at = i; break; }
  transcript.push({ name: g.replace('float ', '').replace('(', ''), same: at < 0, at, a, b });
}

// ---------------------------------------------------------------- rng
function mulberry32(a) {
  return function () {
    a |= 0; a = a + 0x6D2B79F5 | 0;
    let t = Math.imul(a ^ a >>> 15, 1 | a);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}
const rng = mulberry32(0x1D0FF1CE);
const lg = (lo, hi) => Math.exp(Math.log(lo) + (Math.log(hi) - Math.log(lo)) * rng());
const sgn = () => (rng() < 0.5 ? -1 : 1);

// ---------------------------------------------------------------- state building
// Field modes. 'zeroA' puts A exactly 0 in all four quadrants (the degenerate normal
// field for the x pair), 'zeroB' the same for y, 'zeroAll' kills the field entirely,
// 'someZeroA' zeroes A in a random subset so the quadrants disagree about it.
const FIELD_MODES = ['generic', 'strongA', 'strongB', 'weakA', 'zeroA', 'zeroB', 'zeroAll', 'someZeroA', 'huge'];

function field(mode, k) {
  switch (mode) {
    case 'strongA': return [sgn() * lg(3, 30), sgn() * lg(0.01, 0.3)];
    case 'strongB': return [sgn() * lg(0.01, 0.3), sgn() * lg(3, 30)];
    case 'weakA': return [sgn() * lg(1e-7, 1e-3), sgn() * lg(0.1, 3)];
    case 'zeroA': return [0, sgn() * lg(0.05, 3)];
    case 'zeroB': return [sgn() * lg(0.05, 3), 0];
    case 'zeroAll': return [0, 0];
    case 'someZeroA': return [(k & 1) ? 0 : sgn() * lg(0.05, 3), (k & 2) ? 0 : sgn() * lg(0.05, 3)];
    case 'huge': return [sgn() * lg(1, 100), sgn() * lg(1, 100)];
    default: return [sgn() * lg(0.05, 3), sgn() * lg(0.05, 3)];
  }
}

// x target: 'R' = supersonic to the right (SL>0), 'L' = supersonic to the left (SR<0),
// 'S' = the fan straddles the corner. Same three in y, where 'R' means upward (SB>0)
// and 'L' downward (ST<0). The fast speeds do not depend on velocity, so the target is
// hit by construction: SL>0 iff every u exceeds cx.
function velocities(target, c) {
  const u = [];
  for (let k = 0; k < 4; k++) {
    if (target === 'R') u.push(c * (1.05 + 2.0 * rng()));
    else if (target === 'L') u.push(-c * (1.05 + 2.0 * rng()));
    else u.push(c * 0.9 * (2 * rng() - 1));
  }
  return u;
}

function build(opts) {
  const adia = opts.adia, cs2 = 1.0, gam = 5 / 3;
  const U = { adia: adia ? 1 : 0, gam, cs2, smallr: 1e-3 };
  const params = { adia: U.adia, gam, cs2, smallc: opts.smallc };
  const s = [];
  for (let k = 0; k < 4; k++) {
    const r = lg(opts.rlo, opts.rhi);
    const p = adia ? lg(opts.plo, opts.phi) : r * cs2;
    const [A, B] = field(opts.fieldMode, k);
    s.push({ r: F(r), p: F(p), A: F(A), B: F(B) });
  }
  const cx = Math.max(...s.map(q => findSpeedFastRef(q.r, q.p, q.A, q.B, params)));
  const cy = Math.max(...s.map(q => findSpeedFastRef(q.r, q.p, q.B, q.A, params)));
  const u = velocities(opts.tx, cx), v = velocities(opts.ty, cy);
  for (let k = 0; k < 4; k++) { s[k].u = F(u[k]); s[k].v = F(v[k]); }
  return { s, U, params };
}

// The near-degenerate family. SAR-SAL >= max(calfvenL, calfvenR) >= smallc, and
// calfven >= |A|/sqrt(rho*), so driving SAR-SAL towards zero needs three things at
// once: a tiny smallc, a vanishing normal field, and |ustar| below the resulting
// Alfven speed. ustar carries the transverse total-pressure difference, so the four
// quadrants have to agree to within that tolerance too -- hence one base state with
// tiny perturbations rather than four independent ones.
function buildDegenerate(opts) {
  const adia = opts.adia, cs2 = 1.0, gam = 5 / 3;
  const U = { adia: adia ? 1 : 0, gam, cs2, smallr: 1e-3 };
  const smallc = lg(1e-9, 1e-6);
  const params = { adia: U.adia, gam, cs2, smallc };
  const kind = opts.kind;                       // 'x', 'y' or 'both'
  const r0 = lg(0.3, 3), p0 = adia ? lg(0.3, 3) : r0 * cs2;
  const eps = lg(1e-7, 1e-4);                   // quadrant-to-quadrant disagreement
  const tinyA = kind === 'y' ? sgn() * lg(0.1, 2) : (rng() < 0.5 ? 0 : sgn() * lg(1e-9, 1e-6));
  const tinyB = kind === 'x' ? sgn() * lg(0.1, 2) : (rng() < 0.5 ? 0 : sgn() * lg(1e-9, 1e-6));
  const s = [];
  for (let k = 0; k < 4; k++) {
    s.push({
      r: F(r0 * (1 + eps * (2 * rng() - 1))),
      p: F(p0 * (1 + eps * (2 * rng() - 1))),
      A: F(tinyA * (1 + eps * (2 * rng() - 1))),
      B: F(tinyB * (1 + eps * (2 * rng() - 1)))
    });
  }
  const cx = Math.max(...s.map(q => findSpeedFastRef(q.r, q.p, q.A, q.B, params)));
  const cy = Math.max(...s.map(q => findSpeedFastRef(q.r, q.p, q.B, q.A, params)));
  for (let k = 0; k < 4; k++) {
    // a velocity that vanishes on the scale of the Alfven speed in the degenerate
    // direction, and an ordinary subsonic one in the other
    s[k].u = F(kind === 'y' ? cx * 0.9 * (2 * rng() - 1) : smallc * (2 * rng() - 1));
    s[k].v = F(kind === 'x' ? cy * 0.9 * (2 * rng() - 1) : smallc * (2 * rng() - 1));
  }
  return { s, U, params };
}

// ---------------------------------------------------------------- comparison
const toPrim = q => ({ h: { x: q.r, y: q.u, z: q.v, w: 0.5 }, f: { x: q.A, y: q.B, z: 0.0, w: q.p } });

const stats = new Map();
function bucket(name) {
  if (!stats.has(name)) stats.set(name, { n: 0, maxAbs: 0, maxRelScale: 0, maxRelRef: 0, worst: null });
  return stats.get(name);
}
const branchCount = new Array(10).fill(0);
let nBranchMismatch = 0, nBranchMismatchFloored = 0, nMirrorBad = 0, nRefBad = 0, nTotal = 0;
let worstBranchMismatch = null, worstMirrorBad = null;

// How close the campaign drives each guarded quantity to the value at which the guard
// would take over. A guard that never binds is a guard that never moved an answer, and
// that is the claim being made about all of them on physical input; the hostile set
// below is where they do bind.
const near = {
  'min rho    (floor smallr)': Infinity,
  'SL-ustar   (floor -1e-8)': Infinity, 'SR-ustar   (floor +1e-8)': Infinity,
  'SB-vstar   (floor -1e-8)': Infinity, 'ST-vstar   (floor +1e-8)': Infinity,
  'min rho*   (floor smallr)': Infinity, 'SAR-SAL    (floor smallc)': Infinity,
  'SAT-SAB    (floor smallc)': Infinity, 'sum rc x   (floor 1e-20)': Infinity,
  'sum rc y   (floor 1e-20)': Infinity
};
let nGuardBound = 0;
// The same four contact distances measured against the fast speed instead of against
// the floor. This is the scale-free question: is (SL - ustar) small *compared with the
// problem*, or only small in absolute terms? If the ratio is bounded away from zero the
// floor is unreachable except on a state with no dynamics in it at all, and its exact
// value between 1e-8 and that bound is unobservable.
const scaleFree = { 'SL-ustar / cx': Infinity, 'SR-ustar / cx': Infinity, 'SB-vstar / cy': Infinity, 'ST-vstar / cy': Infinity };

function run(name, built, expectBranch) {
  const { s, U, params } = built;
  const [LL, LR, RL, RR] = s.map(toPrim);
  const R = riemann2dRef(s[0], s[1], s[2], s[3], params);
  const M = riemann2dMirror(LL, LR, RL, RR, params.smallc, U);
  nTotal++;
  branchCount[M.branch]++;
  if (!Number.isFinite(M.E)) { nMirrorBad++; worstMirrorBad = worstMirrorBad || { name, s, params }; }
  const g = {
    // the engine's prim() has already floored rho before a Prim exists, so the mirror's
    // clamp is redundant there -- but the oracle takes the raw value, so a sample below
    // the floor is a sample where the two are not solving the same problem
    'min rho    (floor smallr)': Math.min(s[0].r, s[1].r, s[2].r, s[3].r) / U.smallr,
    'SL-ustar   (floor -1e-8)': -(R.SL - R.ustar), 'SR-ustar   (floor +1e-8)': R.SR - R.ustar,
    'SB-vstar   (floor -1e-8)': -(R.SB - R.vstar), 'ST-vstar   (floor +1e-8)': R.ST - R.vstar,
    'min rho*   (floor smallr)': R.minRstar / U.smallr, 'SAR-SAL    (floor smallc)': (R.SAR - R.SAL) / Math.max(params.smallc, 1e-20),
    'SAT-SAB    (floor smallc)': (R.SAT - R.SAB) / Math.max(params.smallc, 1e-20),
    'sum rc x   (floor 1e-20)': R.sumRcx / 1e-20, 'sum rc y   (floor 1e-20)': R.sumRcy / 1e-20
  };
  let bound = false;
  for (const k in g) {
    const x = k.startsWith('SL') || k.startsWith('SR') || k.startsWith('SB') || k.startsWith('ST') ? g[k] / 1e-8 : g[k];
    if (Number.isFinite(x)) { if (x < near[k]) near[k] = x; if (x <= 1) bound = true; }
  }
  const sf = {
    'SL-ustar / cx': -(R.SL - R.ustar) / M.diag.cx, 'SR-ustar / cx': (R.SR - R.ustar) / M.diag.cx,
    'SB-vstar / cy': -(R.SB - R.vstar) / M.diag.cy, 'ST-vstar / cy': (R.ST - R.vstar) / M.diag.cy
  };
  for (const k in sf) if (Number.isFinite(sf[k]) && sf[k] < scaleFree[k]) scaleFree[k] = sf[k];
  if (bound) { nGuardBound++; name = 'z [a guard bound] ' + name.slice(0, 1); }
  // The branch is chosen from SL, SR, SB, ST alone, so it can be checked even where the
  // oracle's unguarded arithmetic has already produced a NaN downstream. Where a floor
  // bound, the two sides are not solving the same problem and the branch may legitimately
  // differ -- those are counted apart.
  if (Number.isFinite(R.SL) && Number.isFinite(R.SR) && Number.isFinite(R.SB) && Number.isFinite(R.ST)) {
    const wrong = R.branch !== M.branch || (expectBranch !== undefined && M.branch !== expectBranch);
    if (wrong && bound) nBranchMismatchFloored++;
    else if (wrong) {
      nBranchMismatch++;
      worstBranchMismatch = worstBranchMismatch || { name, s, params, r: R.branch, m: M.branch, want: expectBranch };
    }
  }
  if (!Number.isFinite(R.E)) { nRefBad++; return; }
  const scale = Math.max(
    Math.abs(s[0].u * s[0].B - s[0].v * s[0].A), Math.abs(s[1].u * s[1].B - s[1].v * s[1].A),
    Math.abs(s[2].u * s[2].B - s[2].v * s[2].A), Math.abs(s[3].u * s[3].B - s[3].v * s[3].A));
  const d = Math.abs(M.E - R.E);
  for (const b of [bucket(name), bucket('ALL')]) {
    b.n++;
    if (d > b.maxAbs) b.maxAbs = d;
    if (scale > 0 && d / scale > b.maxRelScale) { b.maxRelScale = d / scale; b.worst = { s, params, ref: R.E, mir: M.E, branch: M.branch }; }
    if (R.E !== 0 && d / Math.abs(R.E) > b.maxRelRef) b.maxRelRef = d / Math.abs(R.E);
  }
}

// ---------------------------------------------------------------- the campaign
const T = ['R', 'L', 'S'];
const pick = a => a[Math.floor(rng() * a.length)];
const base = () => ({ adia: rng() < 0.5, smallc: rng() < 0.7 ? 1e-3 : lg(1e-6, 1e-1), rlo: 0.05, rhi: 20, plo: 0.05, phi: 20 });

// (a) generic subsonic: both fans straddle the corner, ordinary fields
for (let i = 0; i < 6000; i++) run('a generic subsonic', build({ ...base(), fieldMode: 'generic', tx: 'S', ty: 'S' }));

// (b) field extremes: strong normal field, weak field, and A exactly zero
for (let i = 0; i < 6000; i++) {
  const m = FIELD_MODES[i % FIELD_MODES.length];
  run('b field extremes', build({ ...base(), fieldMode: m, tx: pick(T), ty: pick(T) }));
}

// (c) supersonic corners: in x only, in y only, in both
for (let i = 0; i < 6000; i++) {
  const k = i % 3;
  const tx = k === 1 ? 'S' : pick(['R', 'L']);
  const ty = k === 0 ? 'S' : pick(['R', 'L']);
  const label = k === 0 ? 'c supersonic in x' : k === 1 ? 'c supersonic in y' : 'c supersonic in both';
  run(label, build({ ...base(), fieldMode: pick(FIELD_MODES), tx, ty }));
}

// (d) every leaf of the selection, driven on purpose
const LEAF = { 1: ['R', 'R'], 2: ['L', 'R'], 3: ['S', 'R'], 4: ['R', 'L'], 5: ['L', 'L'], 6: ['S', 'L'], 7: ['R', 'S'], 8: ['L', 'S'], 9: ['S', 'S'] };
for (let i = 0; i < 3600; i++) {
  const b = 1 + (i % 9);
  const [tx, ty] = LEAF[b];
  run('d leaf ' + b, build({ ...base(), fieldMode: pick(FIELD_MODES), tx, ty }), b);
}

// (e) near-degenerate Alfven pairs
for (let i = 0; i < 4000; i++) {
  const kind = ['x', 'y', 'both'][i % 3];
  run('e near-degenerate ' + kind, buildDegenerate({ adia: rng() < 0.5, kind }));
}

// (f) exact boundaries of the selection. Every one of the five comparisons in the tree
// is strict, so a state that puts a wave speed exactly on zero distinguishes > from >=
// -- and those states are reachable: ST = max(v) + cy is exactly zero when the topmost
// v is exactly -cy, and cy does not depend on velocity at all, so it can be computed
// first and then cancelled by construction. These are the only samples not rounded to
// float32, because the cancellation has to be exact in whatever precision is in use.
for (let i = 0; i < 1200; i++) {
  const adia = rng() < 0.5, cs2 = 1.0, gam = 5 / 3;
  const U = { adia: adia ? 1 : 0, gam, cs2, smallr: 1e-3 };
  const params = { adia: U.adia, gam, cs2, smallc: 1e-3 };
  const s = [];
  for (let k = 0; k < 4; k++) {
    const r = lg(0.3, 3);
    s.push({ r, p: adia ? lg(0.3, 3) : r * cs2, ...(([A, B]) => ({ A, B }))(field(pick(FIELD_MODES), k)) });
  }
  const cx = Math.max(...s.map(q => findSpeedFastRef(q.r, q.p, q.A, q.B, params)));
  const cy = Math.max(...s.map(q => findSpeedFastRef(q.r, q.p, q.B, q.A, params)));
  const kind = i % 4;
  // one quadrant sits exactly on the edge, the other three strictly inside it
  const edge = Math.floor(rng() * 4);
  for (let k = 0; k < 4; k++) {
    const off = k === edge ? 0 : lg(0.01, 2);
    if (kind === 0) { s[k].v = -cy - off; s[k].u = cx * 0.9 * (2 * rng() - 1); }        // ST = 0
    else if (kind === 1) { s[k].v = cy + off; s[k].u = cx * 0.9 * (2 * rng() - 1); }    // SB = 0
    else if (kind === 2) { s[k].u = cx + off; s[k].v = cy * 0.9 * (2 * rng() - 1); }    // SL = 0
    else { s[k].u = -cx - off; s[k].v = cy * 0.9 * (2 * rng() - 1); }                   // SR = 0
  }
  run('f exact boundary ' + ['ST=0', 'SB=0', 'SL=0', 'SR=0'][kind], { s, U, params });
}

// (g) wide dynamic range: densities and pressures over ten decades, every field mode
// and every velocity target, which is where the guarded quantities come closest to
// their floors without an unphysical state to help them.
for (let i = 0; i < 6000; i++) {
  run('g wide dynamic range', build({
    adia: rng() < 0.5, smallc: rng() < 0.5 ? 1e-3 : lg(1e-6, 1e-1),
    rlo: 1e-4, rhi: 1e4, plo: 1e-6, phi: 1e6,
    fieldMode: pick(FIELD_MODES), tx: pick(T), ty: pick(T)
  }));
}

// (h) hostile states, mirror only. These are the states the Fortran would trap on:
// zero and negative density, zero and negative pressure, a pressure imbalance large
// enough to push the contact outside the fan so that (SL - ustar) changes sign, and
// four identical states at rest so that every rc vanishes at once. There is no oracle
// for these -- the Fortran's answer is a floating-point exception -- so the only
// question asked is whether the guards hold, which is the question that matters,
// because one non-finite texel takes the whole frame with it.
let nHostile = 0, nHostileRefBad = 0, nHostileMirrorBad = 0, worstHostile = null;
const HOSTILE = [
  () => ({ r: 0, p: 0, A: 0, B: 0, u: 0, v: 0 }),
  () => ({ r: 0, p: -1, A: sgn() * lg(0.1, 10), B: sgn() * lg(0.1, 10), u: 0, v: 0 }),
  () => ({ r: -lg(1e-6, 1), p: -lg(1e-6, 1), A: sgn() * lg(1e-6, 10), B: sgn() * lg(1e-6, 10), u: sgn() * lg(1e-6, 10), v: sgn() * lg(1e-6, 10) }),
  () => ({ r: lg(1e-6, 1e-3), p: lg(1e6, 1e9), A: sgn() * lg(1e-6, 1e3), B: sgn() * lg(1e-6, 1e3), u: sgn() * lg(1e-3, 1e3), v: sgn() * lg(1e-3, 1e3) }),
  () => ({ r: lg(1e3, 1e6), p: lg(1e-9, 1e-6), A: sgn() * lg(1e-9, 1e-6), B: sgn() * lg(1e-9, 1e-6), u: sgn() * lg(1e-6, 1e-3), v: sgn() * lg(1e-6, 1e-3) }),
  () => ({ r: lg(0.5, 2), p: lg(0.5, 2), A: 0, B: 0, u: 0, v: 0 }),
  // The configuration that reaches the max() on cfast2's sound speed: A exactly zero, so
  // the inner root is sqrt(d2*d2), and a pressure negative enough to make d2 negative too.
  // Without the max the outer root's argument would be d2 + |d2|. It is worth knowing
  // that this is still safe -- IEEE gives sqrt(x*x) = |x| exactly, so the argument is
  // exactly zero and not one ulp below it -- which is why removing that max is not
  // detectable here. It is kept in the shader as a value guard, not a finiteness one.
  () => { const r = lg(0.1, 10), B = sgn() * lg(0.1, 10);
          return { r, p: -(B * B) / (5 / 3) * lg(1.0, 1.5), A: 0, B, u: sgn() * lg(0.1, 10), v: sgn() * lg(0.1, 10) }; }
];
for (let i = 0; i < 3000; i++) {
  const adia = rng() < 0.5;
  const U = { adia: adia ? 1 : 0, gam: 5 / 3, cs2: 1.0, smallr: 1e-3 };
  const params = { adia: U.adia, gam: 5 / 3, cs2: 1.0, smallc: rng() < 0.5 ? 1e-3 : 0.0 };
  const same = rng() < 0.25, gen = HOSTILE[i % HOSTILE.length];
  const g0 = gen();
  const s = [];
  for (let k = 0; k < 4; k++) {
    const q = same ? { ...g0 } : gen();
    s.push({ r: F(q.r), p: F(q.p), A: F(q.A), B: F(q.B), u: F(q.u), v: F(q.v) });
  }
  const [LL, LR, RL, RR] = s.map(toPrim);
  const M = riemann2dMirror(LL, LR, RL, RR, params.smallc, U);
  const R = riemann2dRef(s[0], s[1], s[2], s[3], params);
  nHostile++;
  if (!Number.isFinite(R.E)) nHostileRefBad++;
  if (!Number.isFinite(M.E) || Math.abs(M.E) > 3.4e38) {
    nHostileMirrorBad++;
    worstHostile = worstHostile || { s, params, E: M.E, branch: M.branch };
  }
}

// (i) transposition antisymmetry, the mirror against itself. Exchange the two
// directions -- u with v, Bx with By, and the LR quadrant with RL, since the first
// letter of a quadrant is its x side and the second its y side -- and Ez = u By - v Bx
// must change sign and nothing else. No oracle is involved, which is the point: this
// tests the quadrant naming and the LL-RL / LR-RR pressure pairing on their own, and one
// misplaced slot anywhere in the function breaks it.
const Tq = q => ({ h: { x: q.h.x, y: q.h.z, z: q.h.y, w: q.h.w }, f: { x: q.f.y, y: q.f.x, z: q.f.z, w: q.f.w } });
let symAbs = 0, symRel = 0, nSym = 0;
for (let i = 0; i < 20000; i++) {
  const adia = rng() < 0.5;
  const U = { adia: adia ? 1 : 0, gam: 5 / 3, cs2: 1.0, smallr: 1e-3 };
  const mk = () => { const r = F(lg(0.05, 20)); return { h: { x: r, y: F(8 * (rng() - 0.5)), z: F(8 * (rng() - 0.5)), w: 0.5 },
    f: { x: F(sgn() * lg(1e-3, 5)), y: F(sgn() * lg(1e-3, 5)), z: 0, w: adia ? F(lg(0.05, 20)) : r } }; };
  const [LL, LR, RL, RR] = [mk(), mk(), mk(), mk()];
  const a = riemann2dMirror(LL, LR, RL, RR, 1e-3, U).E;
  const b = riemann2dMirror(Tq(LL), Tq(RL), Tq(LR), Tq(RR), 1e-3, U).E;
  const d = Math.abs(a + b), sc2 = Math.max(Math.abs(a), Math.abs(b));
  nSym++;
  if (d > symAbs) symAbs = d;
  if (sc2 > 0 && d / sc2 > symRel) symRel = d / sc2;
}

// ---------------------------------------------------------------- report
const EPS32 = 1.1920929e-7;
const pad = (s, n) => String(s).padEnd(n);
const num = x => x.toExponential(3);

console.log('riemann2d: GLSL mirror vs Fortran oracle');
console.log('samples ' + nTotal + ', inputs rounded to float32, both sides evaluated in doubles\n');

console.log(pad('category', 26) + pad('n', 8) + pad('max |dE|', 12) + pad('max rel(scale)', 16) + 'max rel(|E_ref|)');
const order = [...stats.keys()].filter(k => k !== 'ALL').sort();
for (const k of order) {
  const b = stats.get(k);
  console.log(pad(k, 26) + pad(b.n, 8) + pad(num(b.maxAbs), 12) + pad(num(b.maxRelScale), 16) + num(b.maxRelRef));
}
const all = stats.get('ALL');
console.log(pad('ALL', 26) + pad(all.n, 8) + pad(num(all.maxAbs), 12) + pad(num(all.maxRelScale), 16) + num(all.maxRelRef));

// well conditioned = everything except the deliberately degenerate family and the
// samples where one of the mirror's floors took over, which are deviations by design
let wc = { n: 0, maxAbs: 0, maxRelScale: 0, maxRelRef: 0 };
for (const k of order) {
  if (k.startsWith('e ') || k.startsWith('z ')) continue;
  const b = stats.get(k);
  wc.n += b.n;
  wc.maxAbs = Math.max(wc.maxAbs, b.maxAbs);
  wc.maxRelScale = Math.max(wc.maxRelScale, b.maxRelScale);
  wc.maxRelRef = Math.max(wc.maxRelRef, b.maxRelRef);
}
let dg = { n: 0, maxAbs: 0, maxRelScale: 0, maxRelRef: 0 };
for (const k of order) {
  if (!k.startsWith('e ')) continue;
  const b = stats.get(k);
  dg.n += b.n;
  dg.maxAbs = Math.max(dg.maxAbs, b.maxAbs);
  dg.maxRelScale = Math.max(dg.maxRelScale, b.maxRelScale);
  dg.maxRelRef = Math.max(dg.maxRelRef, b.maxRelRef);
}
console.log('\nwell conditioned (a-d): n ' + wc.n + '  max |dE| ' + num(wc.maxAbs) +
  '  max rel ' + num(wc.maxRelScale) + '   = ' + (wc.maxRelScale / EPS32).toExponential(2) + ' x float32 eps');
console.log('near-degenerate  (e) : n ' + dg.n + '  max |dE| ' + num(dg.maxAbs) +
  '  max rel ' + num(dg.maxRelScale) + '   = ' + (dg.maxRelScale / EPS32).toExponential(2) + ' x float32 eps');

console.log('\nleaves of the EMF selection, as taken by the mirror:');
const LEAFNAME = ['', 'ELL  (SB>0,SL>0)', 'ERL  (SB>0,SR<0)', 'x-fan bottom (SB>0)', 'ELR  (ST<0,SL>0)',
  'ERR  (ST<0,SR<0)', 'x-fan top    (ST<0)', 'y-fan left   (SL>0)', 'y-fan right  (SR<0)', 'full 2-D star'];
for (let b = 1; b <= 9; b++) console.log('  ' + b + '  ' + pad(LEAFNAME[b], 22) + branchCount[b] + (branchCount[b] >= 200 ? '' : '   <-- UNDER 200'));

console.log('\nref non-finite (unguarded divide in the Fortran): ' + nRefBad);
console.log('mirror non-finite (must be zero):                 ' + nMirrorBad);
console.log('leaf disagreements, no floor active (must be zero):' + nBranchMismatch);
console.log('leaf disagreements where a floor bound:           ' + nBranchMismatchFloored);

// How sharp the reference's own selection is at a boundary. This is not a check on the
// port; it is the reason the port cannot promise more than bit-identical branch
// selection. Walk one quadrant's v through ST = 0 and ask the ORACLE alone how much E
// moves: the two formulas either side average over different sub-fans and are under no
// obligation to meet, so E jumps, and a corner within one ulp of the boundary is
// decided by the last bit of the fast speed.
let maxJump = 0, maxJumpRel = 0;
for (let i = 0; i < 400; i++) {
  const adia = rng() < 0.5, cs2 = 1.0, gam = 5 / 3;
  const params = { adia: adia ? 1 : 0, gam, cs2, smallc: 1e-3 };
  const s = [];
  for (let k = 0; k < 4; k++) {
    const r = lg(0.3, 3);
    s.push({ r, p: adia ? lg(0.3, 3) : r * cs2, ...(([A, B]) => ({ A, B }))(field('generic', k)) });
  }
  const cx = Math.max(...s.map(q => findSpeedFastRef(q.r, q.p, q.A, q.B, params)));
  const cy = Math.max(...s.map(q => findSpeedFastRef(q.r, q.p, q.B, q.A, params)));
  for (let k = 0; k < 4; k++) { s[k].u = cx * 0.9 * (2 * rng() - 1); s[k].v = -cy - lg(0.01, 2); }
  const top = Math.floor(rng() * 4);
  const at = v => { const t = s.map(q => ({ ...q })); t[top].v = v; return riemann2dRef(t[0], t[1], t[2], t[3], params); };
  const a = at(-cy), b = at(-cy - 1e-12);
  if (Number.isFinite(a.E) && Number.isFinite(b.E) && a.branch !== b.branch) {
    const j = Math.abs(a.E - b.E);
    if (j > maxJump) maxJump = j;
    const rel = j / Math.max(Math.abs(a.E), Math.abs(b.E));
    if (rel > maxJumpRel) maxJumpRel = rel;
  }
}
console.log('\nthe oracle\'s own discontinuity across ST = 0, over 400 states:');
console.log('  max |E(ST=0) - E(ST=-1e-12)| ' + num(maxJump) + ',  as a fraction of |E| ' + num(maxJumpRel));

console.log('\nhow close the campaign drove each guarded quantity to its floor');
console.log('(1.0 means the guard took over; >1 means it never did):');
for (const k in near) console.log('  ' + pad(k, 28) + (Number.isFinite(near[k]) ? near[k].toExponential(3) : 'n/a'));
console.log('  samples in which some guard bound: ' + nGuardBound);
console.log('\nthe same contact distances measured against the fast speed, not the floor:');
for (const k in scaleFree) console.log('  ' + pad(k, 28) + (Number.isFinite(scaleFree[k]) ? scaleFree[k].toExponential(3) : 'n/a'));

console.log('\nh hostile states, mirror only (no oracle exists -- the Fortran faults):');
console.log('  states                                          ' + nHostile);
console.log('  ref non-finite, i.e. guards that earned their keep: ' + nHostileRefBad);
console.log('  mirror non-finite or out of float32 range:          ' + nHostileMirrorBad);

console.log('\ni transposition antisymmetry of the mirror, over ' + nSym + ' corners:');
console.log('  max |E + E_transposed| ' + num(symAbs) + ',  as a fraction of |E| ' + num(symRel));

console.log('\ntranscription check, mirror against riemann2d-hlld.glsl:');
for (const t of transcript) console.log('  ' + pad(t.name, 12) + (t.same ? 'token-identical' : 'DIFFERS at char ' + t.at + '\n    glsl: ' + t.a.slice(Math.max(0, t.at - 40), t.at + 60) + '\n    js  : ' + t.b.slice(Math.max(0, t.at - 40), t.at + 60)));

const checks = [
  ['mirror is token-identical to the GLSL', transcript.every(t => t.same)],
  // Two of the floors cannot be reached by a state with positive density and pressure,
  // which is a property to measure rather than a gap to apologise for: their exact
  // values are unobservable, and they exist only for the hostile set.
  ['contact floor unreachable: |SL-ustar| > 1e-3 cx', Math.min(scaleFree['SL-ustar / cx'], scaleFree['SR-ustar / cx'], scaleFree['SB-vstar / cy'], scaleFree['ST-vstar / cy']) > 1e-3],
  ['rc floor unreachable: sum rc > 1e6 x 1e-20', Math.min(near['sum rc x   (floor 1e-20)'], near['sum rc y   (floor 1e-20)']) > 1e6],
  ['sample count >= 20000', nTotal >= 20000],
  ['transposition antisymmetry holds to 1e-10', symRel < 1e-10],
  ['mirror finite and in float32 range on hostile states', nHostileMirrorBad === 0],
  ['mirror never non-finite', nMirrorBad === 0],
  ['leaf selection identical to the Fortran', nBranchMismatch === 0],
  ['every leaf exercised >= 200 times', branchCount.slice(1).every(c => c >= 200)],
  ['well-conditioned agreement < 1 float32 eps', wc.maxRelScale < EPS32],
  ['near-degenerate agreement < 1 float32 eps', dg.maxRelScale < EPS32]
];
console.log('');
let ok = true;
for (const [name, pass] of checks) { console.log((pass ? '  ok   ' : '  FAIL ') + name); ok = ok && pass; }
if (!ok) {
  if (worstMirrorBad) console.log('\nfirst non-finite mirror: ' + JSON.stringify(worstMirrorBad));
  if (worstBranchMismatch) console.log('\nfirst leaf disagreement: ' + JSON.stringify(worstBranchMismatch));
  for (const k of order) {
    const b = stats.get(k);
    if (b.maxRelScale >= EPS32 && b.worst) console.log('\nworst in ' + k + ': ' + JSON.stringify(b.worst));
  }
}
console.log('\n' + (ok ? 'PASS' : 'FAIL'));
process.exit(ok ? 0 : 1);

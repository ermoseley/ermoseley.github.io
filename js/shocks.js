/* =========================================================================
   shocks.js — the page around the compressible solver.

   Wires the controls, draws the density PDF against its log-normal prediction,
   and prints the diagnostics. Nothing here is physics; the physics is all in
   euler.js.
   ========================================================================= */

(function () {
  'use strict';

  const $ = (s) => document.querySelector(s);
  const canvas = $('#euler');
  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;

  const sim = window.Euler ? window.Euler.create(canvas, { n: 288, mach: 4.0, zeta: 0.35 }) : null;

  const yearEl = $('#year');
  if (yearEl) yearEl.textContent = String(new Date().getFullYear());

  if (!sim) {
    $('#fallback').hidden = false;
    $('#hint').hidden = true;
    return;
  }
  if (reduced) sim.setSpeed(1);

  // ------------------------------------------------------------- controls

  function slider(id, out, fmt, apply) {
    const el = $(id), lab = $(out);
    function sync() {
      const v = parseFloat(el.value);
      lab.textContent = fmt(v);
      apply(v);
    }
    el.addEventListener('input', sync);
    sync();
  }

  slider('#mach', '#mach-v', (v) => v.toFixed(1), (v) => sim.setMach(v));
  slider('#zeta', '#zeta-v', (v) => v.toFixed(2), (v) => sim.setZeta(v));

  // Resolution reallocates every buffer, so only act when the drag settles.
  let resT = null;
  const resEl = $('#res');
  resEl.addEventListener('input', function () {
    $('#res-v').textContent = resEl.value;
    clearTimeout(resT);
    resT = setTimeout(function () { sim.setResolution(parseInt(resEl.value, 10)); }, 220);
  });

  const MODES = [
    'density + shocks',
    'convergence only',
    'Mach number',
    'divergence, signed'
  ];
  let mode = 0;
  const modeBtn = $('#mode');
  modeBtn.addEventListener('click', function () {
    mode = (mode + 1) % MODES.length;
    sim.setMode(mode);
    modeBtn.textContent = 'view: ' + MODES[mode];
  });

  // HLL is the default; LLF is one wave speed instead of two and correspondingly
  // more diffusive. Backend-selectable and worth exposing, since which one you
  // pick is exactly the kind of thing this page is about.
  let llf = false;
  const rieBtn = $('#riemann');
  if (rieBtn) rieBtn.addEventListener('click', function () {
    llf = !llf;
    sim.setRiemann(llf ? 'llf' : 'hll');
    rieBtn.textContent = 'flux: ' + (llf ? 'LLF' : 'HLL');
  });

  const pauseBtn = $('#pause');
  pauseBtn.addEventListener('click', function () {
    pauseBtn.textContent = sim.toggle() ? 'pause' : 'play';
  });

  $('#reset').addEventListener('click', function () { sim.reset(); });

  document.addEventListener('keydown', function (e) {
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    const t = e.target;
    if (t && (t.isContentEditable || /^(INPUT|TEXTAREA|SELECT|BUTTON)$/.test(t.tagName))) return;
    if (e.key === ' ') { e.preventDefault(); pauseBtn.click(); }
    else if (e.key === 'v') modeBtn.click();
    else if (e.key === 'r') sim.reset();
  });

  // ------------------------------------------------------------- pointer

  (function pointer() {
    let px = null, py = null, drag = false, moved = false, downX = 0, downY = 0;
    let hinted = false;

    function uv(e) {
      const r = canvas.getBoundingClientRect();
      return {
        x: (e.clientX - r.left) / Math.max(r.width, 1),
        y: 1 - (e.clientY - r.top) / Math.max(r.height, 1)
      };
    }

    function hideHint() {
      if (hinted) return;
      hinted = true;
      const h = $('#hint');
      if (h) h.classList.add('gone');
    }

    canvas.addEventListener('pointerdown', function (e) {
      const p = uv(e);
      px = downX = p.x; py = downY = p.y;
      drag = true; moved = false;
      hideHint();
      // Capture is a nicety -- it keeps the stroke alive past the canvas edge --
      // and it is allowed to fail. It used to sit ahead of the rest of this
      // handler, so when it threw (which it does for a synthetic pointer, and
      // can for a real one) it took the whole interaction down with it.
      try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* not fatal */ }
    });

    canvas.addEventListener('pointermove', function (e) {
      if (!drag) return;
      const p = uv(e);
      if (Math.hypot(p.x - downX, p.y - downY) > 0.012) moved = true;
      if (px !== null) {
        const dx = p.x - px, dy = p.y - py;
        // Scaled against the sound speed: a stroke should push the gas at a
        // Mach number you can actually see. The previous factor left the stir
        // an order of magnitude below the detonation it always followed, which
        // is why it read as doing nothing at all.
        // Clamped: a fast mouse can report a tenth of the box in one event, and
        // unclamped that is a Mach 14 kick injected between two CFL measurements.
        const sp = Math.hypot(dx, dy) * 140, cap = 5;
        const k = sp > cap ? cap / sp : 1;
        if (sp > 0.07) sim.push(p.x, p.y, dx * 140 * k, dy * 140 * k);
      }
      px = p.x; py = p.y;
    });

    function end(e) {
      // A press that never moved is a click, and a click detonates. Separating
      // the two means a stir is a stir: it no longer begins with a blast wave
      // that swamps the thing you were trying to see.
      if (drag && !moved) sim.blast(downX, downY, 7.0);
      drag = false; moved = false; px = py = null;
      try {
        if (e && e.pointerId != null && canvas.hasPointerCapture(e.pointerId)) {
          canvas.releasePointerCapture(e.pointerId);
        }
      } catch (err) { /* not fatal */ }
    }
    canvas.addEventListener('pointerup', end);
    canvas.addEventListener('pointercancel', end);
  })();

  // ------------------------------------------------------------- the plot

  const pdfCanvas = $('#pdf');
  const ctx = pdfCanvas.getContext('2d');

  function sizePlot() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = pdfCanvas.clientWidth || 560;
    const h = pdfCanvas.clientHeight || 220;
    pdfCanvas.width = Math.round(w * dpr);
    pdfCanvas.height = Math.round(h * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { w, h };
  }
  let plotSize = sizePlot();
  window.addEventListener('resize', function () { plotSize = sizePlot(); });

  const INK = '#e8ebf4', DIM = 'rgba(232,235,244,0.34)', GRID = 'rgba(232,235,244,0.10)';
  const MEAS = 'rgba(120,158,255,0.85)', PRED = '#ff9f45';

  function drawPDF(s) {
    const { w, h } = plotSize;
    const padL = 40, padR = 10, padT = 12, padB = 26;
    const iw = w - padL - padR, ih = h - padT - padB;
    ctx.clearRect(0, 0, w, h);
    if (!s.pdf) return;

    const bins = s.pdf, NB = bins.length;
    const lo = s.pdfLo - s.meanLn, hi = s.pdfHi - s.meanLn;   // centre on the mean

    // the prediction: log-normal in rho, so Gaussian in ln rho, with the mean
    // fixed by mass conservation at -sigma^2/2 relative to ln(mean rho)
    const sp = Math.max(s.sigmaPred, 1e-3);   // from b(mixture) and the running-mean Mach
    const pred = (x) => Math.exp(-0.5 * Math.pow((x + 0.5 * sp * sp) / sp, 2)) / (sp * Math.sqrt(2 * Math.PI));

    let peak = 0;
    for (let i = 0; i < NB; i++) peak = Math.max(peak, bins[i]);
    for (let i = 0; i <= 60; i++) peak = Math.max(peak, pred(lo + (hi - lo) * i / 60));
    peak = Math.max(peak, 1e-6) * 1.12;

    const X = (v) => padL + (v - lo) / Math.max(hi - lo, 1e-6) * iw;
    const Y = (v) => padT + ih - (v / peak) * ih;

    // frame + zero line
    ctx.strokeStyle = GRID;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padL, padT); ctx.lineTo(padL, padT + ih); ctx.lineTo(padL + iw, padT + ih);
    ctx.stroke();

    ctx.strokeStyle = GRID;
    ctx.setLineDash([2, 3]);
    ctx.beginPath(); ctx.moveTo(X(0), padT); ctx.lineTo(X(0), padT + ih); ctx.stroke();
    ctx.setLineDash([]);

    // measured histogram
    ctx.fillStyle = MEAS;
    const bw = iw / NB;
    for (let i = 0; i < NB; i++) {
      const x = padL + i * bw;
      const y = Y(bins[i]);
      ctx.fillRect(x, y, Math.max(bw - 1, 1), padT + ih - y);
    }

    // prediction
    ctx.strokeStyle = PRED;
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i <= 120; i++) {
      const x = lo + (hi - lo) * i / 120;
      const y = pred(x);
      if (i === 0) ctx.moveTo(X(x), Y(y)); else ctx.lineTo(X(x), Y(y));
    }
    ctx.stroke();

    // axis labels
    ctx.fillStyle = DIM;
    ctx.font = '10px "JetBrains Mono", monospace';
    ctx.textAlign = 'center';
    for (const v of [-4, -2, 0, 2, 4]) {
      if (v < lo || v > hi) continue;
      ctx.fillText(String(v), X(v), padT + ih + 15);
    }
    ctx.textAlign = 'right';
    ctx.fillText('p', padL - 6, padT + 8);
    ctx.textAlign = 'left';
    ctx.fillStyle = INK;
  }

  // ------------------------------------------------------------- readouts

  const R = {
    grid: $('#r-grid'), mach: $('#r-mach'), machmax: $('#r-machmax'),
    sigma: $('#r-sigma'), sigmap: $('#r-sigmap'), bfit: $('#r-bfit'),
    shock: $('#r-shock'), mass: $('#r-mass'), flux: $('#r-flux'),
    resets: $('#r-resets'),
    dt: $('#r-dt'), steps: $('#r-steps'), fps: $('#r-fps')
  };

  // Written defensively: the readout table and this file drift apart easily, and
  // a missing row used to throw inside the stats callback on every tick, which
  // took every other number down with it.
  function put(el, text) { if (el) el.textContent = text; }

  sim.onStats(function (s) {
    put(R.grid, s.grid + '  (' + (s.gridW * s.gridH).toLocaleString() + ' cells)');
    put(R.mach, s.machRms.toFixed(2));
    put(R.machmax, s.machMax.toFixed(2));
    put(R.sigma, s.sigma.toFixed(3) + '   (mean ' + (s.sigmaBar || s.sigma).toFixed(3) + ')');
    put(R.sigmap, s.sigmaPred.toFixed(3) + '   (b = ' + s.bPred.toFixed(2) + ')');
    put(R.bfit, isFinite(s.bFit) ? s.bFit.toFixed(2) : '—');
    put(R.shock, (100 * s.compress).toFixed(1) + ' %');
    put(R.mass, s.massErr < 1e-12 ? '0 (exact)' : s.massErr.toExponential(1) + '  relative');
    put(R.flux, (s.llf ? 'LLF (Rusanov)' : 'HLL') + ', piecewise constant');
    put(R.resets, s.resets ? String(s.resets) + ' — the state went non-finite' : 'none');
    put(R.dt, s.dt.toExponential(2) + ' · ' + s.cfl.toFixed(2));
    put(R.steps, s.steps.toLocaleString());
    put(R.fps, s.fps.toFixed(0));
    drawPDF(s);
  });
})();

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
    let px = null, py = null, drag = false;
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
      px = p.x; py = p.y; drag = true;
      canvas.setPointerCapture(e.pointerId);
      // an overdensity plus the outward momentum to match: in an isothermal gas
      // piling up density *is* piling up pressure, so this relaxes into a real
      // outward-running shock rather than a painted-on ring
      sim.blast(p.x, p.y, 7.0, 5.0);
      hideHint();
    });

    canvas.addEventListener('pointermove', function (e) {
      if (!drag) return;
      const p = uv(e);
      if (px !== null) {
        const dx = p.x - px, dy = p.y - py;
        if (Math.hypot(dx, dy) > 0.002) sim.push(p.x, p.y, dx * 26, dy * 26);
      }
      px = p.x; py = p.y;
    });

    function end(e) {
      drag = false; px = py = null;
      if (e && e.pointerId != null && canvas.hasPointerCapture(e.pointerId)) {
        canvas.releasePointerCapture(e.pointerId);
      }
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
    shock: $('#r-shock'), mass: $('#r-mass'),
    dt: $('#r-dt'), steps: $('#r-steps'), fps: $('#r-fps')
  };

  sim.onStats(function (s) {
    R.grid.textContent = s.grid + '  (' + (s.gridW * s.gridH).toLocaleString() + ' cells)';
    R.mach.textContent = s.machRms.toFixed(2);
    R.machmax.textContent = s.machMax.toFixed(2);
    R.sigma.textContent = s.sigma.toFixed(3) + '   (mean ' + (s.sigmaBar || s.sigma).toFixed(3) + ')';
    R.sigmap.textContent = s.sigmaPred.toFixed(3) + '   (b = ' + s.bPred.toFixed(2) + ')';
    R.bfit.textContent = isFinite(s.bFit) ? s.bFit.toFixed(2) : '—';
    R.shock.textContent = (100 * s.compress).toFixed(1) + ' %';
    R.mass.textContent = s.massErr < 1e-12 ? '0 (exact)' : s.massErr.toExponential(1) + '  relative';
    R.dt.textContent = s.dt.toExponential(2) + ' · ' + s.cfl.toFixed(2);
    R.steps.textContent = s.steps.toLocaleString();
    R.fps.textContent = s.fps.toFixed(0);
    drawPDF(s);
  });
})();

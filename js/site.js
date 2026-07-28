/* =========================================================================
   site.js — navigation, reveals, and the coupling between the page and the
   simulation running behind it.
   ========================================================================= */

(function () {
  'use strict';

  const doc = document;
  const $ = (s, r) => (r || doc).querySelector(s);
  const $$ = (s, r) => Array.prototype.slice.call((r || doc).querySelectorAll(s));
  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;

  // ------------------------------------------------------------ simulation

  const canvas = $('#field');
  const field = window.Field ? window.Field.create(canvas) : null;

  // --------------------------------------------------------------- sections

  const sections = $$('[data-nav]');

  // ------------------------------------------------------------------- boot

  (function boot() {
    const el = $('#boot');
    if (!el) return;
    if (!field || sessionStorage.getItem('wiem-booted')) {
      el.classList.add('done');
      return;
    }
    const s = field.stats();
    $('#boot-grid').textContent = s.gridW ? s.gridW + ' × ' + s.gridH : 'canvas fallback';
    $('#boot-part').textContent = s.nPart.toLocaleString();
    $('#boot-jac').textContent = s.jacobi ? 'jacobi × ' + s.jacobi : 'curl-noise';
    setTimeout(function () {
      el.classList.add('done');
      sessionStorage.setItem('wiem-booted', '1');
    }, 1250);
  })();

  // -------------------------------------------------------------- nav build

  (function nav() {
    const rail = $('#rail');
    const sheet = $('#sheet-list');
    sections.forEach(function (sec) {
      const n = sec.dataset.nav, label = sec.dataset.label;
      if (rail) {
        const a = doc.createElement('a');
        a.href = '#' + sec.id;
        a.dataset.for = sec.id;
        a.innerHTML = '<u>' + n + '</u><i></i><span>' + label + '</span>';
        rail.appendChild(a);
      }
      if (sheet) {
        const li = doc.createElement('li');
        li.innerHTML = '<a href="#' + sec.id + '"><u>' + n + '</u>' + label + '</a>';
        sheet.appendChild(li);
      }
    });

    const btn = $('#menu-btn'), panel = $('#sheet');
    if (!btn || !panel) return;
    function close() {
      panel.classList.remove('open');
      btn.setAttribute('aria-expanded', 'false');
      btn.textContent = 'Index';
      doc.body.classList.remove('is-locked');
    }
    btn.addEventListener('click', function () {
      const open = panel.classList.toggle('open');
      btn.setAttribute('aria-expanded', String(open));
      btn.textContent = open ? 'Close' : 'Index';
      doc.body.classList.toggle('is-locked', open);
    });
    panel.addEventListener('click', function (e) { if (e.target.closest('a')) close(); });
    doc.addEventListener('keydown', function (e) { if (e.key === 'Escape') close(); });
  })();

  // ---------------------------------------------------------------- reveals

  (function reveals() {
    const items = $$('.rv');
    if (reduced || !('IntersectionObserver' in window)) {
      items.forEach(function (n) { n.classList.add('in'); });
      return;
    }
    const io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) { e.target.classList.add('in'); io.unobserve(e.target); }
      });
    }, { rootMargin: '0px 0px -10% 0px', threshold: 0.04 });

    // Anything already on screen at load must reveal regardless of the
    // bottom rootMargin, or the fold of the hero stays invisible.
    const h = window.innerHeight;
    items.forEach(function (n) {
      const r = n.getBoundingClientRect();
      if (r.top < h && r.bottom > 0) n.classList.add('in');
      else io.observe(n);
    });
  })();

  // -------------------------------------------- active section + accent + sim

  const railLinks = {};
  $$('#rail a').forEach(function (a) { railLinks[a.dataset.for] = a; });

  let current = null;

  // sections that get an impulse when you arrive
  const PUNCH = { phrike: 'blast', skate: 'blast', iron: 'shear' };

  function activate(sec) {
    if (!sec || sec === current) return;
    current = sec;

    doc.documentElement.style.setProperty('--accent', sec.dataset.accent || '#b9c6ff');

    Object.keys(railLinks).forEach(function (k) { railLinks[k].classList.remove('on'); });
    if (railLinks[sec.id]) railLinks[sec.id].classList.add('on');

    if (field) {
      field.setPreset(sec.dataset.preset || 'hero');
      const p = PUNCH[sec.dataset.preset];
      if (p === 'blast' && !reduced) {
        field.blast(0.5, 0.5, sec.dataset.preset === 'phrike' ? 430 : 300, 0.012);
      } else if (p === 'shear' && !reduced) {
        field.shear(120);
      }
    }
  }

  const prog = $('#prog');
  let lastY = window.scrollY;
  let queued = false;

  function onScrollFrame() {
    queued = false;
    const y = window.scrollY;
    const dy = y - lastY;
    lastY = y;

    // scrolling shears the fluid: a tanh layer across mid-screen, which is a
    // genuine Kelvin-Helmholtz setup, so it rolls up as you read.
    if (field && Math.abs(dy) > 0.5 && !reduced) {
      field.shear(Math.max(-280, Math.min(280, dy * 3.0)));
    }

    const max = doc.documentElement.scrollHeight - window.innerHeight;
    if (prog) prog.style.width = (max > 0 ? (y / max) * 100 : 0) + '%';

    // the section nearest the middle of the viewport wins
    const mid = y + window.innerHeight * 0.42;
    let best = null, bestD = Infinity;
    for (let i = 0; i < sections.length; i++) {
      const r = sections[i].getBoundingClientRect();
      const top = r.top + y, c = top + r.height / 2;
      const d = Math.abs(c - mid);
      if (d < bestD) { bestD = d; best = sections[i]; }
    }
    activate(best);
  }

  window.addEventListener('scroll', function () {
    if (!queued) { queued = true; requestAnimationFrame(onScrollFrame); }
  }, { passive: true });

  // ------------------------------------------------------------- pointer

  (function pointer() {
    if (!field) return;
    let px = null, py = null;
    let curX = 0, curY = 0, accX = 0, accY = 0, touch = false, queued = false;

    // Each splat is two full-grid passes in the solver, and a fast mouse can
    // deliver dozens of pointermove events between frames. Accumulate them and
    // emit at most one splat per animation frame; visually it is the same
    // stroke, and it removes the worst cost spike on the page.
    function flush() {
      queued = false;
      const w = window.innerWidth, h = window.innerHeight;
      const sp = Math.hypot(accX, accY);
      if (sp > 0.0004) {
        field.splat(curX / w, 1 - curY / h, accX * 2600, accY * 2600, touch ? 0.012 : 0.0085);
      }
      accX = accY = 0;
    }

    function push(x, y, isTouch) {
      const w = window.innerWidth, h = window.innerHeight;
      if (px !== null) {
        accX += (x - px) / w;
        accY += -(y - py) / h;
      }
      px = x; py = y; curX = x; curY = y; touch = isTouch;
      if (!queued) { queued = true; requestAnimationFrame(flush); }
    }

    window.addEventListener('pointermove', function (e) {
      if (reduced) return;
      push(e.clientX, e.clientY, e.pointerType === 'touch');
    }, { passive: true });

    window.addEventListener('pointerdown', function (e) {
      if (reduced) return;
      px = e.clientX; py = e.clientY;
      field.blast(e.clientX / window.innerWidth, 1 - e.clientY / window.innerHeight, 300, 0.006);
    }, { passive: true });

    window.addEventListener('pointerleave', function () { px = py = null; });
  })();

  // ----------------------------------------------------------------- HUD

  (function hud() {
    if (!field) return;
    const gridEl = $('#hud-grid'), partEl = $('#hud-part'), tauEl = $('#hud-tau');
    const clockEl = $('#hud-clock');
    const rmsEl = $('#hud-rms'), maxEl = $('#hud-max'), divEl = $('#hud-div'), fpsEl = $('#hud-fps');

    function sci(v) {
      if (!isFinite(v)) return '—';
      if (v === 0) return '0';
      const e = Math.floor(Math.log10(Math.abs(v)));
      const m = v / Math.pow(10, e);
      return m.toFixed(1) + 'e' + (e < 0 ? '' : '+') + e;
    }

    function tick() {
      const s = field.stats();
      if (gridEl) gridEl.textContent = s.gridW ? s.gridW + '×' + s.gridH : 'canvas2d';
      if (partEl) partEl.textContent = (s.nPart / 1000).toFixed(0) + 'k';
      if (clockEl) clockEl.textContent = '1:' + Math.round(1 / (s.timeScale || 1));
      if (tauEl) tauEl.textContent = s.tau.toFixed(3) + ' s';
      if (rmsEl) rmsEl.textContent = s.rms.toFixed(1);
      if (maxEl) maxEl.textContent = s.max.toFixed(1);
      if (divEl) divEl.textContent = sci(s.div);
      if (fpsEl) fpsEl.textContent = s.fps.toFixed(0);
      setTimeout(tick, 420);
    }
    tick();

    // colophon reports the same measured configuration
    const s = field.stats();
    const cg = $('#colo-grid'), cp = $('#colo-part'), cj = $('#colo-jac'), fb = $('#foot-backend');
    if (cg) cg.textContent = s.gridW ? s.gridW + ' × ' + s.gridH : 'fallback';
    if (cp) cp.textContent = s.nPart.toLocaleString();
    if (cj && s.jacobi) cj.textContent = String(s.jacobi);
    if (fb) fb.textContent = field.kind === 'webgl2' ? 'WebGL2' : 'Canvas 2D fallback';
  })();

  // ------------------------------------------------- drop-in image slots

  (function slots() {
    function card(img) {
      const wrap = img.parentNode;
      if (!wrap || wrap.querySelector('.slot')) return;
      const file = img.dataset.file || img.getAttribute('src') || '';
      const ratio = img.dataset.ratio || '';
      const hint = img.dataset.hint || '';
      const el = doc.createElement('div');
      el.className = 'slot';
      el.innerHTML =
        '<b>drop an image here</b>' +
        '<code>' + file + '</code>' +
        (ratio ? '<i>' + ratio + '</i>' : '') +
        (hint ? '<i>' + hint + '</i>' : '');
      img.style.display = 'none';
      wrap.appendChild(el);
    }

    $$('img[data-slot]').forEach(function (img) {
      if (img.getAttribute('data-failed') || (img.complete && img.naturalWidth === 0)) {
        card(img);
      } else if (!img.complete) {
        img.addEventListener('error', function () { card(img); });
      } else {
        // loaded fine — leave the real photograph alone
      }
    });

    // any non-slot image that fails should not leave an empty hole
    $$('img:not([data-slot])').forEach(function (img) {
      img.addEventListener('error', function () { img.style.opacity = '0'; });
    });
  })();

  // ------------------------------------------------------------ video hover

  (function video() {
    $$('video').forEach(function (v) {
      const frame = v.closest('.frame') || v;
      let playing = false;
      function play() { if (reduced) return; v.play().then(function () { playing = true; }).catch(function () {}); }
      function stop() { v.pause(); playing = false; }
      frame.addEventListener('mouseenter', play);
      frame.addEventListener('mouseleave', stop);
      frame.addEventListener('click', function () { playing ? stop() : play(); });
    });
  })();

  // ------------------------------------------------------------------- misc

  const yr = $('#year');
  if (yr) yr.textContent = String(new Date().getFullYear());

  if (window.Blog) window.Blog.render($('#post-list'), 'blog/', 3);

  onScrollFrame();
})();

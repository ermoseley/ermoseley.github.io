/* =========================================================================
   site.js — the tab deck.

   Thirteen chapters, exactly one of which is in the document at a time. This
   file is the router between them: it owns the history entries, replays the
   reveal cascade whenever a chapter opens, and turns each page-turn into a
   physical event in the simulation running behind the type.
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

  // ----------------------------------------------------------------- deck

  const panels = $$('#shell [data-nav]');
  const byId = {};
  panels.forEach(function (el, i) { byId[el.id] = { id: el.id, el: el, i: i, tabs: [] }; });
  const HOME = panels.length ? panels[0].id : null;

  const rail = $('#rail');
  const strip = $('#strip');
  const prog = $('#prog');

  // ------------------------------------------------------------------- boot

  function bootThen(fn) {
    const el = $('#boot');
    if (!el) { fn(); return; }
    if (!field || sessionStorage.getItem('wiem-booted')) {
      el.classList.add('done');
      fn();
      return;
    }
    const s = field.stats();
    $('#boot-grid').textContent = s.gridW ? s.gridW + ' × ' + s.gridH : 'canvas fallback';
    $('#boot-part').textContent = s.nPart.toLocaleString();
    $('#boot-jac').textContent = s.solver || 'curl-noise';
    setTimeout(function () {
      el.classList.add('done');
      sessionStorage.setItem('wiem-booted', '1');
      // the first chapter opens as the readout clears, so its cascade is
      // actually seen rather than played behind an opaque panel
      fn();
    }, 1150);
  }

  // -------------------------------------------------------------- nav build

  (function build() {
    const sheet = $('#sheet-list');

    function tab(host, kind, rec, n, label) {
      const a = doc.createElement('a');
      a.href = '#' + rec.id;
      a.id = 'tab-' + kind + '-' + rec.id;
      a.setAttribute('role', 'tab');
      a.setAttribute('aria-controls', rec.id);
      a.setAttribute('aria-selected', 'false');
      a.tabIndex = -1;
      a.innerHTML = kind === 'rail'
        ? '<u>' + n + '</u><i></i><span>' + label + '</span>'
        : '<u>' + n + '</u><span>' + label + '</span>';
      host.appendChild(a);
      rec.tabs.push(a);
      return a;
    }

    panels.forEach(function (el) {
      const rec = byId[el.id];
      const n = el.dataset.nav, label = el.dataset.label;
      if (rail) el.setAttribute('aria-labelledby', tab(rail, 'rail', rec, n, label).id);
      if (strip) tab(strip, 'strip', rec, n, label);
      if (sheet) {
        const li = doc.createElement('li');
        li.innerHTML = '<a href="#' + rec.id + '"><u>' + n + '</u>' + label + '</a>';
        sheet.appendChild(li);
      }
    });
  })();

  // the full-screen index, still available on narrow screens
  const sheetPanel = $('#sheet');
  const menuBtn = $('#menu-btn');

  function closeSheet() {
    if (!sheetPanel || !menuBtn) return;
    sheetPanel.classList.remove('open');
    menuBtn.setAttribute('aria-expanded', 'false');
    menuBtn.textContent = 'Index';
    doc.body.classList.remove('is-locked');
  }

  if (menuBtn && sheetPanel) {
    menuBtn.addEventListener('click', function () {
      const open = sheetPanel.classList.toggle('open');
      menuBtn.setAttribute('aria-expanded', String(open));
      menuBtn.textContent = open ? 'Close' : 'Index';
      doc.body.classList.toggle('is-locked', open);
    });
    doc.addEventListener('keydown', function (e) { if (e.key === 'Escape') closeSheet(); });
  }

  // ---------------------------------------------------------------- reveals

  const io = (!reduced && 'IntersectionObserver' in window)
    ? new IntersectionObserver(function (entries) {
        entries.forEach(function (e) {
          if (e.isIntersecting) { e.target.classList.add('in'); io.unobserve(e.target); }
        });
      }, { rootMargin: '0px 0px -8% 0px', threshold: 0.04 })
    : null;

  /* Called while a chapter is hidden, so nothing animates backwards. */
  function resetReveals(root) {
    if (!io) return;
    $$('.rv', root).forEach(function (n) {
      io.unobserve(n);
      n.classList.remove('in');
      n.style.removeProperty('--ds');
    });
  }

  /* Called once a chapter is on screen and scrolled to the top. Everything in
     the opening fold cascades in; everything below it waits for the observer,
     exactly as it did when the site was one long scroll. */
  function reveal(root) {
    const items = $$('.rv', root);
    if (!io) { items.forEach(function (n) { n.classList.add('in'); }); return; }

    // resolve the just-displayed subtree's style before we add the class,
    // or the browser coalesces display:none -> visible -> .in and skips the
    // transition entirely
    void root.offsetHeight;

    const h = window.innerHeight;
    const rects = items.map(function (n) { return n.getBoundingClientRect(); });
    let k = 0;
    items.forEach(function (n, idx) {
      const r = rects[idx];
      if (r.top < h * 0.96 && r.bottom > 0) {
        n.style.setProperty('--ds', Math.min(k * 42, 480) + 'ms');
        k++;
        n.classList.add('in');
      } else {
        io.observe(n);
      }
    });
  }

  // ------------------------------------------------------------ first visit

  /* Images and the blog feed only exist once their chapter has been opened —
     the point of the tab layout is that the other twelve cost nothing. */
  function firstVisit(rec) {
    if (rec.seen) return;
    rec.seen = true;
    wireSlots(rec.el);
    if (rec.id === 'blog' && window.Blog) window.Blog.render($('#post-list'), 'blog/', 3);
  }

  // -------------------------------------------------------------- the field

  // chapters that get an impulse of their own when you arrive
  const PUNCH = { phrike: 'blast', skate: 'blast', iron: 'shear' };

  function drive(rec, dir, origin) {
    if (!field) return;
    field.setPreset(rec.el.dataset.preset || 'hero');
    if (reduced) return;

    // turning to a new chapter is a physical event: a shear layer across the
    // screen in the direction of travel, plus an impulse where you pressed
    field.shear(dir * 190);
    if (origin) field.blast(origin.x, origin.y, 230, 0.0075);

    const p = PUNCH[rec.el.dataset.preset];
    if (p === 'blast') field.blast(0.5, 0.5, rec.el.dataset.preset === 'phrike' ? 430 : 300, 0.012);
    else if (p === 'shear') field.shear(120);
  }

  // --------------------------------------------------------------- routing

  const OUT_MS = reduced ? 0 : 190;
  let cur = null, busy = false, pending = null, lastY = 0;

  function markTabs(rec, on) {
    rec.tabs.forEach(function (a) {
      a.classList.toggle('on', on);
      a.setAttribute('aria-selected', on ? 'true' : 'false');
      a.tabIndex = on ? 0 : -1;
    });
  }

  function centreStrip(rec) {
    if (!strip || !strip.clientWidth) return;
    let a = null;
    for (let i = 0; i < rec.tabs.length; i++) {
      if (rec.tabs[i].parentNode === strip) { a = rec.tabs[i]; break; }
    }
    if (!a) return;
    const max = strip.scrollWidth - strip.clientWidth;
    const want = a.offsetLeft - strip.clientWidth / 2 + a.offsetWidth / 2;
    const left = Math.max(0, Math.min(max, want));
    if (strip.scrollTo) strip.scrollTo({ left: left, behavior: reduced ? 'auto' : 'smooth' });
    else strip.scrollLeft = left;
  }

  function enter(rec, dir, opt) {
    cur = rec;
    doc.documentElement.style.setProperty('--accent', rec.el.dataset.accent || '#b9c6ff');

    resetReveals(rec.el);
    rec.el.classList.add('is-on');
    markTabs(rec, true);
    centreStrip(rec);
    if (prog) prog.style.width = (((rec.i + 1) / panels.length) * 100).toFixed(2) + '%';

    window.scrollTo(0, 0);
    lastY = window.scrollY;          // do not let the jump register as a shear

    reveal(rec.el);
    firstVisit(rec);
    drive(rec, dir, opt.origin);

    if (opt.focus && rec.el.focus) rec.el.focus({ preventScroll: true });
  }

  function go(id, opt) {
    opt = opt || {};
    const rec = byId[id];
    if (!rec) return;
    if (busy) { pending = { id: id, opt: opt }; return; }
    if (rec === cur) { window.scrollTo(0, 0); return; }

    const prev = cur;
    if (!prev) { enter(rec, 1, opt); return; }

    const dir = rec.i > prev.i ? 1 : -1;
    busy = true;
    prev.el.classList.remove('is-on');
    prev.el.classList.add('is-off');
    markTabs(prev, false);

    setTimeout(function () {
      prev.el.classList.remove('is-off');
      resetReveals(prev.el);
      busy = false;
      enter(rec, dir, opt);
      if (pending) { const p = pending; pending = null; go(p.id, p.opt); }
    }, OUT_MS);
  }

  /* One history entry per chapter, so the back button walks the deck. */
  function nav(id, origin, focus) {
    if (!byId[id]) return;
    closeSheet();
    if (cur && byId[id] === cur) { window.scrollTo(0, 0); return; }
    try { history.pushState(null, '', '#' + id); }
    catch (err) { location.hash = id; }
    go(id, { origin: origin, focus: focus });
  }

  function fromUrl() {
    const id = (location.hash || '').slice(1);
    go(byId[id] ? id : HOME, {});
  }

  window.addEventListener('popstate', fromUrl);
  window.addEventListener('hashchange', fromUrl);

  // every in-page anchor is a chapter change, including the ones buried in
  // body copy (#papers from the dust chapter, #top from the wordmark)
  doc.addEventListener('click', function (e) {
    if (e.defaultPrevented || e.button !== 0) return;
    if (e.metaKey || e.ctrlKey || e.shiftKey || e.altKey) return;
    const a = e.target.closest ? e.target.closest('a[href]') : null;
    if (!a) return;
    const href = a.getAttribute('href') || '';
    if (href.charAt(0) !== '#' || href.length < 2) return;
    const id = href.slice(1);
    if (!byId[id]) return;
    e.preventDefault();
    const origin = (e.clientX || e.clientY)
      ? { x: e.clientX / window.innerWidth, y: 1 - e.clientY / window.innerHeight }
      : null;
    nav(id, origin);
  });

  // left/right walk the deck. Up/down and page keys are left alone — a long
  // chapter still scrolls.
  doc.addEventListener('keydown', function (e) {
    if (e.metaKey || e.ctrlKey || e.altKey || !cur) return;
    if (sheetPanel && sheetPanel.classList.contains('open')) return;
    const t = e.target;
    if (t && (t.isContentEditable || /^(INPUT|TEXTAREA|SELECT)$/.test(t.tagName))) return;
    let d = 0;
    if (e.key === 'ArrowRight') d = 1;
    else if (e.key === 'ArrowLeft') d = -1;
    else return;
    const n = cur.i + d;
    if (n < 0 || n >= panels.length) return;
    e.preventDefault();
    nav(panels[n].id, tabOrigin(), true);
  });

  /* Where the index currently sits on screen, so a keyboard page-turn pushes
     the dust from the same place a click would have. */
  function tabOrigin() {
    if (!cur) return null;
    for (let i = 0; i < cur.tabs.length; i++) {
      const a = cur.tabs[i];
      if (!a.offsetParent) continue;               // hidden at this breakpoint
      const r = a.getBoundingClientRect();
      return {
        x: (r.left + r.width / 2) / window.innerWidth,
        y: 1 - (r.top + r.height / 2) / window.innerHeight
      };
    }
    return null;
  }

  // ---------------------------------------------------------------- scroll

  let queued = false;

  function onScrollFrame() {
    queued = false;
    const y = window.scrollY;
    const dy = y - lastY;
    lastY = y;

    // scrolling inside a chapter shears the fluid: a tanh layer across
    // mid-screen, which is a genuine Kelvin-Helmholtz setup, so it rolls up
    if (field && Math.abs(dy) > 0.5 && !reduced) {
      field.shear(Math.max(-280, Math.min(280, dy * 3.0)));
    }
  }

  window.addEventListener('scroll', function () {
    if (!queued) { queued = true; requestAnimationFrame(onScrollFrame); }
  }, { passive: true });

  // ------------------------------------------------------------- pointer

  (function pointer() {
    if (!field) return;
    let px = null, py = null;
    let curX = 0, curY = 0, accX = 0, accY = 0, touch = false, pend = false;

    // Each splat is two full-grid passes in the solver, and a fast mouse can
    // deliver dozens of pointermove events between frames. Accumulate them and
    // emit at most one splat per 28 ms; visually it is the same stroke, because
    // the accumulator carries the whole displacement either way, and it removes
    // the worst cost spike on the page. The gate is a duration rather than one
    // per animation frame so that the field running at the display rate does
    // not lay ink down twice as fast on a 120 Hz screen as on a 60 Hz one.
    let lastFlush = 0;
    function flush() {
      pend = false;
      const now = performance.now();
      if (now - lastFlush < 28) { pend = true; requestAnimationFrame(flush); return; }
      lastFlush = now;
      const w = window.innerWidth, h = window.innerHeight;
      const sp = Math.hypot(accX, accY);
      if (sp > 0.0004) {
        field.splat(curX / w, 1 - curY / h, accX * 4200, accY * 4200, touch ? 0.0042 : 0.0030);
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
      if (!pend) { pend = true; requestAnimationFrame(flush); }
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
    if (cj && s.mgLevels) cj.textContent = String(s.mgLevels);
    if (fb) fb.textContent = field.kind === 'webgl2' ? 'WebGL2' : 'Canvas 2D fallback';
  })();

  // ------------------------------------------------- drop-in image slots

  function wireSlots(root) {
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

    $$('img[data-slot]', root).forEach(function (img) {
      if (img.getAttribute('data-failed') || (img.complete && img.naturalWidth === 0)) {
        card(img);
      } else if (!img.complete) {
        img.addEventListener('error', function () { card(img); });
      }
    });

    // any non-slot image that fails should not leave an empty hole
    $$('img:not([data-slot])', root).forEach(function (img) {
      img.addEventListener('error', function () { img.style.opacity = '0'; });
    });
  }

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

  bootThen(fromUrl);
})();

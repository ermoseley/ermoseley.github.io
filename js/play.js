/* =========================================================================
   play.js — the chapter that hands the solver over.

   The knobs belong to the engine. field.controls() says what exists, this file
   turns each descriptor into a real form control, and every change goes back
   through field.set() — whose return value, not the widget, is what gets
   displayed. So a knob added in field-mhd.js appears here on its own, and the
   panel has no second copy of the state to drift out of step with the box it is
   driving.
   ========================================================================= */

(function () {
  'use strict';

  const doc = document;
  const $ = (s, r) => (r || doc).querySelector(s);

  const mount = $('#play-controls');
  if (!mount) return;              // the chapter is optional, and this file is not load-bearing

  const table = $('#play-readout');
  const resetBtn = $('#play-reset');
  const said = $('#play-said');
  const note = $('#play-note');
  const panel = mount.closest ? mount.closest('.panel') : null;

  // site.js keeps its own instance private, so the engine publishes itself. On a
  // browser without WebGL2 or float render targets it never publishes at all.
  const field = window.__mhdField;

  function notice(msg) {
    mount.textContent = '';
    const p = doc.createElement('p');
    p.className = 'mono faint';
    p.textContent = msg;
    mount.appendChild(p);
  }

  /* Nothing on this page is allowed to depend on the simulation, this chapter
     included: it says so in one line and stops, rather than leaving a row of dead
     switches, a table of dashes and a reset button that resets nothing. */
  if (!field || typeof field.controls !== 'function' || typeof field.set !== 'function') {
    notice('The gas behind the page is not running here, so there is nothing to turn.');
    if (table && table.parentNode) table.parentNode.hidden = true;
    if (resetBtn) resetBtn.hidden = true;
    if (note) note.hidden = true;
    return;
  }

  // ------------------------------------------------------------------ numbers

  function fx(v, n) { return isFinite(v) ? Number(v).toFixed(n) : '—'; }
  function ints(v) { return isFinite(v) ? Math.round(v).toLocaleString() : '—'; }
  function str(v) { return (v === null || v === undefined || v === '') ? '—' : String(v); }

  // the same two-significant-figure exponent site.js prints the divergence in
  function sci(v) {
    if (!isFinite(v)) return '—';
    if (v === 0) return '0';
    const e = Math.floor(Math.log10(Math.abs(v)));
    const m = v / Math.pow(10, e);
    return m.toFixed(1) + 'e' + (e < 0 ? '' : '+') + e;
  }

  // how many decimals a slider should print, taken from its own step rather than
  // guessed, so a knob that arrives with step 0.01 reads to two places on its own
  function places(step) {
    if (!isFinite(step) || step <= 0) return 2;
    if (step >= 1) return 0;
    return Math.min(4, Math.max(1, Math.ceil(-Math.log10(step))));
  }

  // ------------------------------------------------------------------ storage

  /* sessionStorage, and deliberately nothing in the hash: site.js:284 sends any
     hash it does not recognise to the Cover, so a settings string written there
     would throw the reader to the front page mid-drag. One tab, one set of
     settings, gone when the tab closes. */
  const KEY = 'wiem-play';

  function load() {
    try {
      const raw = sessionStorage.getItem(KEY);
      const obj = raw ? JSON.parse(raw) : null;
      return (obj && typeof obj === 'object') ? obj : {};
    } catch (e) { return {}; }      // private mode, or a half-written value
  }

  let store = load();

  // a slider drag is dozens of changes a second and every write is synchronous,
  // so the settings are flushed on a trailing timer instead of per event
  let saveT = null;
  function save() {
    if (saveT) return;
    saveT = setTimeout(function () {
      saveT = null;
      try { sessionStorage.setItem(KEY, JSON.stringify(store)); } catch (e) { /* not fatal */ }
    }, 300);
  }

  function put(key, value) {
    try { return field.set(key, value); }
    catch (e) { return undefined; }  // a knob that refuses a value is its own business
  }

  /* Only what was actually changed is remembered. Persisting the whole list would
     freeze the engine's own defaults into the session, and a default later moved
     in field-mhd.js would never be seen again in this tab. */
  function apply(key, value) {
    const taken = put(key, value);
    if (taken === undefined) return undefined;
    store[key] = taken;
    save();
    return taken;
  }

  // ----------------------------------------------------------------- the CSS

  /* main.css styles no form control — the site has never had one — and a chapter
     does not get to add a stylesheet, so the handful of declarations the widgets
     need in order to sit on a dark page live here instead of being scattered
     through the builders below. Everything else is the page's own classes. */
  const CSS = {
    row:    'margin-top:1.15rem',
    label:  'display:flex;justify-content:space-between;align-items:baseline;gap:1rem;' +
            'font-family:var(--mono);font-size:0.6rem;letter-spacing:0.14em;' +
            'text-transform:uppercase;color:var(--ink-faint);cursor:pointer',
    value:  'color:var(--accent);font-weight:500;letter-spacing:0.06em',
    range:  'width:100%;display:block;margin:0.55rem 0 0;accent-color:var(--accent)',
    check:  'width:0.95rem;height:0.95rem;margin:0;accent-color:var(--accent);cursor:pointer',
    select: 'width:100%;display:block;margin:0.5rem 0 0;padding:0.42rem 0.5rem;' +
            'font-family:var(--mono);font-size:0.66rem;letter-spacing:0.06em;' +
            'color:var(--ink);background:rgba(6,7,12,0.62);border:1px solid var(--rule)',
    note:   'margin:0.5rem 0 0;font-size:0.78rem;line-height:1.5;color:var(--ink-faint)'
  };

  // ---------------------------------------------------------------- controls

  const wired = {};                  // key -> a function that repaints that widget
  const boxes = {};                  // group -> the element its rows go into

  // A name for each group the engine currently has. An unknown group is not an
  // error: it gets its own heading under its own name, which is the point.
  const GROUP = {
    gas: 'The gas',
    solver: 'The solver',
    field: 'The field',
    dust: 'The dust',
    display: 'What is drawn'
  };

  // The groups appear in the order the engine first mentions them, so the display
  // order in controls() is the display order here.
  function groupOf(name) {
    const g = String(name === null || name === undefined ? 'other' : name);
    if (boxes[g]) return boxes[g];
    const box = doc.createElement('div');
    box.className = 'prose';         // .prose h3 is the page's own section rule
    const h = doc.createElement('h3');
    h.textContent = GROUP[g] || g;
    box.appendChild(h);
    mount.appendChild(box);
    boxes[g] = box;
    return box;
  }

  function shell(d, id) {
    const box = doc.createElement('div');
    box.setAttribute('style', CSS.row);
    const lab = doc.createElement('label');
    lab.setAttribute('for', id);
    lab.setAttribute('style', CSS.label);
    const name = doc.createElement('span');
    name.textContent = d.label || d.key;
    lab.appendChild(name);
    box.appendChild(lab);

    // A note can be computed -- the engine's line under the plasma beta slider
    // quotes the Alfven speed the current beta implies -- so it is repainted with
    // the rest of the row rather than written once and left to go stale. It is
    // appended last, after whatever the builder puts in.
    let p = null;
    function note(n) {
      const txt = n.note || '';
      if (!txt) {
        if (p) { box.removeChild(p); p = null; }
        return;
      }
      if (!p) {
        p = doc.createElement('p');
        p.setAttribute('style', CSS.note);
        box.appendChild(p);
      }
      if (p.textContent !== txt) p.textContent = txt;
    }

    return { box: box, lab: lab, note: note };
  }

  function toggle(d, id) {
    const s = shell(d, id);
    const el = doc.createElement('input');
    el.type = 'checkbox';
    el.id = id;
    el.setAttribute('style', CSS.check);
    el.checked = !!d.value;
    s.lab.appendChild(el);
    s.note(d);

    el.addEventListener('change', function () {
      const taken = apply(d.key, el.checked);
      if (taken !== undefined) el.checked = !!taken;
      resync(d.key);
    });

    return {
      box: s.box,
      sync: function (n, self) { if (!self) el.checked = !!n.value; s.note(n); }
    };
  }

  function select(d, id) {
    const s = shell(d, id);
    const el = doc.createElement('select');
    el.id = id;
    el.setAttribute('style', CSS.select);
    s.box.appendChild(el);

    // the option list is rebuilt only when it has actually moved, because
    // replacing it under an open dropdown closes it
    let sig = null;
    function fill(options, value) {
      const opts = options || [];
      const next = opts.map(function (o) { return String(o && o.value); }).join('');
      if (next !== sig) {
        sig = next;
        el.textContent = '';
        opts.forEach(function (o) {
          const op = doc.createElement('option');
          op.value = String(o && o.value);
          op.textContent = (o && o.label) || String(o && o.value);
          el.appendChild(op);
        });
      }
      if (value !== null && value !== undefined) el.value = String(value);
    }

    fill(d.options, d.value);
    s.note(d);

    el.addEventListener('change', function () {
      const taken = apply(d.key, el.value);
      if (taken !== undefined && taken !== null) el.value = String(taken);
      resync(d.key);
    });

    return {
      box: s.box,
      // the option list is still taken, in case the change moved it; the value is
      // not, because on the knob that was just turned set() has already had the
      // last word on that
      sync: function (n, self) { fill(n.options, self ? null : n.value); s.note(n); }
    };
  }

  function range(d, id) {
    const s = shell(d, id);
    const out = doc.createElement('b');
    out.setAttribute('style', CSS.value);
    s.lab.appendChild(out);

    const el = doc.createElement('input');
    el.type = 'range';
    el.id = id;
    el.setAttribute('style', CSS.range);
    s.box.appendChild(el);

    let dp = places(d.step);
    function show(v) { out.textContent = isFinite(v) ? Number(v).toFixed(dp) : '—'; }

    function bounds(n) {
      if (n.min !== null && n.min !== undefined) el.min = String(n.min);
      if (n.max !== null && n.max !== undefined) el.max = String(n.max);
      if (n.step !== null && n.step !== undefined) { el.step = String(n.step); dp = places(n.step); }
    }

    bounds(d);
    el.value = String(d.value);
    show(d.value);
    s.note(d);

    function send() {
      const want = parseFloat(el.value);
      const taken = apply(d.key, want);
      const v = (taken === undefined || !isFinite(taken)) ? want : taken;
      show(v);
      // the thumb is only moved when the engine really did disagree: writing the
      // value back on every event would snap it out from under the finger on a
      // knob whose own step is coarser than the slider's
      if (Math.abs(v - want) > 1e-9) el.value = String(v);
      resync(d.key);
    }

    /* A descriptor that says live: false is a knob whose set() costs something an
       order heavier than a frame -- reallocating the grid, reseeding the field --
       so the number under the label follows the drag and the engine only hears
       about it when the drag stops. Anything that does not say so is applied as it
       moves, which is the point of a slider. */
    const live = d.live !== false;
    el.addEventListener('input', function () {
      if (live) send(); else show(parseFloat(el.value));
    });
    if (!live) el.addEventListener('change', send);

    return {
      box: s.box,
      sync: function (n, self) {
        bounds(n);
        if (!self && isFinite(n.value)) { el.value = String(n.value); show(n.value); }
        s.note(n);
      }
    };
  }

  function read() {
    let l;
    try { l = field.controls(); } catch (e) { return []; }
    return Array.isArray(l) ? l : [];
  }

  /* The value set() handed back is authoritative for the knob that was just
     turned, so that one keeps the number it was given; every other row is
     repainted from a fresh list, because one setting can move another — turning
     the field off has nothing to say about its own switch and everything to say
     about the two that draw it. The turned knob is still repainted in every other
     respect, since its own note may quote what it just did.
     Called with no argument -- after a reset -- nothing is exempt. */
  function resync(turned) {
    read().forEach(function (d) {
      const fn = d && wired[d.key];
      if (fn) fn(d, d.key === turned);
    });
  }

  let list = read();

  // The engine's own values, read before the session's settings go anywhere near
  // it, are what the reset button restores.
  const defaults = {};
  list.forEach(function (d) { if (d && d.key) defaults[d.key] = d.value; });

  (function restore() {
    const known = {};
    list.forEach(function (d) { if (d && d.key) known[d.key] = true; });
    let moved = false;
    Object.keys(store).forEach(function (k) {
      moved = true;
      if (!known[k]) { delete store[k]; return; }  // a knob that no longer exists
      const taken = put(k, store[k]);
      if (taken === undefined) delete store[k];
      else store[k] = taken;                      // a stored value the engine has since clamped
    });
    if (moved) save();                            // a first visit writes nothing at all
    list = read();                                // applying one setting can move another
  })();

  if (!list.length) {
    // the engine is there and measurable, it just has nothing on offer, so the
    // numbers below stay
    notice('The gas is running, but it is not publishing any controls.');
  } else {
    mount.textContent = '';
    list.forEach(function (d) {
      if (!d || !d.key) return;
      // site.js reaches for a lot of single ids by querySelector, so everything
      // this file puts in the document is namespaced out of its way
      const id = 'play-c-' + String(d.key).replace(/[^A-Za-z0-9_-]/g, '-');
      const made = d.kind === 'toggle' ? toggle(d, id)
        : d.kind === 'select' ? select(d, id)
          : d.kind === 'range' ? range(d, id)
            : null;
      if (!made) return;           // an unfamiliar kind is left out rather than guessed at
      groupOf(d.group).appendChild(made.box);
      wired[d.key] = made.sync;
    });
  }

  /* Nothing built here carries .rv. The cascade is armed when a chapter opens and
     a .rv element is invisible until it is reached, so a row created after this
     chapter was first entered would sit at zero opacity for good. The static
     wrappers in the markup are revealed instead, and their contents ride along. */

  let saidT = null;

  /* The engine undoes itself if it knows how, because it can do more than pushing
     the old numbers back would: a value that has gone through set() is an
     override, and an override outlives the reset, so those knobs would stop
     following the chapter presets for the rest of the session. Where there is no
     such method the captured defaults are the honest best. */
  function restoreDefaults() {
    if (typeof field.resetConfig === 'function') {
      try { field.resetConfig(); return; } catch (e) { /* fall through */ }
    }
    Object.keys(defaults).forEach(function (k) { put(k, defaults[k]); });
  }

  if (resetBtn) resetBtn.addEventListener('click', function () {
    restoreDefaults();
    store = {};
    clearTimeout(saveT);             // or a flush still in flight writes the old settings back
    saveT = null;
    try { sessionStorage.removeItem(KEY); } catch (e) { /* not fatal */ }
    resync();
    if (said) {
      said.textContent = 'defaults restored';
      clearTimeout(saidT);
      saidT = setTimeout(function () { said.textContent = ''; }, 2600);
    }
  });

  // ----------------------------------------------------------------- readout

  /* Labels are this file's own strings, so they are written as markup and get the
     subscripts the physics wants; every value is set as text, because it comes
     from the engine. */
  const ROWS = [
    ['grid', function (r) {
      if (!Array.isArray(r.grid) || !isFinite(r.grid[0])) return '—';
      return r.grid[0] + ' × ' + r.grid[1] +
        '  (' + (r.grid[0] * r.grid[1]).toLocaleString() + ' cells)';
    }],
    ['fps', function (r) { return fx(r.fps, 0); }],
    ['steps', function (r) { return ints(r.steps); }],
    ['&#8499;<sub>rms</sub> · &#8499;<sub>max</sub>',
      function (r) { return fx(r.rms, 2) + ' · ' + fx(r.machMax, 2); }],
    ['c<sub>tot</sub>', function (r) { return fx(r.ctot, 3); }],
    ['&Delta;t / &Delta;x', function (r) { return fx(r.dtdx, 4); }],
    ['mass', function (r) { return fx(r.dens, 5); }],
    ['scheme', function (r) { return str(r.solver); }],
    ['reconstruction, last step', function (r) { return String(str(r.recon)).toUpperCase(); }],
    ['plasma &beta;', function (r) { return fx(r.beta, 2); }],
    ['|B<sub>0</sub>| · &lang;|B|&rang;',
      function (r) { return fx(r.b0, 3) + ' · ' + fx(r.bmean, 3); }],
    ['&#8499;<sub>A</sub>', function (r) { return fx(r.alfvenMach, 2); }],
    ['grain <i>q/m</i>', function (r) { return fx(r.charge, 0); }],
    ['&nabla;&middot;B rms · max',
      function (r) { return sci(r.divbRms) + ' · ' + sci(r.divbMax); }]
  ];

  const cells = [];

  function paint() {
    // the deck holds the other twelve chapters out of the document at any moment,
    // and a readout nobody can see has no business reading the solver either
    if (!cells.length || (panel && !panel.classList.contains('is-on'))) return;
    let r;
    try { r = field.readout(); } catch (e) { return; }
    if (!r) return;
    cells.forEach(function (c) {
      let v;
      try { v = c[1](r); } catch (e) { v = '—'; }
      if (c[0].textContent !== v) c[0].textContent = v;
    });
  }

  if (table && typeof field.readout === 'function') {
    table.textContent = '';
    ROWS.forEach(function (spec) {
      const tr = doc.createElement('tr');
      const th = doc.createElement('th');
      th.innerHTML = spec[0];
      const td = doc.createElement('td');
      td.textContent = '—';
      tr.appendChild(th);
      tr.appendChild(td);
      table.appendChild(tr);
      cells.push([td, spec[1]]);
    });

    /* About three times a second, as a setTimeout chain rather than setInterval —
       the pattern at site.js:471. Each tick is scheduled by the one before it, so
       a frame the browser took a long time over cannot leave a queue of ticks
       waiting to fire back to back. */
    (function tick() {
      paint();
      setTimeout(tick, 320);
    })();
  } else if (table) {
    const tr = doc.createElement('tr');
    tr.innerHTML = '<th>readout</th><td>not published by this build</td>';
    table.textContent = '';
    table.appendChild(tr);
  }
})();

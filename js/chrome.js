/* =========================================================================
   chrome.js — two switches in the corner, and the page furniture behind them.

   Both are the reader's, not the page's, and that is the whole design brief:
   the simulation can run behind the type, which is lovely until it is not, and
   the answer to "this is a bit much" should not be to leave. It opens still so
   the reader opts into both the motion and the work.

     freeze — stops the gas. Held in the engine (field.setPaused) rather than
              here, because the play chapter has a view of the same state and
              two booleans that can disagree is a bug waiting to be filed.
     bare   — takes the type away and leaves the box. Everything except the
              canvas and this cluster.

   The cluster stays visible in bare mode. It has to: it is the way back, and a
   full-screen simulation with no affordance is a trap, not a feature.

   Top right, because the scrim in main.css is weighted to the left where the
   type lives, the rail is on the left, and .strip owns the bottom edge on a
   narrow screen. Low contrast until hovered or focused, and a little brighter
   in bare mode where there is nothing else to look at.

   No dependency on the engine: with no WebGL2 the freeze switch hides itself and
   bare mode still works, because bare mode is only ever a class on <body>.
   ========================================================================= */

(function (global) {
  'use strict';

  const doc = global.document;
  const KEY = 'wiem-chrome';
  const TOUR_KEY = 'wiem-tour-v1';
  const reduced = global.matchMedia && global.matchMedia('(prefers-reduced-motion: reduce)').matches;

  // ------------------------------------------------------------------ storage

  /* sessionStorage, for the same reason play.js uses it: site.js sends any hash it
     does not recognise to the Cover, so state in the URL would throw the reader to
     the front page. One tab, one choice, gone when the tab closes. */
  function load() {
    try {
      const raw = sessionStorage.getItem(KEY);
      const o = raw ? JSON.parse(raw) : null;
      return (o && typeof o === 'object') ? o : {};
    } catch (e) { return {}; }
  }

  function save(o) {
    try { sessionStorage.setItem(KEY, JSON.stringify(o)); } catch (e) {}
  }

  function loadTour() {
    try {
      const stage = sessionStorage.getItem(TOUR_KEY);
      return stage === 'play' || stage === 'done' ? stage : '';
    } catch (e) { return ''; }
  }

  function saveTour(stage) {
    try { sessionStorage.setItem(TOUR_KEY, stage); } catch (e) {}
  }

  const store = load();
  let tourStage = loadTour();
  // A fresh tab always opens on the painted initial state. Once the introduction
  // has been answered, reloads in that tab respect the reader's last choice.
  let paused = tourStage ? !!store.paused : true;
  let bare = !!store.bare;
  if (!tourStage) {
    paused = true;
    bare = false;
    store.paused = true;
    store.bare = false;
    save(store);
  }

  // ------------------------------------------------------------------- styling

  const CSS = `
#corner {
  position: fixed;
  top: calc(env(safe-area-inset-top, 0px) + 0.85rem);
  right: calc(env(safe-area-inset-right, 0px) + 0.85rem);
  z-index: 50;
  display: flex;
  gap: 0.3rem;
  font-family: var(--mono);
  font-size: 0.55rem;
  letter-spacing: 0.16em;
  text-transform: uppercase;
}
#corner button {
  appearance: none;
  -webkit-appearance: none;
  cursor: pointer;
  margin: 0;
  padding: 0.36rem 0.5rem;
  min-width: 3.1rem;
  color: var(--ink-faint);
  background: rgba(6, 7, 12, 0.5);
  border: 1px solid var(--rule-soft);
  border-radius: 2px;
  backdrop-filter: blur(10px);
  -webkit-backdrop-filter: blur(10px);
  font: inherit;
  letter-spacing: inherit;
  text-transform: inherit;
  opacity: 0.42;
  transition: opacity 220ms var(--ease), color 220ms var(--ease),
              border-color 220ms var(--ease);
}
#corner button:hover,
#corner button:focus-visible {
  opacity: 1;
  color: var(--ink);
  border-color: var(--rule);
}
#corner button:focus-visible {
  outline: 1px solid var(--accent);
  outline-offset: 1px;
}
/* On when engaged, so the reader can see at a glance that the box is held. */
#corner button[aria-pressed="true"] {
  opacity: 0.92;
  color: var(--accent);
  border-color: var(--rule);
}

/* ------------------------------------------------------------- first-visit tour */
#sim-tour {
  position: fixed;
  inset: 0;
  z-index: 49;
  display: flex;
  align-items: flex-start;
  justify-content: flex-end;
  min-height: 100vh;
  min-height: 100svh;
  padding: calc(env(safe-area-inset-top, 0px) + 4.9rem)
           calc(env(safe-area-inset-right, 0px) + 0.85rem)
           calc(env(safe-area-inset-bottom, 0px) + 4.5rem)
           calc(env(safe-area-inset-left, 0px) + 0.85rem);
  background: rgba(2, 3, 8, 0.72);
  backdrop-filter: blur(4px);
  -webkit-backdrop-filter: blur(4px);
  opacity: 0;
  visibility: hidden;
  overflow-y: auto;
  pointer-events: none;
  transition: opacity 360ms var(--ease), visibility 360ms;
}
#sim-tour.is-visible {
  opacity: 1;
  visibility: visible;
  pointer-events: auto;
}
#sim-tour.tour-play {
  z-index: 44;
  align-items: flex-end;
  justify-content: center;
  padding-bottom: calc(env(safe-area-inset-bottom, 0px) + 5.25rem);
  pointer-events: none;
}
#sim-tour .tour-card {
  position: relative;
  width: min(42rem, calc(100vw - 1.7rem));
  padding: clamp(1.15rem, 3.2vw, 2rem);
  color: #fff;
  background: rgba(5, 6, 11, 0.50);
  border: 1px solid rgba(255, 255, 255, 0.18);
  border-radius: 3px;
  box-shadow: 0 1.4rem 4rem rgba(0, 0, 0, 0.34);
  transform: translateY(10px);
  transition: transform 560ms var(--ease);
}
#sim-tour.is-visible .tour-card { transform: none; }
#sim-tour .tour-eyebrow {
  margin: 0 0 0.8rem;
  color: var(--accent);
  font-family: var(--mono);
  font-size: 0.57rem;
  letter-spacing: 0.2em;
  text-transform: uppercase;
}
#sim-tour .tour-copy {
  margin: 0;
  max-width: 38rem;
  color: #fff;
  font-family: var(--sans);
  font-size: clamp(1rem, 1.25vw, 1.17rem);
  font-weight: 300;
  line-height: 1.65;
}
#sim-tour .tour-copy + .tour-copy { margin-top: 0.8rem; }
#sim-tour .tour-copy strong { font-weight: 500; }
#sim-tour .tour-actions {
  display: flex;
  flex-wrap: wrap;
  gap: 0.65rem;
  margin-top: 1.35rem;
}
#sim-tour .tour-actions button {
  appearance: none;
  -webkit-appearance: none;
  cursor: pointer;
  min-height: 2.7rem;
  padding: 0.72rem 0.95rem;
  color: #fff;
  background: rgba(255, 255, 255, 0.05);
  border: 1px solid rgba(255, 255, 255, 0.24);
  border-radius: 2px;
  font-family: var(--mono);
  font-size: 0.59rem;
  letter-spacing: 0.13em;
  text-transform: uppercase;
  transition: color 180ms var(--ease), background 180ms var(--ease),
              border-color 180ms var(--ease);
}
#sim-tour .tour-actions button[data-tour="start"] {
  color: #05060a;
  background: #fff;
  border-color: #fff;
}
#sim-tour .tour-actions button:hover,
#sim-tour .tour-actions button:focus-visible {
  color: #05060a;
  background: var(--accent);
  border-color: var(--accent);
  outline: none;
}

/* The lines terminate at the live controls rather than at decorative copies. */
#sim-tour.tour-start .tour-card::before {
  content: '';
  position: absolute;
  top: -2.55rem;
  right: 4.2rem;
  width: 1px;
  height: 2rem;
  background: rgba(255, 255, 255, 0.62);
}
#sim-tour.tour-start .tour-card::after {
  content: '';
  position: absolute;
  top: -2.65rem;
  right: calc(4.2rem - 3px);
  width: 7px;
  height: 7px;
  border-top: 1px solid #fff;
  border-left: 1px solid #fff;
  transform: rotate(45deg);
}
#sim-tour.tour-play .tour-card { width: min(34rem, calc(100vw - 1.7rem)); }
#sim-tour.tour-play .tour-card::after {
  content: '';
  position: absolute;
  bottom: -1.9rem;
  left: 50%;
  width: 1px;
  height: 1.35rem;
  background: rgba(255, 255, 255, 0.66);
}

body.is-tour-start #corner button {
  opacity: 0.12;
  pointer-events: none;
}
body.is-tour-start #corner #freeze-button {
  opacity: 1;
  pointer-events: auto;
  color: #fff;
  border-color: var(--accent);
  box-shadow: 0 0 0 1px rgba(185, 198, 255, 0.16),
              0 0 1.6rem rgba(185, 198, 255, 0.48);
}
body.is-tour-play #corner {
  opacity: 0.12;
}
body.is-tour-play .strip a {
  opacity: 0.14;
}
body.is-tour-play .strip #tab-strip-play {
  opacity: 1;
  pointer-events: auto;
  color: #fff;
  background: rgba(255, 255, 255, 0.10);
  box-shadow: inset 0 2px 0 var(--accent), 0 0 1.7rem rgba(185, 198, 255, 0.35);
}

@media (max-width: 620px) {
  #sim-tour {
    padding-top: calc(env(safe-area-inset-top, 0px) + 4.75rem);
  }
  #sim-tour .tour-card {
    padding: 1.15rem;
  }
  #sim-tour .tour-actions { flex-direction: column; }
  #sim-tour .tour-actions button { width: 100%; }
}

/* ------------------------------------------------------------------- bare mode
   Everything but the canvas and the cluster. Visibility, not display, and that
   choice is load-bearing three times over: the page keeps its height so the reader
   comes back to the paragraph they left, site.js's IntersectionObserver keeps
   firing so the router does not lose its place, and scrolling still reaches the gas
   -- the scroll shear is one of the interactions, and it would be a shame to take
   it away from the one view where the box is all there is.

   <footer> lives inside <main id="shell">, so it is covered by the first selector
   rather than needing its own. */
body.is-bare #shell,
body.is-bare .rail,
body.is-bare .strip,
body.is-bare .sheet,
body.is-bare .menu-btn {
  visibility: hidden !important;
}
body.is-bare #shell { pointer-events: none; }
/* The legibility scrim exists for type that is no longer there. */
body.is-bare #veil { opacity: 0; }
#veil { transition: opacity 420ms var(--ease); }

@media (prefers-reduced-motion: reduce) {
  #corner button, #veil, #sim-tour, #sim-tour .tour-card { transition: none; }
}
`;

  function styles() {
    const s = doc.createElement('style');
    s.id = 'corner-css';
    s.textContent = CSS;
    doc.head.appendChild(s);
  }

  // -------------------------------------------------------------------- engine

  // site.js keeps its instance private, so the engine publishes itself. Looked up
  // lazily: this file may well be parsed before field-mhd.js has created one.
  function field() {
    const f = global.__mhdField;
    return (f && typeof f.setPaused === 'function') ? f : null;
  }

  // ------------------------------------------------------------------- the pair

  const listeners = [];
  function announce() { for (const fn of listeners) { try { fn(); } catch (e) {} } }

  let bFreeze = null, bBare = null;

  function paintButtons() {
    if (bFreeze) {
      bFreeze.setAttribute('aria-pressed', paused ? 'true' : 'false');
      bFreeze.textContent = paused ? 'frozen' : 'freeze';
      bFreeze.title = paused ? 'Let the simulation run again' : 'Hold the simulation still';
    }
    if (bBare) {
      bBare.setAttribute('aria-pressed', bare ? 'true' : 'false');
      bBare.textContent = bare ? 'page' : 'bare';
      bBare.title = bare ? 'Bring the page back' : 'Hide everything but the simulation';
    }
  }

  function applyPause() {
    const f = field();
    if (f) f.setPaused(paused);
  }

  function applyBare() {
    doc.body.classList.toggle('is-bare', bare);
  }

  function setPaused(v, quiet) {
    paused = !!v;
    store.paused = paused;
    save(store);
    applyPause();
    paintButtons();
    if (!quiet) announce();
    return paused;
  }

  function setBare(v, quiet) {
    bare = !!v;
    store.bare = bare;
    save(store);
    applyBare();
    paintButtons();
    if (!quiet) announce();
    return bare;
  }

  // ----------------------------------------------------------- first-visit tour

  let tour = null;
  let tourTimer = null;
  let tourStarted = false;
  let bootTries = 0;
  let bootSettled = false;
  let playTarget = null;
  let stripLeft = null;
  let previousFocus = null;

  function tourShell() {
    if (tour) return tour;
    tour = doc.createElement('div');
    tour.id = 'sim-tour';
    tour.setAttribute('role', 'dialog');
    tour.setAttribute('aria-modal', 'true');
    tour.addEventListener('click', function (e) {
      const b = e.target.closest ? e.target.closest('[data-tour]') : null;
      if (!b) return;
      if (b.dataset.tour === 'start') {
        setPaused(false);
        advanceTour();
      } else if (b.dataset.tour === 'decline') {
        setPaused(true);
        advanceTour();
      }
    });
    doc.body.appendChild(tour);
    return tour;
  }

  function revealTour(focus) {
    requestAnimationFrame(function () {
      if (!tour) return;
      tour.classList.add('is-visible');
      if (focus) setTimeout(function () {
        if (focus.isConnected) focus.focus({ preventScroll: true });
      }, reduced ? 0 : 280);
    });
  }

  function showStartTour() {
    const el = tourShell();
    previousFocus = doc.activeElement;
    el.className = 'tour-start';
    el.setAttribute('aria-label', 'Start the interactive background simulation');
    el.innerHTML =
      '<div class="tour-card">' +
        '<p class="tour-eyebrow">Interactive background</p>' +
        '<p class="tour-copy" id="sim-tour-copy">Click here to begin an interactive background simulation with real ISM physics. ' +
          'The physics you see here is the same <strong>MUSCL-Hancock/PLM/HLLD</strong> algorithm as you see in ' +
          '<strong>RAMSES</strong> and <strong>mini-RAMSES</strong>.</p>' +
        '<p class="tour-copy">If at any time you decide it\'s all a bit much, just click here again.</p>' +
        '<div class="tour-actions">' +
          '<button type="button" data-tour="start">Start the simulation</button>' +
          '<button type="button" data-tour="decline">No, thank you</button>' +
        '</div>' +
      '</div>';
    doc.body.classList.remove('is-tour-play');
    doc.body.classList.add('is-tour-start');
    if (bFreeze) bFreeze.setAttribute('aria-describedby', 'sim-tour-copy');
    revealTour(el.querySelector('[data-tour="start"]'));
  }

  function onTourClickAway() { finishTour(false); }

  function showPlayTour() {
    const el = tourShell();
    playTarget = doc.getElementById('tab-strip-play');
    if (!playTarget) { finishTour(false); return; }

    el.className = 'tour-play';
    el.setAttribute('aria-modal', 'false');
    el.setAttribute('aria-label', 'Find the simulation controls');
    el.innerHTML =
      '<div class="tour-card">' +
        '<p class="tour-eyebrow">The solver is yours</p>' +
        '<p class="tour-copy" id="sim-tour-copy">Go here to modify the physics and the details of the solver.</p>' +
      '</div>';
    doc.body.classList.remove('is-tour-start');
    doc.body.classList.add('is-tour-play');
    if (bFreeze) bFreeze.removeAttribute('aria-describedby');

    const strip = playTarget.parentNode;
    if (stripLeft === null) stripLeft = strip.scrollLeft;
    const max = Math.max(0, strip.scrollWidth - strip.clientWidth);
    const left = Math.max(0, Math.min(max,
      playTarget.offsetLeft - strip.clientWidth / 2 + playTarget.offsetWidth / 2));
    if (strip.scrollTo) strip.scrollTo({ left: left, behavior: reduced ? 'auto' : 'smooth' });
    else strip.scrollLeft = left;

    playTarget.setAttribute('aria-describedby', 'sim-tour-copy');
    // This stage is a pointer, not a gate. Capture the reader's next click only
    // to dismiss the cue; the click itself still reaches whatever they chose.
    doc.addEventListener('click', onTourClickAway, { capture: true, once: true });
    revealTour(playTarget);
  }

  function advanceTour() {
    if (tourStage === 'play' || tourStage === 'done') return;
    tourStage = 'play';
    saveTour(tourStage);
    if (tour) tour.classList.remove('is-visible');
    doc.body.classList.remove('is-tour-start');
    if (bFreeze) bFreeze.removeAttribute('aria-describedby');
    clearTimeout(tourTimer);
    tourTimer = setTimeout(showPlayTour, reduced ? 0 : 400);
  }

  function finishTour(restoreStrip) {
    tourStage = 'done';
    saveTour(tourStage);
    clearTimeout(tourTimer);
    doc.body.classList.remove('is-tour-start', 'is-tour-play');
    if (bFreeze) bFreeze.removeAttribute('aria-describedby');
    if (playTarget) {
      const strip = playTarget.parentNode;
      playTarget.removeAttribute('aria-describedby');
      if (restoreStrip && stripLeft !== null) {
        if (strip.scrollTo) strip.scrollTo({ left: stripLeft, behavior: reduced ? 'auto' : 'smooth' });
        else strip.scrollLeft = stripLeft;
      }
    }
    doc.removeEventListener('click', onTourClickAway, true);
    const old = tour;
    tour = null;
    if (old) {
      old.classList.remove('is-visible');
      setTimeout(function () { if (old.parentNode) old.parentNode.removeChild(old); }, reduced ? 0 : 380);
    }
    if (restoreStrip && previousFocus && previousFocus.focus && previousFocus.isConnected) {
      previousFocus.focus({ preventScroll: true });
    }
  }

  function beginTour() {
    if (tourStarted || tourStage === 'done') return;
    const boot = doc.getElementById('boot');
    if (boot && !boot.classList.contains('done') && ++bootTries < 36) {
      setTimeout(beginTour, 100);
      return;
    }
    // Let the opaque boot readout finish fading before the introduction itself
    // begins to appear; otherwise the gentle entrance happens behind it.
    if (boot && !bootSettled) {
      bootSettled = true;
      setTimeout(beginTour, reduced ? 0 : 640);
      return;
    }
    tourStarted = true;
    if (tourStage === 'play') showPlayTour();
    else showStartTour();
  }

  function toggleFreeze() {
    const introducing = doc.body.classList.contains('is-tour-start');
    setPaused(!paused);
    if (introducing) advanceTour();
  }

  function wireFreeze(b) {
    b.id = 'freeze-button';
    b.addEventListener('click', toggleFreeze);
  }

  // ---------------------------------------------------------------------- build

  function button(label) {
    const b = doc.createElement('button');
    b.type = 'button';
    b.textContent = label;
    b.setAttribute('aria-pressed', 'false');
    return b;
  }

  function build() {
    styles();

    const box = doc.createElement('div');
    box.id = 'corner';

    bFreeze = button('freeze');
    wireFreeze(bFreeze);

    bBare = button('bare');
    bBare.id = 'bare-button';
    bBare.addEventListener('click', function () { setBare(!bare); });

    // With no engine there is nothing to freeze, so that switch does not appear at
    // all rather than sitting there doing nothing.
    if (field()) box.appendChild(bFreeze); else bFreeze = null;
    box.appendChild(bBare);
    doc.body.appendChild(box);

    applyBare();
    paintButtons();

    /* The engine may not exist yet. Poll briefly for it -- the same pattern the
       footer's B/S hint uses -- and adopt a stored freeze the moment it turns up, so
       a reload with the switch on comes back still. */
    let tries = 0;
    (function wait() {
      if (field()) {
        if (!bFreeze) {
          bFreeze = button('freeze');
          wireFreeze(bFreeze);
          box.insertBefore(bFreeze, box.firstChild);
        }
        applyPause();
        paintButtons();
        beginTour();
        return;
      }
      if (++tries < 60) setTimeout(wait, 100);
    })();

    /* Escape declines the motion, then dismisses the pointer if it is pressed again.
       Outside the introduction it keeps its original job of leaving bare mode. */
    doc.addEventListener('keydown', function (e) {
      if (e.key !== 'Escape') return;
      if (doc.body.classList.contains('is-tour-start')) {
        e.preventDefault();
        setPaused(true);
        advanceTour();
      } else if (doc.body.classList.contains('is-tour-play')) {
        e.preventDefault();
        finishTour(true);
      } else if (bare) {
        setBare(false);
      }
    });
  }

  if (doc.readyState === 'loading') doc.addEventListener('DOMContentLoaded', build);
  else build();

  /* The play chapter drives the same two switches, so it needs a handle rather than
     its own copy. onChange lets it repaint when the corner is used instead. */
  global.__pageChrome = {
    paused: function () { return paused; },
    setPaused: setPaused,
    bare: function () { return bare; },
    setBare: setBare,
    onChange: function (fn) { if (typeof fn === 'function') listeners.push(fn); }
  };

})(window);

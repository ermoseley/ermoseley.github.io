/* =========================================================================
   chrome.js — two switches in the corner, and the page furniture behind them.

   Both are the reader's, not the page's, and that is the whole design brief:
   the simulation is always running behind the type, which is lovely until it
   is not, and the answer to "this is a bit much" should not be to leave.

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

  const store = load();
  let paused = !!store.paused;
  let bare = !!store.bare;

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
  #corner button, #veil { transition: none; }
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
    bFreeze.addEventListener('click', function () { setPaused(!paused); });

    bBare = button('bare');
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
          bFreeze.addEventListener('click', function () { setPaused(!paused); });
          box.insertBefore(bFreeze, box.firstChild);
        }
        applyPause();
        paintButtons();
        return;
      }
      if (++tries < 60) setTimeout(wait, 100);
    })();

    /* Escape leaves bare mode. site.js binds Escape too, to close the nav sheet, but
       that sheet is hidden in bare mode and closing a closed sheet does nothing, so
       the two handlers do not fight. */
    doc.addEventListener('keydown', function (e) {
      if (e.key === 'Escape' && bare) setBare(false);
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

# whoisericmoseley

Personal site for **Eric R. Moseley** — numerical astrophysicist, KIPAC Fellow at Stanford/SLAC.

**Live:** https://ermoseley.github.io/ — this is the main page.

Served from the `site` branch of `ermoseley/ermoseley.github.io`. The Quarto project on that
repo's `main` branch and its build on `gh-pages` are both untouched and still buildable; the
site it produced is preserved verbatim under `legacy/` and reachable at
https://ermoseley.github.io/legacy/. To revert, point Pages back at `gh-pages`.

The other agents' attempts remain live at their own project-page URLs and are unaffected:
https://ermoseley.github.io/whoisericmoseley/ and
https://ermoseley.github.io/whoisericmoseley-opus/ (a staging copy of this site).

Push `main` here to update the source repo; push `main:site` to
`ermoseley/ermoseley.github.io` to publish.

No build step, no dependencies, no framework. Static HTML, one stylesheet, three scripts.
Push to `main` and GitHub Pages redeploys in under a minute.

---

## The background is a real simulation

`js/field.js` is a two-dimensional incompressible Navier–Stokes solver coupled to a
population of Lagrangian dust superparticles, written from scratch in WebGL2. Per frame:

| Stage | Method |
|---|---|
| Gas advection | Semi-Lagrangian backtrace, bilinear, periodic |
| Small-scale swirl | Fedkiw vorticity confinement, `f = ε ω (N_y, −N_x)`, `N = ∇\|ω\|/\|∇\|ω\|\|` |
| Incompressibility | Divergence → 16-iteration Jacobi pressure solve → gradient subtract, with a smooth speed governor at `\|u\| = 420` cells/s |
| Dust | `dv/dt = (u_gas − v)/τ_s`, advanced by **first-order implicit (backward Euler)**: `v ← (v + a·u)/(1+a)`, `a = dt/τ`. The population is **polydisperse** — each grain draws τ from a log-uniform spectrum 1.4 decades wide, hashed out of its own texel rather than stored — so the median grain's `dt/τ` is the only scalar the CPU supplies and a grain costs one divide, one multiply-add and one multiply, still L-stable |
| Dust rendering | One additive soft sprite per grain, at native resolution, straight to the screen. Brightness and size are set by the **drift** `\|v − u_gas\|` — the only quantity in the drag law that does work on a grain — so tightly coupled dust stays dim and decoupled grains light up where they pile up |
| Gas rendering | Dye through a squared density response, so faint gas stays dark and only dense wisps register. The gas sits behind the dust rather than competing with it |
| Guide field | Optional relaxation of `u` toward `(u·B̂)B̂`, standing in for magnetic tension |
| Diagnostics | Per-cell `\|u\|²`, `\|u\|`, `\|∇·u\|` reduced 4×4 three times, then read back — the HUD numbers are measured, not decorative |

Interaction:

- **Pointer** injects momentum and dye.
- **Scroll** drives a hyperbolic-tangent shear layer with a sinusoidal transverse seed —
  the canonical Kelvin–Helmholtz setup. The roll-up you see as you read is a genuine
  KH instability, not a canned animation.
- **Turning to a chapter** drives the same shear layer in the direction of travel, plus an
  impulse at the point on screen where you pressed — so the page-turn is a physical event
  in the fluid rather than a cut.
- **Each chapter retunes the stopping time.** Small `τ_s` (PIC dust) → grains trace the
  gas. Large `τ_s` (DFMM, cacti) → grains decouple, lag the flow, and concentrate in the
  strain field. Cosmic rays add a field-aligned streaming velocity. Entering `phrike`
  fires a radial blast.

**Pacing.** The simulation clock runs at **one tenth of wall time**, so the field drifts
rather than churns. Every rate-like term (driving, Brownian kicks, reseeding, shear decay) is
scaled by the step, so that choice does not change the physics, the Reynolds number, or the driving/dissipation balance — it only slows the
playback, and it buys a Courant number of a few tenths of a cell per step. Because injected
energy now also lingers ten times longer in wall time, the gradient-subtract pass carries a
smooth governor that asymptotes `|u|` to 420 cells/s; without it, a minute of enthusiastic
mouse movement would pump the field straight back up.

The frame rate is **not** capped at 30 any more. It was, on the argument that the field is
slow enough that a higher rate buys nothing — which held while a grain was smeared into a
trail, because the smear covered the gap between frames. A hard little sprite does not, and a
hundred thousand of them strobe at 30. The gate is now the display rate, and quality is judged
against a **floor** of 34 fps rather than against the target, so a machine holding a steady 40
is left alone instead of being ratcheted down for missing 60.

To retune: `TIME_SCALE`, `TARGET_FPS` and `FPS_FLOOR` sit together near the animation loop in
`js/field.js`; `VMAX` is just above the preset table.

**Cost.** This is a background for a personal site, not a solver anyone will publish from,
and it is tuned that way: 16 Jacobi sweeps (the residual is not visible), a device-pixel ratio
capped at 1.35, and pointer input coalesced to at most one splat per 28 ms — each splat is two
full passes, and a fast mouse can otherwise deliver dozens between frames.

The grains used to be smeared into a persistent trail buffer held at a fraction of the canvas
and then blurred back up on composite: two full-canvas float targets, a full-canvas
read-modify-write every frame, and a five-tap upsample — and it still looked soft. Drawing the
grains where they actually are removed all of that, and the budget it freed bought **3.3× the
grains** (102,400, up from 30,976) and a dye buffer at 1.8× the solver grid instead of 1×, at
equal measured cost per frame. Tiers now shed solver resolution and device pixels before they
shed grains, because the dust is where all the structure is.

Degradation, in order: WebGL2 + float render targets → four quality tiers chosen by measured
frame rate → if the cheapest tier still cannot hold ~14 fps, the solver stops for good and
leaves its last field on screen → a divergence-free canvas-2D curl-noise field where WebGL2 or
float render targets are missing → under `prefers-reduced-motion`, a couple of seconds of
evolution and then freeze. A lost graphics context is caught and the whole GL state is rebuilt.

Render targets are tracked and deleted on every reallocation. They previously were not, and
because mobile browsers fire `resize` on each URL-bar show/hide — i.e. continuously while
scrolling — each scroll leaked a full set of float buffers until the context died.

---

## Where to put pictures

Every image slot below is wired so that **the moment you drop a correctly named file
into `assets/img/`, it appears.** Until then the page shows a dashed card naming the
exact file and aspect ratio it wants. Nothing breaks either way.

### Already filled in (from your own files)

| File | Section | Source |
|---|---|---|
| `bio-portrait.jpg` | 00 Cover | `~/Desktop/profile_pic.jpg` — **worth replacing**, see below |
| `dust-turb-column.jpg` | 01 PIC dust | gas column + plane-of-sky field, ℳ ≈ 23 |
| `dust-ot-coldens.jpg` | 01 PIC dust | dusty Orszag–Tang column density |
| `dust-size-cube.jpg` | 01 PIC dust | grain-size volume render (white background keyed out) |
| `dust-along-b.jpg` | 01 PIC dust | dust/gas ratio and grain size viewed along ⟨B⟩ |
| `dfmm-caustic.jpg` | 02 DFMM | near-caustic `log₁₀(γ₁)` panel |
| `mhd-orszag-tang.jpg` | 03 Non-ideal MHD | UCT-HLLD Orszag–Tang density |
| `cr-bounce.jpg` | 04 Cosmic rays | guiding-centre vs full gyro-orbit in a magnetic bottle |
| `phrike-khi.jpg` | 05 phrike | Kelvin–Helmholtz roll-up |
| `phrike-orszag-tang.jpg` | 05 phrike | Orszag–Tang, N = 1024 |
| `iron-squat.jpg` `iron-walkout.jpg` `iron-deadlift.jpg` `iron-chalk-bw.jpg` `iron-bench.jpg` `iron-podium.jpg` | 08 Iron | Heartbreak Open, from `~/Desktop/Meet photos` |
| `assets/video/iron-reel.mp4` | 08 Iron | your competition reel, re-encoded to 480 px / ~4 MB, plays on hover |

### Still empty — drop these in

| File | Aspect | What it wants |
|---|---|---|
| `bio-portrait.jpg` | 4:5 | The current one is a daylight selfie. A dark-background portrait with hard side light would land far better against this design. |
| `dust-drift-pdf.jpg` | 3:2 | A drift-velocity PDF / survival-function figure from Paper II |
| `dfmm-hero.jpg` | 3:2 | The best DFMM figure you have — a Zel'dovich pancake or a dust-trap field. A colour map, not line plots. |
| `mhd-convergence.jpg` | 3:2 | The Alfvén-wave convergence / ambipolar damping-rate figure |
| `cr-growth.jpg` | 3:2 | The CRSI growth diagnostic — σ_ρ, σ_B, σ_u against Ω_c t |
| `iron-strongman.jpg` | 3:2 | Log, yoke, farmers, stones — anything strongman |
| `cactus-01.jpg` | 4:5 | The hero plant. Dark wall, hard side light, spines rim-lit. |
| `cactus-02.jpg` | 1:1 | Macro of spines or areoles |
| `cactus-03.jpg` | 1:1 | A flower, if you have one |
| `cactus-04.jpg` | 3:2 | The whole garden or shelf, wide |
| `skate-01.jpg` | 3:2 | Best action shot. Mid-air beats rolling. |
| `skate-02.jpg` | 4:5 | The board itself — griptape, worn deck, scuffed rails |
| `skate-03.jpg` | 4:5 | A spot, empty. Concrete as landscape. |
| `skate-04.jpg` | 16:9 | Optional. A video still, or swap the frame for a clip. |
| `contact-wide.jpg` | 3:2 | Optional closing image — a dome, a chalked bar, a desert |

Keep files under ~600 kB. To resize:

```sh
sips -Z 1900 in.jpg --out assets/img/name.jpg
sips -s format jpeg -s formatOptions 72 assets/img/name.jpg --out assets/img/name.jpg
```

To swap the iron video for a different clip (`libx264` is missing from the local
ffmpeg, so this uses hardware VideoToolbox):

```sh
ffmpeg -i clip.mov -an -vf scale=480:-2 -c:v h264_videotoolbox -b:v 620k \
  -pix_fmt yuv420p -movflags +faststart assets/video/iron-reel.mp4
ffmpeg -ss 5 -i clip.mov -frames:v 1 -vf scale=480:-2 assets/video/iron-reel-poster.jpg
```

---

## Where to write text

Every place that needs your words is marked in the page itself with a dashed
**WRITE HERE** block containing `[this is where you talk about …]`. Search the repo for
`this is where you talk about` to find all of them. Delete the `class="ph"` wrapper and
write a normal `<p>`.

Placeholders that are plain `[fill in]` rather than WRITE HERE blocks:

- **09 Iron** — federation, weight class, best squat/bench/deadlift/total, strongman events
- **10 Cacti** — the six genus chips
- **11 Skate** — deck, trucks, wheels, bearings, terrain, and the trick/spot chips

Everything else — publications, education, collaborators, the numbers in the research
sections — is already real, taken from your CV, your papers and your repositories.

---

## Blog

```
blog/
  index.html            all posts
  posts.json            the manifest the lists are built from
  posts/0001-hello.html a post; copy it for new ones
```

To add a post: copy `blog/posts/0001-hello.html` to `blog/posts/0003-your-slug.html`,
write it, then add an entry at the top of the `posts` array in `blog/posts.json`:

```json
{ "file": "0003-your-slug.html", "title": "…", "date": "2026-08-04",
  "summary": "…", "tags": ["dust", "numerics"] }
```

Set `"draft": true` to mark one as unfinished (it still shows — delete the entry to hide
it). The three most recent appear on the home page; the full list lives at `/blog/`.

---

## Structure

```
index.html              the whole site — one document, twelve tabbed chapters
css/main.css            one stylesheet
js/field.js             the WebGL2 dust-gas solver + canvas fallback
js/site.js              the router, reveals, chapter→simulation coupling, HUD, image slots
js/blog.js              renders posts.json
assets/img/             photographs and figures
assets/video/           the lifting reel
assets/favicon.svg
.nojekyll               tells Pages to serve the files as-is
```

Chapters, in order: **00 Cover** (which also carries *Who*) · 01 PIC dust · 02 DFMM dust ·
03 Non-ideal MHD · 04 Cosmic rays · 05 phrike · 06 Publications · 07 Blog · 08 Iron ·
09 Cacti · 10 Skate · 11 Contact.

### The tab deck

The site is **tabbed, not scrolled**. Each chapter is an independent `<section class="panel">`
and exactly one of them is in the document at a time; the other eleven are `display: none`.
That is a layout decision and a performance one — the browser composites one panel's worth of
`backdrop-filter` panes, and the lazy images and the blog feed in the other eleven are never
fetched until you open them.

- The **left rail** (≥ 1080 px) and the **bottom strip** (below it) are two views of the same
  tablist; the strip scrolls the active chapter into its centre on every change. The
  full-screen **Index** overlay is still there as an overview.
- **One history entry per chapter**, so the back button walks the deck and every chapter is
  a real URL (`/#phrike`). Hash links buried in body copy — `#papers`, `#top` — are chapter
  changes too.
- **← / →** turn the page. Up/down and the page keys are left alone, so a long chapter still
  scrolls normally.
- Opening a chapter **replays the reveal cascade**: everything in the opening fold staggers
  in, everything below it waits for the intersection observer, exactly as it did when the
  site was one long scroll.

Each `<section>` carries `data-nav` / `data-label` (its place in the index), `data-preset`
(which simulation tuning to switch to) and `data-accent` (the colour the whole page ramps to).
Presets live at the top of `js/field.js`; add a chapter by adding a `<section class="chapter
panel" role="tabpanel" tabindex="-1" data-nav data-label data-preset data-accent>` and, if you
want new physics, a preset entry. The router picks it up automatically — nothing is hardcoded.

## Local preview

```sh
cd whoisericmoseley-opus
python3 -m http.server 8731
# then open http://127.0.0.1:8731/
```

Serve over HTTP rather than opening `index.html` directly — the blog list uses `fetch`,
which `file://` blocks.

## Accuracy note

The research copy is deliberately conservative and states limits explicitly: the Hall
term is *not* claimed for the non-ideal MHD work, the DFMM methods paper is labelled in
preparation, and the cosmic-ray section says plainly that the first-author CR paper does
not exist yet. If you loosen any of that, loosen it on purpose.

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

No build step, no dependencies, no framework. Static HTML, one stylesheet and plain scripts.
Push to `main` and GitHub Pages redeploys in under a minute.

---

## The three backgrounds

The front page runs the **magnetised** solver, `js/field-mhd.js` — isothermal MHD, two-dimensional
HLLD constrained transport, charged dust — and carries a **Play** chapter that exposes it. Dedner
cleaning remains selectable there. The two schemes it grew out of are kept
as backups: the compressible hydrodynamics at `hydro/` (`js/field-fv.js`) and the incompressible
original at `projection/` (`js/field.js`). `mhd/` and `fv-test/`, where the last two front pages were
developed, are pointer stubs. Every page is the same deck differing only in which engine it loads and
sharing `css/main.css`, `js/site.js`, `assets/` and `blog/`, so they cannot drift; the sub-pages are
`noindex` and carry a banner.

### The Play chapter

`js/play.js` renders a control panel from `field.controls()` and writes back through `field.set()`.
The engine owns the manifest — 28 knobs across gas, solver, field, dust and display — so a knob added
in `js/field-mhd.js` appears on the page with no change to the panel, and the panel cannot disagree
with the engine because it displays the value `set()` returns rather than the one it sent. Two knobs
are not live and say so: the grid reallocates every texture, and the plasma beta re-seeds the field.
Settings persist in `sessionStorage`, never in the URL — the router reads an unrecognised hash as the
Cover, so a settings string there would throw you back to the front page mid-drag.

The Riemann solver is one of those knobs, which is the most interesting thing to play with: **LLF**
(one wave), **HLL** (two), **HLLC** (three, contact resolved), **HLLD** (five, Alfvén waves resolved).
All four cost the same ~10 ms per frame, measured, so the switch is free; what changes is what
survives. HLL has no contact, so it averages the transverse velocity and the passive scalar across
every interface — and a shear layer is a discontinuity in exactly those and nothing else.

### The display path

Two things live between the grid and the screen and neither is part of the solver.

**Fronts.** The convergence channel is a centred difference, so a captured front is one cell wide and
effectively all of its power sits at the grid's Nyquist frequency. Catmull-Rom is right for the dye —
it interpolates, so a wisp keeps its amplitude — but its kernel has negative lobes, and negative lobes
on a one-cell feature ring: an undershoot ridge either side, axis-aligned because the kernel is
separable. Stretched over five to eight screen pixels per cell that reads as grid facets on every
shock. The fix is per channel: the dye keeps Catmull-Rom, the convergence gets a **cubic B-spline**,
which is C2 and strictly positive so it cannot ring at all, for four bilinear taps against nine. Plus
a deliberately small grid-level blur (`SHOCK_BLUR`, 0.4 cells) so the front is something an
interpolant has any information about, and a soft shoulder (`1 - exp(-x)`) in place of the clamp,
whose saturation plateau had a grid-following contour of its own. Measured: frame cost unchanged,
frame brightness unchanged.

**The pointer.** A splat is a velocity increment, not a force, and site.js emits one every 28 ms, so
the per-event ceiling never bounded a stroke. Thirty-six a second at the old ceiling was 86 sound
speeds a second into a spot a twentieth of the box across: peak Mach 8 on the front page and 22 on the
magnetised one, ctot past the point where the Courant condition takes over from the default timestep,
and still throttled six seconds after the mouse stopped. Now the momentum draws on a leaky budget and
the ink does not — the ink is a passive scalar that costs the solver nothing and is what makes a stroke
visible.

### The front page — isothermal MHD, plasma beta 5 by default

Magnetised gas, charged dust. Measured at `gridH = 96`, 172 × 96:

| | |
|---|---|
| flux | HLLD (Miyoshi & Kusano), the reference's `hlld_mhd_fluxes` reduced to isothermal and 2D |
| reconstruction | PLM (`slope_type = 2` MonCen), matching mini-RAMSES's CT trace. PPM (CW84) remains selectable with the Dedner branch |
| div B | staggered constrained transport driven by mini-RAMSES's `riemann2d_hlld` corner EMF; Dedner GLM remains selectable |
| targets | CT stores `Bx` and `By` on their left and bottom faces; Dedner stores cell-centred `(Bx, By, psi)` |
| dust | Lorentz on the drift by mini-ramses's Cayley rotation, then **one** backward-Euler drag stage; charge-to-mass 100; three velocity components |
| pace | 0.055 box heights per wall second, against the live site's 0.057 |
| field, drawn | line integral convolution along **B**, signed and accent-tinted, normalised to unit variance; **B** toggles it |
| fronts, not drawn | the convergence overlay is off here — long-lived at beta = 1, and a second overlay competing with the field. **S** puts it back; its prepare pass does not run while it is off |
| cost | ~10 ms/frame, 1.7x headroom at 60 Hz; boot 0.7 s; mass 1.00003 over 100 s |

Two things came out of building it that are worth keeping in mind.

**Two dimensions have no energy sink.** Enstrophy cascades to the grid but energy cascades the other
way, into box-scale structures that nothing dissipates. Measured: with the drive cut by a factor of a
hundred the rms Mach held at 0.85 for half a minute while the peaks decayed. The servo had no plant
to act on and wound itself to the floor. `FRIC` is the large-scale drag that gives the box a real
equilibrium, which is what a 2D driven-turbulence run does for the same reason, and the servo is PD
rather than P because the magnetised plant answers slowly.

**PPM has one MHD-specific trap.** The PPM face state is the parabola averaged over the slab the
flow crosses in half a step, which carries the `-u dq/dx` term for you at one speed for every
variable. That is fine for the hydrodynamic variables, which are all advected along the sweep normal.
It is wrong for two of the three field variables: `Bn` has no such term at all — the normal component
is advected only transversely, which is the fact constrained transport is built around — and `psi` is
not advected by anything. Shifting them anyway injects a spurious `-u dBn/dx` of exactly the order of
the terms being added, straight into the divergence the cleaning is damping. So the field takes the
slab average in its transverse component and the plain parabola edge in the other two.

**A scroll does not roll up here.** The mean field is vertical, the layer a scroll deposits is
horizontal, and field-line tension holds a layer against Kelvin–Helmholtz until `delta_u > 2 vA` —
2.8 sound speeds at beta = 1, against the 1.4 one gesture deposits. What you get instead is
Alfven waves and current sheets. `BETA` is the knob if you want the vortices back.

`field-fv.js` is a drop-in `window.Field`: same public surface, same dust, same palette, same
per-chapter presets, entirely different gas. Isothermal Euler, piecewise constant, HLL with RAMSES's
clamped wave speeds, unsplit conservative update, `cmpdt` per step reduced on the GPU. The dye rides
along as a **passive scalar in the fourth conserved slot** — which is exactly what `scalar_flux`
reduces to in the HLL branch, componentwise HLL on `rho*Y` with flux `u*(rho*Y)`, so it is free.

**What it shows that the incompressible build cannot:** real shock fronts, and dust piling up
against them. Preferential concentration in a divergence-free field happens in strain; here it
happens at discontinuities the Riemann solver captured.

**What it costs, which is the point of running the experiment.** An explicit compressible scheme is
CFL-limited by the *acoustic* speed, so it needs roughly `ctot/(C cs) * gridH ~ 20 * gridH` steps per
box crossing time. A projection method needs none of them — its pressure is elliptic, so it steps on
the advective timescale alone. Measured: at `gridH = 176` this background reached **t = 0.33 crossing
times in 25 s of wall clock at 22 fps**, which is not a background, it is a still image with a frame
counter. Getting to 60 fps and `t = 1.28` in 22 s took three concessions:

1. **Coarser gas** — `gridH` 112 down from 192. dt improves linearly, pass cost quadratically, and
   the dust is drawn at native resolution regardless, so the visible loss is fatter shocks.
2. **One Courant reduction per step**, ahead of the source terms, rather than one before and one
   after. Three passes saved for a correction the drive servo keeps small.
3. **A seeded initial condition** — a few large-scale modes and a smooth dye pattern, because
   spinning up from rest costs thousands of steps and the page would open on a still frame.

Also lowered the target Mach to ~2.75: at Mach 5 the rarefactions span five decades and pin the
density floor over a sixth of the box.

## Cache busting

Asset URLs in `index.html` and `shocks.html` carry a `?v=<token>` version token. GitHub Pages serves
`css/` and `js/` with `Cache-Control: max-age=600` and no fingerprinting, so a browser holding the
previous stylesheet alongside new HTML renders a layout that never existed — old CSS showing the left
rail while the new router runs the deck, for one real example, which looks exactly like a deploy that
did not take.

**Bump the token whenever anything under `css/` or `js/` changes.** There is no build step to do it
automatically.

## One layout

The chapters are a deck at every width: exactly one panel is in the document, the tab strip along
the bottom is the navigation, and arrow keys turn the page. On a phone the strip scrolls
horizontally to keep the active tab visible. The retired upper-corner index and its sheet stay
dormant; the bottom strip is the visible, accessible tablist on both pointer and touch devices.

## The compressible sub-page

`shocks.html` is the other half of the subject. The front-page background is incompressible by
construction — a projection method with an elliptic pressure has no shocks in it anywhere. That page
solves the **isothermal Euler equations** conservatively instead:

| | |
|---|---|
| State | `(rho, rho u, rho v)`, periodic, square-celled. Isothermal closure `P = rho cs^2` with `cs = 1`, so every velocity on the page *is* a Mach number |
| Reconstruct | **none** — piecewise constant, so the interface states are the cell averages. First-order Godunov: no slopes, no limiter, no predictor half-step |
| Flux | HLL with RAMSES's wave-speed estimate, `SL = min(min(uL,uR) - cs, 0)`, `SR = max(max(uL,uR) + cs, 0)`. Clamping through zero means `SL <= 0 <= SR` always, so the central formula is always correct, the supersonic branches vanish, and `SR - SL >= 2cs` needs no divide guard. **LLF** (Rusanov) selectable: one speed, cheaper, more diffusive |
| Update | unsplit and conservative, `U += dt/dx [(Fx_i - Fx_i+1) + (Fy_j - Fy_j+1)]`, all four faces in one pass on a five-point stencil. Each face is solved twice, once per adjoining cell — wasteful on a CPU, right here, because a second pass to store the flux costs more than re-solving |
| Timestep | RAMSES `cmpdt` verbatim: `ctot = sum_dim(|u_dim| + cs)`, `dt = C dx / ctot`, ceiling `C dx / smallc`. `C = 0.8` (the RAMSES default) with **no** safety margin. Recomputed from the current state **before every step**, the way `courant_fine` runs ahead of `amr_step`; the max is reduced on the GPU to one texel and read by the solver as a texture, so there is no readback in the loop |
| Precision | single throughout; WebGL2 `highp float` is IEEE binary32 and there is no double path. No `-ffast-math` switch exists in GLSL ES, so the intent is in the code: reciprocal multiplies, **no sqrt at all** (an isothermal sound speed is a uniform), branch-free Riemann solver |
| Provenance | the scheme from `mini-ramses-ism`'s GPU hydro path (`gpu_hydro.cuf`) — same Riemann solver, same wave speeds, same floors in the same places, same Courant condition — reduced to 2D and to its cheapest configuration. Missing relative to the reference: AMR, MHD, constrained transport, and the slope/trace machinery piecewise-constant reconstruction makes unnecessary |

**What it demonstrably gets right:** total mass drifts by ~5e-8 relative over twenty thousand steps
— float32 round-off and nothing else, which is the check that the flux differencing is genuinely
conservative. And the density PDF goes log-normal on its own from a uniform initial condition.

### The instability, and two wrong fixes

Symptom: after minutes, an instability swept the box in ~1 s and left it uniform.

Cause: **a stale CFL**. `dt` was sized from a max signal speed read back every 16th step. At high
Mach, deep rarefactions put ~1/6 of cells on the density floor, and a floor on `rho` is not a floor
on `u = m/rho` — so the true signal speed could multiply several-fold inside one measurement
interval. An explicit scheme past its CFL limit spreads error ~1 cell/step; across 460 cells that is
about the second it took. The uniform box afterwards was my own guard calling the initial
conditions, which destroyed the evidence each time.

Two fixes that were wrong, recorded so they don't get tried again:

1. **A safety margin on `dt`** — compensating for staleness instead of removing it.
2. **Scaling the margin by the observed growth in signal speed** — *worse than the bug*. Ordinary
   fluctuation reads as growth, so `dt` was throttled nearly permanently and the turbulence never
   developed: sigma collapsed from ~1.4 to ~0.1. A dead box is not a stable one.

The fix was to do what RAMSES does and compute the Courant condition from the current state every
step. Same abuse (Mach 12 + sustained detonations and fast strokes): peak `|u|` 79 → 30, peak `ctot`
105 → 37, no velocity limiter anywhere, no restarts. Costs 3 extra passes per step (~65 fps → ~45).

**What it does not:** the value of `b` in `sigma^2 = ln(1 + b^2 M^2)`. Measured sigma from this box
scatters by tens of per cent — it is a second moment of a heavy-tailed field over about three box
crossing times — and the relation is calibrated in 3D anyway. The page says so.

Two bugs worth remembering from building it, both of which produced a *silent* wrong answer rather
than an error:

- Setting `LINEAR` filtering on an `RGBA32F` texture makes it **incomplete** unless
  `OES_texture_float_linear` is present, and then every read returns zero — `texelFetch` included.
  A display-only change killed the entire solver. Interpolation is now written out by hand, and
  `reset()` reads one pixel back and throws if the state is not there.
- The drive servo was multiplicative *per step*. At ~500 steps/s a 2%-per-step gain compounds to
  x148 per second, saturates the amplitude between two measurements, and violates the CFL condition.
  It now runs once per measurement with a hard slew limit.

## The background is a real simulation

`js/field.js` is a two-dimensional incompressible Navier–Stokes solver coupled to a
population of Lagrangian dust superparticles, written from scratch in WebGL2. Per frame:

| Stage | Method |
|---|---|
| Gas advection | Semi-Lagrangian backtrace, bilinear, periodic |
| Small-scale swirl | Fedkiw vorticity confinement, `f = ε ω (N_y, −N_x)`, `N = ∇\|ω\|/\|∇\|ω\|\|` |
| Incompressibility | Divergence → 16-iteration Jacobi pressure solve → gradient subtract, with a smooth speed governor at `\|u\| = 420` cells/s |
| Dust | `dv/dt = (u_gas − v)/τ_s`, advanced by **first-order implicit (backward Euler)**: `v ← (v + a·u)/(1+a)`, `a = dt/τ`. The population is **polydisperse** — each grain draws τ from a log-uniform spectrum 1.4 decades wide, hashed out of its own texel rather than stored — so the median grain's `dt/τ` is the only scalar the CPU supplies and a grain costs one divide, one multiply-add and one multiply, still L-stable |
| Grain recycling | **Off.** No grain is ever deleted and reborn, so the clustering on screen is the converged attractor of the drag law, not an average over a population that keeps being reset — preferential concentration takes many eddy turnovers to build and recycling truncates exactly that tail. `RESEED_RATE` in `js/field.js` restores it; `0.06` is the old value |
| Dust rendering | One additive soft sprite per grain, at native resolution, straight to the screen |
| Grain brightness | The **drift** `\|v − u_gas\|`, and nothing else — the only quantity in the drag law that does work on a grain. Tightly coupled dust stays dim; decoupled grains light up where they pile up. A large and a small grain drifting at the same speed are equally bright per unit area |
| Grain apparent size | Literally proportional to grain radius. In the Epstein regime `tau ∝ a`, so the τ multiplier *is* the size multiplier, and because τ is drawn log-uniformly the population is already sampled evenly per log interval in size. Normalised so the largest grain matches the previous largest apparent size, which makes the range 25:1 across the 1.4 decades of τ. The bottom of that range is sub-pixel, so the smallest grains sit at a 1 px floor and are slightly over-represented in area — a limit of point rasterisation, not of the distribution. Exact *within* a chapter; between chapters `pointSize` tracks the median τ only weakly, or the cosmic-ray grains would be a fiftieth the size of the cactus ones and invisible |
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
rather than churns. Every rate-like term (driving, Brownian kicks, shear decay) is
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

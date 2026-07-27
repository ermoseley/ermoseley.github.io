# whoisericmoseley

Personal site for **Eric R. Moseley** — numerical astrophysicist, KIPAC Fellow at Stanford/SLAC.

**Live:** https://ermoseley.github.io/whoisericmoseley-opus/

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
| Incompressibility | Divergence → 30-iteration Jacobi pressure solve → gradient subtract |
| Dust | `dv/dt = (u_gas − v)/τ_s`, advanced with the **exact exponential integrator** `v ← u + (v−u)e^(−dt/τ)`, so stiff tightly coupled grains never go unstable |
| Dust rendering | Grains splatted into a fading trail buffer, so the streaks are real particle paths |
| Guide field | Optional relaxation of `u` toward `(u·B̂)B̂`, standing in for magnetic tension |
| Diagnostics | Per-cell `\|u\|²`, `\|u\|`, `\|∇·u\|` reduced 4×4 three times, then read back — the HUD numbers are measured, not decorative |

Interaction:

- **Pointer** injects momentum and dye.
- **Scroll** drives a hyperbolic-tangent shear layer with a sinusoidal transverse seed —
  the canonical Kelvin–Helmholtz setup. The roll-up you see as you read is a genuine
  KH instability, not a canned animation.
- **Each section retunes the stopping time.** Small `τ_s` (PIC dust) → grains trace the
  gas. Large `τ_s` (DFMM, cacti) → grains decouple, lag the flow, and concentrate in the
  strain field. Cosmic rays add a field-aligned streaming velocity. Entering `phrike`
  fires a radial blast.

Degradation, in order: WebGL2 + float render targets → four quality tiers chosen by
measured frame rate → a divergence-free canvas-2D curl-noise field → under
`prefers-reduced-motion`, three seconds of evolution and then freeze.

---

## Where to put pictures

Every image slot below is wired so that **the moment you drop a correctly named file
into `assets/img/`, it appears.** Until then the page shows a dashed card naming the
exact file and aspect ratio it wants. Nothing breaks either way.

### Already filled in (from your own files)

| File | Section | Source |
|---|---|---|
| `bio-portrait.jpg` | 01 Who | `~/Desktop/profile_pic.jpg` — **worth replacing**, see below |
| `dust-turb-column.jpg` | 02 PIC dust | gas column + plane-of-sky field, ℳ ≈ 23 |
| `dust-ot-coldens.jpg` | 02 PIC dust | dusty Orszag–Tang column density |
| `dust-size-cube.jpg` | 02 PIC dust | grain-size volume render (white background keyed out) |
| `dust-along-b.jpg` | 02 PIC dust | dust/gas ratio and grain size viewed along ⟨B⟩ |
| `dfmm-caustic.jpg` | 03 DFMM | near-caustic `log₁₀(γ₁)` panel |
| `mhd-orszag-tang.jpg` | 04 Non-ideal MHD | UCT-HLLD Orszag–Tang density |
| `cr-bounce.jpg` | 05 Cosmic rays | guiding-centre vs full gyro-orbit in a magnetic bottle |
| `phrike-khi.jpg` | 06 phrike | Kelvin–Helmholtz roll-up |
| `phrike-orszag-tang.jpg` | 06 phrike | Orszag–Tang, N = 1024 |
| `iron-squat.jpg` `iron-walkout.jpg` `iron-deadlift.jpg` `iron-chalk-bw.jpg` `iron-bench.jpg` `iron-podium.jpg` | 09 Iron | Heartbreak Open, from `~/Desktop/Meet photos` |
| `assets/video/iron-reel.mp4` | 09 Iron | your competition reel, re-encoded to 480 px / ~4 MB, plays on hover |

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
index.html              the whole single-page site
css/main.css            one stylesheet
js/field.js             the WebGL2 dust-gas solver + canvas fallback
js/site.js              nav, reveals, section→simulation coupling, HUD, image slots
js/blog.js              renders posts.json
assets/img/             photographs and figures
assets/video/           the lifting reel
assets/favicon.svg
.nojekyll               tells Pages to serve the files as-is
```

Sections, in order: hero · 01 Who · 02 PIC dust · 03 DFMM dust · 04 Non-ideal MHD ·
05 Cosmic rays · 06 phrike · 07 Publications · 08 Blog · 09 Iron · 10 Cacti ·
11 Skate · 12 Contact.

Each `<section>` carries `data-preset` (which simulation tuning to switch to) and
`data-accent` (the colour the whole page ramps to). Presets live at the top of
`js/field.js`; add a section by adding a `<section data-nav data-label data-preset
data-accent>` and, if you want new physics, a preset entry.

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

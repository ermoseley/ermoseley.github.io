# Constrained transport, parked

Not wired into the engine. Nothing under `js/` references anything in this directory,
and the Play chapter has no control for it. This is work in progress kept where it can
be picked up, not a feature that is switched off.

Both pieces here are verified. What is missing is not correctness of the parts but the
scheme that would use them: face-centred storage, the predictor pass, the corner
reconstruction, and the update itself. See "What is still missing" at the end — that is
the remaining work, and it touches every consumer of `B` outside the solver.

The front page cleans `div B` with Dedner's generalised Lagrange multiplier: `B` is
cell-centred, `psi` is an eighth variable that carries divergence error out of the box
at speed `ch` and is damped as it goes. Constrained transport is the alternative — a
staggered field and a corner EMF — and the reason to want it is that CT keeps `div B`
at zero identically rather than damping it toward zero. The reference
(mini-ramses, `codex/stage0-9-yldh-production`) runs `riemann='hlld'`,
`riemann2d='hlld'`, `slope_mag_type=2`.

## What is here, and what state it is in

### `riemann2d-hlld.glsl` — verified, and ready to be wired

A GLSL port of `riemann2d_hlld` (Miyoshi & Kusano) from
`godunov_utils.f90:1995-2212`: the 2-D Riemann problem at a cell corner, returning
`Ez = u By - v Bx`, which is the code convention and minus the physical `Ez`. Simplified
for two dimensions only where that is exact — `Bz` and `vz` are homogeneous here, so if
they start at zero they stay there, and the `C` terms vanish. All six branches of the
EMF selection and both Alfvén-speed pairs are kept.

`test-riemann2d.mjs` compares a mechanical JavaScript mirror of this GLSL against
`ref-riemann2d.mjs`, an independent transcription of the Fortran, over 20 000 randomised
quadrant states. All eleven checks pass:

```
ok  mirror is token-identical to the GLSL       ok  leaf selection identical to the Fortran
ok  transposition antisymmetry to 1e-10         ok  every leaf exercised >= 200 times
ok  well-conditioned agreement < 1 float32 eps  ok  near-degenerate agreement < 1 eps
ok  contact floor unreachable                   ok  rc floor unreachable
ok  never non-finite on hostile states
```

Three of those are worth naming, because they are what makes the rest mean anything. The
mirror is *token-identical* to the GLSL rather than a second reading of the Fortran, so
the comparison cannot pass by two transcriptions sharing a mistake. Every leaf is
exercised at least 200 times, so branch coverage is measured and not hoped for. And the
guard epsilons are shown unreachable on the sampled states, so they cannot be quietly
changing a well-conditioned answer.

It has also been run on real hardware, which for parked shader code is the doubt that
matters most: `gpu-riemann2d.mjs` compiles and links this GLSL in WebGL2 against the
current EOS block and pushes 8192 corners through it, in both EOS modes.

```
iso   max rel 1.048e-6   (8.8 x float32 eps)
adi   max rel 9.982e-7   (8.4 x float32 eps)
```

Mutation testing, 13 seeded defects: 10 killed. The three survivors are the two floors
the test proves unreachable and a redundant `sound2` clamp, each now asserted as a
measured property rather than left as a silent gap.

Two caveats, both inherited rather than introduced:

- `riemann2d_hlld` is genuinely discontinuous across its branch boundaries — measured,
  it jumps by up to 190% of |E| across `ST = 0`. The port matches that bit-for-bit at
  equal precision, but in float32 a corner within ~1e-7 of a boundary is a coin flip and
  E can differ by O(|E|). No port can fix that.
- Where the star-density floor binds, the port stops tracking the Fortran. That happened
  only below the engine's own density floor (2822 of 32800 samples, all from a
  deliberately unphysical eight-decade sweep; none inside physical ranges). Consistency
  with `riemann()`'s own floor was judged worth more than following the reference into a
  density the rest of the file has already clamped.

One deliberate non-optimisation in the GLSL: `cfast2` keeps two divides by `r` rather
than one reciprocal. Those eight numbers set SL/SR/SB/ST, and the branch selection is
discontinuous where one crosses zero — with `* ri` the port sat one ulp from the Fortran
and picked the other leaf on boundary states; with `/ r` it is bit-identical.

Run it:

```sh
cd notes/constrained-transport && node test-riemann2d.mjs
node gpu-riemann2d.mjs          # needs Chrome; drives it over CDP
```

### `divb-project.mjs` — verified, and the piece that makes the switch safe

CT *preserves* whatever divergence it starts from. Seeding it by interpolating the
cell-centred Dedner field onto faces would therefore freeze that field's 0.6% error in
place permanently, which is worse than the cleaning it replaced.

The way out is the vector potential. In 2-D any in-plane divergence-free field is
`B = curl(A z)`, and with `A` at cell corners

```
Bx(i-1/2, j) = [A(i,j+1) - A(i,j)] / dy
By(i, j-1/2) = -[A(i+1,j) - A(i,j)] / dx
```

the discrete divergence telescopes the same four corner values with opposite signs, so
it vanishes to round-off *whatever `A` is*. Convergence therefore buys fidelity, never
correctness, which is what makes the conversion safe to do mid-run. The same reason
`vecpotentialinit.f90` exists in the reference.

The remaining question is only how closely the rebuilt field resembles the one on
screen, i.e. solving `grad^2 A = -curl B`. This does it exactly instead of iteratively:
the box is periodic and the composite map from `A` to the cell-centred field is diagonal
in Fourier space, so the least-squares solution is one division per mode. `k = 0` is the
uniform field, which is not representable by a periodic `A` and is carried separately as
a mean — divergence-free on its own. The single (Nyquist, Nyquist) checkerboard mode is
also in the null space and is genuinely not the face average of any divergence-free
field, so it is dropped rather than fudged.

Measured (`divb-project.test.mjs`), at 152x96:

```
naive face interpolation   normalised div B  8.0e-2
after projection           normalised div B  4.2e-17     <- round-off
exactly-representable field round-trips to   1.5e-13
solve time                                   19 ms, one-time
```

19 ms is a one-time cost behind a toggle, so an exact direct solve was worth more than
avoiding the pipeline stall of a `readPixels`. `js/field.js` has a tuned multigrid
V-cycle that could do this on the GPU without the stall, but it is shaped as a smoother
— 0.118 relative residual at V(2,1)x3 — which is right for a pressure *correction* and
not for solving `A` itself.

Run it:

```sh
node notes/constrained-transport/divb-project.test.mjs
```

## What is still missing

Beyond fixing the corner solver:

- **Face-centred storage.** In CT mode the field target's `.xy` become `Bx` on the left
  face and `By` on the bottom face. Cell-centred `B` is then the average of a cell's own
  two faces, exactly as `ctoprim` does it, and every consumer outside the solver needs
  it: the flow texture the dust samples, the display pass, the Lorentz force, the
  divergence diagnostic, the metrics reduction.
- **Two passes, not one.** A cheap predictor EMF at each corner
  (`Ez = u By - v Bx` from the four-cell average, which `trace2d` uses to time-centre
  the face fields by half a step), then the real corner solve. Computing the predictor
  inline would need a 4x4 cell stencil per corner instead of a 3x3 texel read.
- **The corner states.** `trace2d`'s bilinear reconstruction, with the crucial asymmetry
  that the component *normal* to a face is single-valued and taken from the staggered
  array with no normal slope, while transverse components are reconstructed like hydro
  variables. `slope_mag_type=2` is the same MonCen limiter as `slope_type=2` applied to
  the face-averaged field, in transverse directions only, with the multiplier clamped by
  `MIN(slope_mag_type, 2)`.
- **The update.** `Bx += (dt/dy)(Ez(i,j+1) - Ez(i,j))`,
  `By -= (dt/dx)(Ez(i+1,j) - Ez(i,j))`, and with CT owning the normal component the
  `psi` flux and the Powell source term both go away — a CT scheme needs neither.

The reference measured its `uct-hlld` variant, which builds the edge EMF from the 1-D
HLLD fans instead of solving a second Riemann problem, at 8% *slower* than hlld/hlld,
so there is no cheaper route worth taking here.

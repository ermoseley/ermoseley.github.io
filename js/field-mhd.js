/* ============================================================================
   field-mhd.js — the same background, magnetised.

   A third gas for the same page. field.js is incompressible (projection);
   field-fv.js is compressible hydrodynamics (the live site); this is compressible
   *magneto*hydrodynamics at plasma beta = 1, with charged dust.

   Everything above the gas -- the palette, the chapters, the dye, the grain size
   spectrum, the interaction surface site.js drives -- is unchanged. What differs
   is what the grains are falling through, and that they now feel it twice: once
   through drag, once through the Lorentz force.

   ---------------------------------------------------------------- the gas

     variables   : (rho, rho u, rho v, rho Y) and (Bx, By, psi), two RGBA32F
                   targets written in one pass through MRT. Isothermal,
                   P = rho cs^2, cs = 1, so every velocity is a Mach number and
                   every field strength is in units of rho^1/2 cs.
     div B       : Dedner et al. (2002) generalised Lagrange multiplier. B is
                   cell-centred and psi is an eighth variable that carries
                   divergence error out of the box at speed ch and is damped as it
                   goes. This is the alternative to the reference's constrained
                   transport, and it is the right one here: CT needs a staggered
                   field and a corner EMF solver, which is two more passes and a
                   second stencil, where GLM is two extra numbers in a texture that
                   the same 1D solver handles as part of its own flux.
     flux        : HLLD (Miyoshi & Kusano), the reference's hlld_mhd_fluxes reduced
                   to isothermal and to two dimensions -- five waves, so the
                   Alfven waves and the contact are all resolved and a shear layer
                   in the transverse field survives.
     reconstruct : PLM. slope_type = 2, the MonCen limiter from uslope, on all
                   seven primitives, plus a half-step predictor in the primitive
                   variables. No PPM on this page: a parabola fitted through a
                   fresh discontinuity is monotonised variable by variable, which
                   is fragile enough in hydrodynamics and worse when two of the
                   variables are a magnetic field whose divergence is supposed to
                   vanish.
     update      : unsplit and conservative, all four faces in one pass, plus the
                   one non-conservative term GLM needs (below).
     timestep    : dt = min(dt_default, dt_Courant), as on the live site. ctot is
                   now the sum over directions of |u_dim| + cf_dim with cf the fast
                   magnetosonic speed, so it is roughly twice the hydrodynamic
                   value at the same Mach number and the default step is set lower
                   to keep the fixed branch the binding one.

   Three places where this departs from the reference, all deliberate:

     1. Isothermal HLLD. hlld_mhd_fluxes is adiabatic and carries an energy
        equation; there is none here. P = rho cs^2 goes in wherever gamma*P/rho
        appeared, the etot star states and the v.B terms drop out entirely, and
        what is left is the same wave structure and the same algebra. The star
        region keeps a density on each side of the contact, exactly as the
        reference's does. Strictly, isothermal MHD has no entropy wave and the two
        should be equal; the same approximation is in the isothermal HLLC the live
        site runs, where it costs nothing visible and buys a resolved tangential
        discontinuity.

     2. GLM rather than CT, as above. The normal-field flux is not zero here
        (hlld_mhd_fluxes sets flux%Bx = 0 because CT owns that component): it is
        psi, and psi's own flux is ch^2 Bn. That pair is a linear 2x2 hyperbolic
        system with eigenvalues +-ch, decoupled from the rest, so it is solved
        exactly and added to the HLLD flux rather than approximated by it.

     3. The Powell source term -(div B) B in the momentum equation, which GLM needs
        and a CT scheme does not. Without it a residual divergence exerts a force
        along B that no physics put there. It is evaluated from the same
        face-averaged normal fields the psi flux uses, so the divergence the
        momentum equation feels is the one the cleaning is damping.

   -------------------------------------------------------- the thermal state

   The gas above is isothermal and that is still the default, but it is now the
   default of a switch rather than the only thing the file can do. adia = 1 carries
   a total energy and a gamma law -- one more channel in the field target, one more
   flux out of the same HLLD, one more primitive through the same limiter -- and on
   top of that a second operator: the atomic heating and cooling of the interstellar
   medium, split after the hydrodynamic step exactly where amr_step puts it.

   The two go together and are offered as one switch, because a gas that cannot
   change its temperature cannot have two phases and two phases are the whole point.
   The inventory is mini-ramses's ism_net_cooling_h2 with the molecular hydrogen set
   to zero: grain photoelectric heating with the charge-efficiency fit, cosmic rays
   on atomic H, C+ and O fine structure, Lyman-alpha, grain-assisted recombination,
   and a blended CIE metal term. Lambda changes sign twice along an isobar, so the
   equilibrium pressure turns over twice and the branch between the turning points
   is isobarically unstable -- and that is not a feature that had to be added. The
   cold phase is the instability doing its own work.

   What had to be arranged is that switching it on changes nothing. The pressure of
   the isothermal page is rho cs^2 with cs = 1; the equilibrium pressure of the
   cooling function at its unstable branch's centre is a number in erg/cm^3. Putting
   rho_code = 1 at that centre and taking the velocity unit to be that state's own
   sound speed makes the two the same number identically, not approximately -- see
   UNITS. So the switch changes how the gas answers being compressed and nothing
   else, and the two phases that follow (49 K at 15.6 rho_code, 7800 K at 0.096) are
   a consequence of it rather than of a second set of choices.

   The one free unit left over is the length, and it is a control: it sets how many
   cooling times fit into a box crossing, which is the only thing that decides
   whether the cold phase is resolved or is a picture of the mesh.

   ---------------------------------------------------------------- the dust

   Charged grains, charge-to-mass ratio 100 for the whole population. The update is
   mini-ramses's, from pm/move_fine.f90:

     - the force acts on the *drift* velocity w = v - u_gas, which is correct and
       not a shortcut: in ideal MHD the electric field in the lab frame is -u x B,
       so q(E + v x B) = q(w x B). A grain moving with the gas feels nothing.
     - compute_lorentz, the implicit Cayley form
         w' = Inverse[e - (dteff/2) m] . [e + (dteff/2) m] . w,   dteff = -dt*charge
       with m the cross-product matrix of B. That is an exact rotation about B for
       any step -- it preserves |w| to roundoff and cannot heat the population --
       and it is the reference's default (analytic_dust_force = .false.). The
       matrix is theirs, with Bz = 0 substituted; the sign of dteff is theirs too.
     - then drag. Here is the one exception the brief asked for: a single backward
       Euler stage after the Lorentz update, w <- w/(1 + dt/tau), rather than the
       reference's half-drag / Lorentz / half-drag splitting. One stage, first
       order, L-stable at any stopping time.

   The grains are 2.5-dimensional and the gas is not. B lies in the plane, so
   w x B for an in-plane drift points out of it: a charged grain cannot stay in the
   plane even though the flow it lives in is planar. So the population carries
   (vx, vy, vz) -- a second particle target -- and gyrates in the plane containing
   the field normal. The gas needs no such thing: with Bz = vz = 0 initially, their
   equations are homogeneous and they stay zero exactly.

   At charge 100 and |B| ~ 1.4 the gyrofrequency is 140 and the gyration is
   resolved at about thirty steps a turn (the reference's gyro_factor limits the
   same angle to 0.63 rad; this runs at 0.2). The gyroradius is a twelfth of a
   cell, so the grains are strongly magnetised: they slide along field lines and
   are held across them, and the dust draws the field rather than the flow.

   That anisotropy is also why the Kelvin-Helmholtz layer a scroll deposits behaves
   differently here. The mean field is vertical and the layer is horizontal, so the
   field threads it at right angles and its tension resists the roll-up -- the
   classical result is that a layer is stable to KH once delta_u < 2 v_A, and at
   beta = 1 that is delta_u < 2.8 cs. A scroll here makes Alfven waves and current
   sheets where the live site makes vortices. That is the physics of beta = 1 with a
   mean field, not a bug; the plasma beta is a control on the page.

   Written from scratch. No libraries.
   ========================================================================= */

(function (global) {
  'use strict';

  // ---------------------------------------------------------------- shaders

  const VERT = `#version 300 es
precision highp float;
in vec2 aPos;
out vec2 vUv;
void main(){ vUv = aPos * 0.5 + 0.5; gl_Position = vec4(aPos, 0.0, 1.0); }`;

  const HEAD = `#version 300 es
precision highp float;
precision highp sampler2D;
in vec2 vUv;
out vec4 outColor;
`;

  // Two attachments: the hydrodynamic state and the field. Explicit locations,
  // because ES 3.0 has no bindFragDataLocation.
  const HEAD2 = `#version 300 es
precision highp float;
precision highp sampler2D;
in vec2 vUv;
layout(location = 0) out vec4 outColor;
layout(location = 1) out vec4 outColor1;
`;

  // The cooling-off solver is compiled from a preprocessor branch containing the
  // shipped arithmetic verbatim. A uniform branch is numerically equivalent but makes
  // the driver choose different fused operations in the Godunov kernels, enough for a
  // one-ulp drift that grows under a chaotic flow.
  function isothermalShader(src) {
    return src.replace('#version 300 es\n', '#version 300 es\n#define ISOTHERMAL_ONLY 1\n');
  }

  // The primitive state, in two vec4s so one limiter serves both. The field slot
  // carries psi as its third component -- limited and predicted like any other
  // primitive -- and the pressure as its fourth.
  //
  // Pressure is in the state rather than derived at each use because two different
  // equations of state now feed it, and because the constrained-transport corner
  // states need it reconstructed alongside everything else. Under the isothermal law
  // it is redundant, being exactly rho cs^2; under a gamma law it is the primitive
  // form of the conserved energy the field target's fourth channel carries. Same slot
  // in both, which is what lets one trace and one Riemann solver serve both.
  const EOS = `
uniform float cs, cs2, smallr;

// The equation of state, on a switch. adia = 0 is the isothermal gas this page has
// always run: P = rho cs^2, cs = 1, so every velocity is already a Mach number and
// every field strength is in units of rho^1/2 cs. adia = 1 carries a total energy and
// a gamma law, which is what cooling requires -- a gas that cannot change its
// temperature cannot have two phases, and the whole point of the cooling option is
// that it does.
uniform float adia, gam;

struct Prim { vec4 h; vec4 f; };      // h = (rho,u,v,Y), f = (Bx,By,psi,P)
struct Flux { vec4 h; float bn; float bt; float ps; float en; };

vec4 toCons(vec4 q){ return vec4(q.x, q.x * q.y, q.x * q.z, q.x * q.w); }

// The sound speed squared, from whichever law is in force. Every wave speed in the
// file goes through this one function: find_speed_fast, the Courant reduction and the
// corner solver all need it, and three copies of an EOS branch is three places for one
// of them to be left behind.
float sound2(float r, float p){ return adia > 0.5 ? gam * p / max(r, smallr) : cs2; }

// Total energy from the primitives, and back again. b2 is |B|^2 at the cell, so the
// magnetic energy is in the conserved variable exactly as the reference keeps it --
// which is why the cooling pass has to subtract it again before it can talk about a
// temperature. Meaningless when adia = 0, where the slot carries nothing and every
// update leaves it alone.
float etot(float r, float u, float v, float p, float b2){
  return p / max(gam - 1.0, 1e-6) + 0.5 * r * (u * u + v * v) + 0.5 * b2;
}
float pfromE(float E, float r, float u, float v, float b2){
  return max((gam - 1.0) * (E - 0.5 * r * (u * u + v * v) - 0.5 * b2), smallr * 1e-4);
}

`;


  // ---- the flux, four ways -------------------------------------------------
  //
  // Shared by both godunov programs, because there is no reason for the PLM and PPM
  // paths to carry their own copy of a Riemann solver -- they did, identically, and
  // that is two places for a fix to be applied to one of.
  const RIEMANN = `
uniform float solver;

// find_speed_fast, isothermal: d2 is the half sum of the total Alfven speed squared
// and cs^2, and the normal field enters a second time because it is the component
// that can drive the discriminant to zero and make the fast and slow speeds meet.
float cfast(float r, float A2, float bt, float p){
  float ri = 1.0 / r;
  float c2 = sound2(r, p);
  float d2 = 0.5 * ((A2 + bt * bt) * ri + c2);
  return sqrt(d2 + sqrt(max(d2 * d2 - c2 * A2 * ri, 0.0)));
}

// The physical flux of one state, in the rotated frame, and the only place in the
// four where it is written down. Nothing magnetic is optional here and none of it
// depends on how many waves the solver believes in: the normal momentum carries the
// total pressure less A*A and the transverse momentum carries -A*Bt, which is the
// tension in the field lines threading the face. bn and ps come back empty because
// that pair is exact and is filled in once, above the branch.
// E arrives as an argument rather than being rebuilt from pt, and that is not
// fastidiousness: outside the fan E follows from the primitives, but in the star
// regions it does not. Miyoshi & Kusano's starred energies are their own expressions,
// carrying the work done by the rotational discontinuities, and recovering a gas
// pressure from pstar to feed back through etot() would quietly replace them with
// something that is not a solution of anything. Isothermal mode passes zero and the
// update ignores the slot.
Flux physFlux(float r, float u, float v, float y, float bt, float A, float pt, float E){
  float fd = r * u;
  return Flux(vec4(fd, fd * u + pt - A * A, fd * v - A * bt, fd * y),
              0.0, bt * u - A * v, 0.0,
              (E + pt) * u - A * (u * A + v * bt));
}

// Four fluxes behind one signature, chosen by a uniform so that the wave structure
// is a knob and not a build, and so that two of them can be compared with nothing
// else in the box moving. Each is named for what it gives up:
//
//   0  LLF / Rusanov. One wave, at the fastest speed anywhere in the problem, and
//      the same one for every variable: the mean of the two physical fluxes plus a
//      diffusion. It cannot tell that the dye is advected and the field is not.
//   1  HLL. Two waves, so the fan is bounded and a shock is captured, and one state
//      between them that serves rho, the transverse momentum, the transverse field
//      and the dye alike. With no contact in it, the transverse velocity and the
//      passive scalar are averaged across the interface rather than upwinded -- and a
//      shear layer is a discontinuity in exactly those and nothing else, so this
//      solver dissipates the one structure the page is about, at a rate the grid
//      picks rather than the physics.
//   2  HLLC. Three waves: the contact is back, so the transverse pair and the dye
//      come from whichever side of it the interface fell on. What is gone is the two
//      Alfven waves and the rotation of (u_t, B_t) they carry, so a rotational
//      discontinuity arrives as a step and is smeared at the rate the fast waves set.
//   3  HLLD. Five waves, hlld_mhd_fluxes reduced to isothermal, the default and what
//      the rest of the page is written around.
//
// They share more than the signature. The GLM pair is exact and is computed once,
// above the branch; the fast speed, the physical flux and the conserved vector are
// one function each; and HLLC here is HLLD with the two Alfven branches struck out
// and nothing else changed, which is what makes the degenerate limits exact -- at
// A = 0 the Alfven waves sit on the contact and those two return the same numbers to
// the last bit, and with the field gone entirely all four return what the live site's
// own hll and hllc return for the same states.
Flux riemann(Prim L, Prim R, float ch){
  float rl = max(L.h.x, smallr), rr = max(R.h.x, smallr);
  float ul = L.h.y, ur = R.h.y;
  float vl = L.h.z, vr = R.h.z;
  float bl = L.f.y, br = R.f.y;                 // transverse field
  float A  = 0.5 * (L.f.x + R.f.x);             // normal field: continuity enforced
  float A2 = A * A;
  // The gas pressure. Under a gamma law it is a reconstructed primitive and arrives in
  // the state; under the isothermal law it is re-derived here from the floored density
  // instead, and deliberately not read from the slot. Isothermal pressure is not an
  // independent variable -- it is cs^2 times the density at the same point -- and
  // letting it be limited, predicted and slab-averaged on its own would let the two
  // drift apart by exactly the difference between their source terms. Deriving it here
  // also means neither trace nor either of the two reconstructions had to change, and
  // that this path is the arithmetic it always was, to the bit.
  float pl  = adia > 0.5 ? L.f.w : rl * cs2;
  float pr_ = adia > 0.5 ? R.f.w : rr * cs2;
  float ptl = pl  + 0.5 * (A2 + bl * bl);       // total pressure, gas + magnetic
  float ptr = pr_ + 0.5 * (A2 + br * br);
  // and the conserved energies of the two physical states, which the fan's edges and
  // both diffusive solvers need
  float EL = etot(rl, ul, vl, pl,  A2 + bl * bl);
  float ER = etot(rr, ur, vr, pr_, A2 + br * br);

  // GLM, first and once. (Bn, psi) is a linear 2x2 system with eigenvalues +-ch and
  // no coupling to anything else, so its Riemann problem has an exact solution -- an
  // upwind average, which is what these two lines are -- and because it is exact
  // there is nothing here for a choice of flux to improve on or to spoil: all four
  // get the same pair. This is the term that makes flux%Bx nonzero, where a
  // constrained-transport scheme sets it to zero and lets the EMF own that component.
  Flux fx;
  fx.bn = 0.5 * (L.f.z + R.f.z) - 0.5 * ch * (R.f.x - L.f.x);
  fx.ps = ch * (0.5 * ch * (L.f.x + R.f.x) - 0.5 * (R.f.z - L.f.z));

  // the fast magnetosonic speed on each side. All four need it: it is the edge of the
  // fan, and for the first two it is the only speed in the problem.
  float cm = max(cfast(rl, A2, bl, pl), cfast(rr, A2, br, pr_));

  // LLF. One speed for the whole fan, so the dissipation is isotropic in state space
  // and every variable is damped at the rate the fastest wave sets, whether or not it
  // has anything to do with that wave.
  if (solver < 0.5) {
    float sm = max(abs(ul), abs(ur)) + cm;
    Flux FL = physFlux(rl, ul, vl, L.h.w, bl, A, ptl, EL);
    Flux FR = physFlux(rr, ur, vr, R.h.w, br, A, ptr, ER);
    vec4 UL = toCons(vec4(rl, ul, vl, L.h.w)), UR = toCons(vec4(rr, ur, vr, R.h.w));
    // the transverse field is a conserved variable in its own right, so the jump the
    // diffusion acts on is just br - bl. So is the total energy.
    fx.h  = 0.5 * (FL.h  + FR.h)  - 0.5 * sm * (UR - UL);
    fx.bt = 0.5 * (FL.bt + FR.bt) - 0.5 * sm * (br - bl);
    fx.en = 0.5 * (FL.en + FR.en) - 0.5 * sm * (ER - EL);
    return fx;
  }

  // Davis estimates, and deliberately not clamped through zero: HLL below wants only
  // that the fan contain the interface, and the two solvers after it need the real
  // signs to know which side of the contact it fell on.
  float SL = min(ul, ur) - cm;
  float SR = max(ul, ur) + cm;

  // HLL. The clamp is what makes one expression cover the supersonic cases too --
  // where the whole fan is on one side, the formula returns that side's physical flux.
  if (solver < 1.5) {
    float sl = min(SL, 0.0), sr = max(SR, 0.0);
    Flux FL = physFlux(rl, ul, vl, L.h.w, bl, A, ptl, EL);
    Flux FR = physFlux(rr, ur, vr, R.h.w, br, A, ptr, ER);
    vec4 UL = toCons(vec4(rl, ul, vl, L.h.w)), UR = toCons(vec4(rr, ur, vr, R.h.w));
    float iw = 1.0 / (sr - sl);
    fx.h  = (sr * FL.h  - sl * FR.h  + sr * sl * (UR - UL)) * iw;
    fx.bt = (sr * FL.bt - sl * FR.bt + sr * sl * (br - bl)) * iw;
    fx.en = (sr * FL.en - sl * FR.en + sr * sl * (ER - EL)) * iw;
    return fx;
  }

  // The contact, which the remaining two share. ustar and pstar follow from mass and
  // normal-momentum conservation across the whole fan, so they know nothing about the
  // Alfven waves and need no energy equation: the same two numbers serve three waves
  // and five.
  float rcl = rl * (ul - SL), rcr = rr * (SR - ur);
  float inv = 1.0 / (rcr + rcl);
  float ustar = (rcr * ur + rcl * ul + (ptl - ptr)) * inv;
  float pstar = (rcr * ptl + rcl * ptr + rcl * rcr * (ul - ur)) * inv;

  // left star region. The guard is the reference's: where estar vanishes the Alfven
  // wave has met the contact and the transverse state is simply continuous.
  float eps = max(1e-4 * A2, 1e-9);
  float rsl = max(rl * (SL - ul) / min(SL - ustar, -1e-8), smallr);
  float esl = rl * (SL - ul) * (SL - ustar) - A2;
  float ell = rl * (SL - ul) * (SL - ul)    - A2;
  float vsl = vl, bsl = bl;
  if (abs(esl) >= eps) {
    float ie = 1.0 / esl;
    vsl = vl - A * bl * (ustar - ul) * ie;
    bsl = bl * ell * ie;
  }

  // right star region
  float rsr = max(rr * (SR - ur) / max(SR - ustar, 1e-8), smallr);
  float esr = rr * (SR - ur) * (SR - ustar) - A2;
  float err = rr * (SR - ur) * (SR - ur)    - A2;
  float vsr = vr, bsr = br;
  if (abs(esr) >= eps) {
    float ie = 1.0 / esr;
    vsr = vr - A * br * (ustar - ur) * ie;
    bsr = br * err * ie;
  }

  // The starred energies, which exist only when there is an energy equation to put
  // them in. Miyoshi & Kusano: across the fast wave the energy jump is fixed by the
  // work done on the interface, and the term A*(v.B - v*.B*) is the Poynting flux
  // carried by the field threading it -- the piece that has no hydrodynamic analogue
  // and the reason this cannot be recovered from pstar alone. The divides are the same
  // guarded denominators the star densities above already use, so a degenerate fan
  // cannot produce a different sign here than it did there.
  float vbl = ul * A + vl * bl, vbr = ur * A + vr * br;
  float vbsl = ustar * A + vsl * bsl, vbsr = ustar * A + vsr * bsr;
  float Esl = ((SL - ul) * EL - ptl * ul + pstar * ustar + A * (vbl - vbsl))
            / min(SL - ustar, -1e-8);
  float Esr = ((SR - ur) * ER - ptr * ur + pstar * ustar + A * (vbr - vbsr))
            / max(SR - ustar, 1e-8);

  // and the one thing the two disagree about: what lies between the fast wave and the
  // contact on the side the interface is on.
  float ro, uo, vo, bo, pto, Eo;
  if (solver < 2.5) {
    // HLLC: outside the fan on either side, or the star state of whichever side of the
    // contact the interface fell on, taken whole. Upwinding the transverse pair rather
    // than averaging it is what the third wave buys; taking it whole is what is still
    // missing, since (u_t, B_t) does turn between the fast wave and the contact and the
    // waves that turn it are not in this fan. It is at least the pair the fast wave
    // leaves behind and not the far-field one, these being the same star states HLLD
    // builds.
    if (SL > 0.0)         { ro = rl;  uo = ul;    vo = vl;  bo = bl;  pto = ptl;   Eo = EL; }
    else if (ustar > 0.0) { ro = rsl; uo = ustar; vo = vsl; bo = bsl; pto = pstar; Eo = Esl; }
    else if (SR > 0.0)    { ro = rsr; uo = ustar; vo = vsr; bo = bsr; pto = pstar; Eo = Esr; }
    else                  { ro = rr;  uo = ur;    vo = vr;  bo = br;  pto = ptr;   Eo = ER; }
  } else {
    // HLLD: the two Alfven waves, and between them one rotated transverse state. That
    // is the pair of branches HLLC is missing and the reason a shear layer in B_t
    // survives here. |A|/sqrt(rho*) is the Alfven speed in the star region, so as the
    // normal field goes away the two waves close on the contact and these six branches
    // become the four above, term for term.
    float sgnm = A >= 0.0 ? 1.0 : -1.0;
    float sql = sqrt(rsl), SAL = ustar - abs(A) / sql;
    float sqr = sqrt(rsr), SAR = ustar + abs(A) / sqr;
    float den = 1.0 / (sql + sqr);
    float vss = (sql * vsl + sqr * vsr + sgnm * (bsr - bsl)) * den;
    float bss = (sql * bsr + sqr * bsl + sgnm * sql * sqr * (vsr - vsl)) * den;
    // The rotational discontinuities do work too, and it is the same Poynting term
    // once more -- now the jump in v.B across the Alfven wave, weighted by the star
    // density's root, with the sign of the normal field deciding which way it goes.
    float vbss = ustar * A + vss * bss;
    float Essl = Esl - sgnm * sql * (vbsl - vbss);
    float Essr = Esr + sgnm * sqr * (vbsr - vbss);
    if (SL > 0.0)         { ro = rl;  uo = ul;    vo = vl;  bo = bl;  pto = ptl;   Eo = EL; }
    else if (SAL > 0.0)   { ro = rsl; uo = ustar; vo = vsl; bo = bsl; pto = pstar; Eo = Esl; }
    else if (ustar > 0.0) { ro = rsl; uo = ustar; vo = vss; bo = bss; pto = pstar; Eo = Essl; }
    else if (SAR > 0.0)   { ro = rsr; uo = ustar; vo = vss; bo = bss; pto = pstar; Eo = Essr; }
    else if (SR > 0.0)    { ro = rsr; uo = ustar; vo = vsr; bo = bsr; pto = pstar; Eo = Esr; }
    else                  { ro = rr;  uo = ur;    vo = vr;  bo = br;  pto = ptr;   Eo = ER; }
  }

  // sampled at x/t = 0, and the flux is the physical one evaluated there. scalar_flux:
  // the dye rides on the mass flux and is upwinded on its sign, which in the star
  // region is the side of the contact the interface lies on.
  Flux f = physFlux(ro, uo, vo, ro * uo >= 0.0 ? L.h.w : R.h.w, bo, A, pto, Eo);
  fx.h = f.h; fx.bt = f.bt; fx.en = f.en;
  return fx;
}
`;

  // ---- Courant condition ---------------------------------------------------
  //
  // cmpdt, with the fast magnetosonic speed in place of cs. Two numbers come off
  // this reduction, not one:
  //
  //   x  ctot = sum over dim of (|u_dim| + cf_dim), which sets dt
  //   y  ch   = max over dim of (|u_dim| + cf_dim), the largest signal speed in
  //             the cell, which is what the GLM cleaning wave should travel at
  //
  // ch has to be the real maximum rather than ctot: it is a wave speed the
  // timestep must respect, and using the sum would either waste half the step or
  // -- worse, if the step were fixed regardless -- put the psi waves outside the
  // stability limit while every physical wave stayed inside it.
  const F_CMAX_STATE = HEAD + EOS + `
uniform sampler2D uSrc, uSrcB;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 8;
  vec2 m = vec2(0.0);
  for (int j = 0; j < 8; j++) {
    for (int i = 0; i < 8; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      vec4 U = texelFetch(uSrc, p, 0);
      vec4 B = texelFetch(uSrcB, p, 0);
      float r  = max(U.x, smallr);
      float ri = 1.0 / r;
      vec2  u  = U.yz * ri;
      float b2 = dot(B.xy, B.xy);
      // The real sound speed, which under a gamma law means recovering the pressure
      // from the energy channel before anything else happens. This reduction is the
      // Courant condition: an under-estimated wave speed is not a slightly wrong
      // number, it is a step outside the stability limit. And the under-estimate would
      // be large -- the warm phase sits at cs = 3.2 code units against the isothermal
      // 1, so ctot roughly doubles the moment the gas goes two-phase.
      float pg = adia > 0.5 ? pfromE(B.w, r, u.x, u.y, b2) : r * cs2;
      float c2 = sound2(r, pg);
      // find_speed_fast: d2 = (b^2/rho + c^2)/2, cf^2 = d2 + sqrt(d2^2 - c^2 bn^2/rho)
      float d2 = 0.5 * (b2 * ri + c2);
      float cx = sqrt(d2 + sqrt(max(d2 * d2 - c2 * B.x * B.x * ri, 0.0)));
      float cy = sqrt(d2 + sqrt(max(d2 * d2 - c2 * B.y * B.y * ri, 0.0)));
      float sx = abs(u.x) + cx, sy = abs(u.y) + cy;
      m = max(m, vec2(sx + sy, max(sx, sy)));
    }
  }
  outColor = vec4(m, 0.0, 1.0);
}`;

  const F_CMAX_DOWN = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 8;
  vec2 m = vec2(0.0);
  for (int j = 0; j < 8; j++) {
    for (int i = 0; i < 8; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      m = max(m, texelFetch(uSrc, p, 0).xy);
    }
  }
  outColor = vec4(m, 0.0, 1.0);
}`;

  // dt = min(dt_default, dt_adaptive), as on the live site: a fixed step chosen
  // inside the stability limit, with the Courant condition as a floor that engages
  // only when the flow gets faster than the state it was chosen for. The default is
  // lower here than in the hydrodynamic build because ctot is about twice as large
  // at the same Mach number -- at beta = 1 the fast speed is sqrt(cs^2 + vA^2) =
  // 1.73 across the field and vA = 1.41 along it, against cs = 1.
  const DTDX = `
uniform sampler2D uCmax;
uniform float cfl, smallc, dtdxMax;
vec2 waves(){ return max(texelFetch(uCmax, ivec2(0, 0), 0).xy, vec2(1e-20)); }
float stepDtDx(){ return min(dtdxMax, min(cfl / smallc, cfl / waves().x)); }
float chSpeed(){ return waves().y; }
`;

  // ---- the solver: one unsplit pass, two targets ---------------------------
  const F_GODUNOV = HEAD2 + EOS + RIEMANN + DTDX + `
uniform sampler2D uU, uB;
uniform ivec2 size;
uniform float dyeDiss, dxCell, slopeType, psiDamp, powell, fric;

Prim prim(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  vec4 U = texelFetch(uU, p, 0);
  vec4 B = texelFetch(uB, p, 0);
  float r = max(U.x, smallr), ri = 1.0 / r;
  // The field target's fourth channel is the conserved total energy under a gamma law,
  // and the primitive slot wants a pressure, so it is recovered here once per cell
  // rather than at every use. Isothermal mode leaves the channel alone and riemann()
  // derives the pressure from the density itself.
  float pg = adia > 0.5 ? pfromE(B.w, r, U.y * ri, U.z * ri, dot(B.xy, B.xy)) : r * cs2;
  return Prim(vec4(r, U.y * ri, U.z * ri, U.w * ri), vec4(B.xyz, pg));
}

// uslope, slope_type = 2: MonCen. dlft/drgt are the one-sided differences scaled
// by slope_type, dcen the central one; the limiter takes the smaller one-sided
// slope and zero where they disagree in sign.
vec4 uslope(vec4 qm, vec4 q0, vec4 qp){
  vec4 dlft = slopeType * (q0 - qm);
  vec4 drgt = slopeType * (qp - q0);
  vec4 dcen = 0.5 * (dlft + drgt) / max(slopeType, 1.0);
  vec4 dlim = min(abs(dlft), abs(drgt)) * step(vec4(0.0), dlft * drgt);
  return sign(dcen) * min(dlim, abs(dcen));
}

// The half-step predictor, trace2d's job, for isothermal GLM-MHD. Derived rather
// than transcribed: the reference's MHD trace predicts a staggered field through
// EMFs, which a cell-centred scheme has no use for. In primitive form,
//
//   D rho/Dt = -rho div u
//   D u  /Dt = (-P_x - By Jz) / rho               Jz = By_x - Bx_y
//   D v  /Dt = (-P_y + Bx Jz) / rho               so the force is (J x B)/rho
//   D Y  /Dt = 0
//   D P  /Dt = -gamma P div u
//   dBx/dt   = u_y By + u By_y - v_y Bx - v Bx_y - psi_x
//   dBy/dt   = -u_x By - u By_x + v_x Bx + v Bx_x - psi_y
//   dpsi/dt  = -ch^2 div B
//
// P_x is cs^2 rho_x under the isothermal law and the pressure's own slope under a
// gamma law, and the pressure line is the one the isothermal page
// did not have. Its first two terms are the advection every primitive gets and the
// third is the work of compression, which is the whole reason a gamma-law gas can go
// two-phase: squeeze it and it heats, and then the cooling has something to argue
// with. Isothermal mode carries a pressure slot as well, but riemann() re-derives it
// from the density there and never reads what is written here.
//
// The magnetic pressure gradient and the field-line tension have been combined
// into J x B, which is exact and cheaper: the two cancel along B. hx/hy and gx/gy
// are half slopes, as drx/dry are in the reference, so the whole source is
// multiplied by dt/dx once.
void trace(inout Prim q, vec4 hx, vec4 hy, vec4 gx, vec4 gy, float dtdx, float ch2){
  float r = max(q.h.x, smallr), u = q.h.y, v = q.h.z, ri = 1.0 / r;
  float bx = q.f.x, by = q.f.y;
  float jz = gx.y - gy.x;
  vec4 sh, sf;
  sh.x = -u * hx.x - v * hy.x - (hx.y + hy.z) * r;
#ifdef ISOTHERMAL_ONLY
  sh.y = -u * hx.y - v * hy.y + (-cs2 * hx.x - by * jz) * ri;
  sh.z = -u * hx.z - v * hy.z + (-cs2 * hy.x + bx * jz) * ri;
#else
  sh.y = -u * hx.y - v * hy.y + (-gx.w - by * jz) * ri;
  sh.z = -u * hx.z - v * hy.z + (-gy.w + bx * jz) * ri;
#endif
  sh.w = -u * hx.w - v * hy.w;
  sf.x =  hy.y * by + u * gy.y - hy.z * bx - v * gy.x - gx.z;
  sf.y = -hx.y * by - u * gx.y + hx.z * bx + v * gx.x - gy.z;
  sf.z = -ch2 * (gx.x + gy.y);
#ifdef ISOTHERMAL_ONLY
  sf.w = 0.0;
#else
  sf.w = -u * gx.w - v * gy.w - gam * q.f.w * (hx.y + hy.z);
#endif
  vec4 oh = q.h + sh * dtdx;
  // the reference falls back to the unpredicted value if the predictor undershoots
  // the density floor, rather than clamping to it
  oh.x = oh.x < smallr ? q.h.x : oh.x;
  q.h = oh;
  q.f = q.f + sf * dtdx;
}



vec4 swapv(vec4 a){ return vec4(a.x, a.z, a.y, a.w); }
// the y sweep, in the same frame the x sweep uses: normal velocity in .y, normal
// field in .x, exactly as the reference always puts the normal in velocity_x
Prim rot(Prim q){ return Prim(swapv(q.h), vec4(q.f.y, q.f.x, q.f.z, q.f.w)); }

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  const ivec2 ex = ivec2(1, 0), ey = ivec2(0, 1);

  Prim q0 = prim(c);
  Prim qE = prim(c + ex),      qW = prim(c - ex);
  Prim qN = prim(c + ey),      qS = prim(c - ey);
  Prim qEE = prim(c + 2 * ex), qWW = prim(c - 2 * ex);
  Prim qNN = prim(c + 2 * ey), qSS = prim(c - 2 * ey);
  Prim qNE = prim(c + ex + ey), qSE = prim(c + ex - ey);
  Prim qNW = prim(c - ex + ey), qSW = prim(c - ex - ey);

  // half slopes on the five-cell cross, both directions and both halves of the
  // state: the predictor couples all of it
  vec4 h0x = 0.5 * uslope(qW.h, q0.h, qE.h),   h0y = 0.5 * uslope(qS.h, q0.h, qN.h);
  vec4 g0x = 0.5 * uslope(qW.f, q0.f, qE.f),   g0y = 0.5 * uslope(qS.f, q0.f, qN.f);
  vec4 hEx = 0.5 * uslope(q0.h, qE.h, qEE.h),  hEy = 0.5 * uslope(qSE.h, qE.h, qNE.h);
  vec4 gEx = 0.5 * uslope(q0.f, qE.f, qEE.f),  gEy = 0.5 * uslope(qSE.f, qE.f, qNE.f);
  vec4 hWx = 0.5 * uslope(qWW.h, qW.h, q0.h),  hWy = 0.5 * uslope(qSW.h, qW.h, qNW.h);
  vec4 gWx = 0.5 * uslope(qWW.f, qW.f, q0.f),  gWy = 0.5 * uslope(qSW.f, qW.f, qNW.f);
  vec4 hNx = 0.5 * uslope(qNW.h, qN.h, qNE.h), hNy = 0.5 * uslope(q0.h, qN.h, qNN.h);
  vec4 gNx = 0.5 * uslope(qNW.f, qN.f, qNE.f), gNy = 0.5 * uslope(q0.f, qN.f, qNN.f);
  vec4 hSx = 0.5 * uslope(qSW.h, qS.h, qSE.h), hSy = 0.5 * uslope(qSS.h, qS.h, q0.h);
  vec4 gSx = 0.5 * uslope(qSW.f, qS.f, qSE.f), gSy = 0.5 * uslope(qSS.f, qS.f, q0.f);

  float dtdx = stepDtDx();
  float ch = chSpeed(), ch2 = ch * ch;

  trace(q0, h0x, h0y, g0x, g0y, dtdx, ch2);
  trace(qE, hEx, hEy, gEx, gEy, dtdx, ch2);
  trace(qW, hWx, hWy, gWx, gWy, dtdx, ch2);
  trace(qN, hNx, hNy, gNx, gNy, dtdx, ch2);
  trace(qS, hSx, hSy, gSx, gSy, dtdx, ch2);

  Prim e0 = Prim(q0.h + h0x, q0.f + g0x), e1 = Prim(qE.h - hEx, qE.f - gEx);
  Prim w0 = Prim(qW.h + hWx, qW.f + gWx), w1 = Prim(q0.h - h0x, q0.f - g0x);
  Prim n0 = Prim(q0.h + h0y, q0.f + g0y), n1 = Prim(qN.h - hNy, qN.f - gNy);
  Prim s0 = Prim(qS.h + hSy, qS.f + gSy), s1 = Prim(q0.h - h0y, q0.f - g0y);

  Flux Fe = riemann(e0, e1, ch);
  Flux Fw = riemann(w0, w1, ch);
  Flux Fn = riemann(rot(n0), rot(n1), ch);
  Flux Fs = riemann(rot(s0), rot(s1), ch);

  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);

  vec4 dU = (Fw.h - Fe.h) + swapv(Fs.h - Fn.h);
  // the transverse flux of one component is the normal (psi) flux of the other
  vec2 dB = (vec2(Fw.bn, Fw.bt) - vec2(Fe.bn, Fe.bt))
          + (vec2(Fs.bt, Fs.bn) - vec2(Fn.bt, Fn.bn));
  float dpsi = (Fw.ps - Fe.ps) + (Fs.ps - Fn.ps);
  vec4 Un = U + dU * dtdx;
#ifdef ISOTHERMAL_ONLY
  float divb = (0.5 * (e0.f.x + e1.f.x) - 0.5 * (w0.f.x + w1.f.x))
             + (0.5 * (n0.f.y + n1.f.y) - 0.5 * (s0.f.y + s1.f.y));
  Un.yz -= powell * divb * B.xy * dtdx;
  Un.yz = Un.yz / (1.0 + fric * dtdx * dxCell);
  Un.x = max(Un.x, smallr);
  Un.w = Un.w / (1.0 + dyeDiss * dtdx * dxCell);
  float psi = (B.z + dpsi * dtdx) * exp(-psiDamp * ch * dtdx);
  outColor  = Un;
  outColor1 = vec4(B.xy + dB * dtdx, psi, 0.0);
#else
  // the energy flux is a scalar, so the rotated sweep's contributes as it stands
  float dEn = (Fw.en - Fe.en) + (Fs.en - Fn.en);
  vec2 Bn = B.xy + dB * dtdx;
  // ------------------------------------------------------------------ the energy
  //
  // The fluxes give a total energy at the new time, and then three momentum sinks act
  // below: the Powell term, the large-scale drag, the density floor. How the energy
  // answers them is a physical choice and not bookkeeping, so it is made explicitly.
  // The thermal part is taken out here, before any of them, and put back after -- so
  // the kinetic energy they remove leaves the box rather than turning into heat.
  //
  // For the drag that is the point of it. FRIC stands in for the three-dimensional
  // cascade this plane does not have, and a cascade carries energy away to scales that
  // are not in the box: it does not warm the gas it took the energy from. Letting the
  // total energy ride through the division would have made it a viscosity instead, and
  // a viscosity at that strength -- a box-scale timescale -- would heat the gas faster
  // than the cooling could take it away, which is the two-phase structure gone.
  //
  // For Powell it is not an approximation either, and pleasingly so. Holding the
  // thermal energy fixed across a momentum kick changes the total by v.dm, which for
  // dm = -(div B) B dt is -(div B)(v.B) dt -- exactly the energy source term the
  // reference's Powell scheme carries alongside the momentum one. It comes out of the
  // bookkeeping rather than having to be added to it.
  float rf = max(Un.x, smallr);
  float rfi = 1.0 / rf;
  float pth = pfromE(B.w + dEn * dtdx, rf, Un.y * rfi, Un.z * rfi, dot(Bn, Bn));
  // Powell. The divergence here is the one the psi flux sees -- the same
  // face-averaged normal fields -- so the force the momentum equation loses is
  // exactly the one the cleaning is working on.
  float divb = (0.5 * (e0.f.x + e1.f.x) - 0.5 * (w0.f.x + w1.f.x))
             + (0.5 * (n0.f.y + n1.f.y) - 0.5 * (s0.f.y + s1.f.y));
  Un.yz -= powell * divb * B.xy * dtdx;
  // A large-scale drag, backward Euler so it is stable at any step. Two dimensions
  // do not have three's energy cascade: enstrophy goes to the grid and energy goes
  // the other way, into box-scale coherent structures that nothing dissipates.
  // Measured, that is not a small effect here -- with the drive cut by a factor of a
  // hundred the rms Mach held at 0.85 for half a minute while the peaks decayed, so
  // the box was storing everything it had ever been given at the largest scale and
  // the servo had no plant to act on. This is the sink that gives it one, and it is
  // what a two-dimensional driven-turbulence run does for the same reason. The
  // timescale is longer than a box-scale turnover, so it sets the equilibrium
  // without touching the eddies.
  Un.yz = Un.yz / (1.0 + fric * dtdx * dxCell);
  Un.x = max(Un.x, smallr);
  Un.w = Un.w / (1.0 + dyeDiss * dtdx * dxCell);

  // and the parabolic half of Dedner's operator: psi is damped as it travels, so
  // the divergence error is not merely moved somewhere else in the box
  float psi = (B.z + dpsi * dtdx) * exp(-psiDamp * ch * dtdx);

  // and the energy reassembled: the thermal part the sinks were not allowed to touch,
  // the kinetic part they did, the magnetic part the induction equation just wrote.
  float En = etot(Un.x, Un.y / Un.x, Un.z / Un.x, pth, dot(Bn, Bn));

  outColor  = Un;
  outColor1 = vec4(Bn, psi, En);
#endif
}`;

  // ---- third order in space: PPM for MHD ----------------------------------
  //
  // Colella & Woodward 1984 on all seven primitives, which is not in the reference
  // -- its slope_type 0 through 6 are every one of them piecewise linear -- so the
  // PLM program above, which *is* the reference's configuration, stays as the
  // alternative and as what runs while the box is being driven.
  //
  // Steps 1 to 3 are the hydrodynamic page's, componentwise, and the field rides
  // through them unchanged: limited MonCen slopes, the fourth-order interface value
  // a(i+1/2) = (a_i + a_i+1)/2 + (dq_i - dq_i+1)/6 built on them, and CW84 (1.10)
  // monotonisation -- flatten at an extremum, pull an overshooting edge back until
  // the parabola is monotone across the cell.
  //
  // Step 4 is where MHD differs, and getting it wrong would quietly undo the
  // divergence cleaning. The face state in PPM is the parabola averaged over the
  // slab the flow crosses in half a step, which is an advective time centring: it
  // carries the -u dq/dx term for you, at one speed for every variable. That is an
  // approximation for the hydrodynamic variables, all of which *are* advected along
  // the normal. It is simply wrong for two of the three field variables, because
  //
  //   dBn/dt  has no  -u dBn/dx  term at all -- the normal component is advected
  //           only transversely, which is the same fact constrained transport is
  //           built around
  //   dpsi/dt = -ch^2 div B  is not advected by anything
  //   dBt/dt  does have -u dBt/dx, so the transverse component is the one that
  //           belongs in the slab average
  //
  // Shifting Bn and psi by u dt/2 anyway would inject a spurious -u dBn/dx of
  // exactly the order of the terms being added, straight into the divergence the
  // cleaning is trying to damp. So the field takes the slab average in its
  // transverse component and the plain parabola edge in the other two -- one mix
  // against a constant mask, no extra cost -- and its source term keeps every term
  // the predictor has except the one the average now carries.
  const F_GODUNOV3 = HEAD2 + EOS + RIEMANN + DTDX + `
uniform sampler2D uU, uB;
uniform ivec2 size;
uniform float dyeDiss, dxCell, slopeType, psiDamp, powell, fric;

Prim prim(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  vec4 U = texelFetch(uU, p, 0);
  vec4 B = texelFetch(uB, p, 0);
  float r = max(U.x, smallr), ri = 1.0 / r;
  // The field target's fourth channel is the conserved total energy under a gamma law,
  // and the primitive slot wants a pressure, so it is recovered here once per cell
  // rather than at every use. Isothermal mode leaves the channel alone and riemann()
  // derives the pressure from the density itself.
  float pg = adia > 0.5 ? pfromE(B.w, r, U.y * ri, U.z * ri, dot(B.xy, B.xy)) : r * cs2;
  return Prim(vec4(r, U.y * ri, U.z * ri, U.w * ri), vec4(B.xyz, pg));
}

vec4 uslope(vec4 qm, vec4 q0, vec4 qp){
  vec4 dlft = slopeType * (q0 - qm);
  vec4 drgt = slopeType * (qp - q0);
  vec4 dcen = 0.5 * (dlft + drgt) / max(slopeType, 1.0);
  vec4 dlim = min(abs(dlft), abs(drgt)) * step(vec4(0.0), dlft * drgt);
  return sign(dcen) * min(dlim, abs(dcen));
}

// The monotonised parabola for one cell, from slopes computed by the caller.
// Adjacent cells' windows overlap by two, so the five slopes are computed once per
// sweep and passed in: that saves eight of eighteen uslope evaluations per variable
// group, and with seven variables uslope is the most expensive thing in the pass.
void ppmEdges(vec4 qm, vec4 q0, vec4 qp, vec4 dm, vec4 d0, vec4 dp, out vec4 aL, out vec4 aR){
  aL = 0.5 * (qm + q0) + (dm - d0) * (1.0 / 6.0);
  aR = 0.5 * (q0 + qp) + (d0 - dp) * (1.0 / 6.0);
  vec4 ext = step((aR - q0) * (q0 - aL), vec4(0.0));
  aL = mix(aL, q0, ext);
  aR = mix(aR, q0, ext);
  vec4 d = aR - aL;
  vec4 m = q0 - 0.5 * (aL + aR);
  vec4 c1 = step(d * d * (1.0 / 6.0), d * m);
  vec4 c2 = step(d * m, -d * d * (1.0 / 6.0));
  aL = mix(aL, 3.0 * q0 - 2.0 * aR, c1);
  aR = mix(aR, 3.0 * q0 - 2.0 * aL, c2);
}

// The parabola's average over the slab of width |u| dt/2 adjoining a face.
vec4 ppmFace(vec4 q0, vec4 aL, vec4 aR, float sig, float right){
  vec4 d = aR - aL;
  vec4 a6 = 6.0 * (q0 - 0.5 * (aL + aR));
  float k = 1.0 - (2.0 / 3.0) * sig;
  return right > 0.5 ? aR - 0.5 * sig * (d - k * a6)
                     : aL + 0.5 * sig * (d + k * a6);
}

// Which field components are advected along the sweep normal, and so belong in the
// slab average: the transverse one, and only that, in the shipped isothermal path.
// Pressure joins it under the gamma law.
const vec4 ADVECTS = vec4(0.0, 1.0, 0.0, 0.0);
const vec4 ADVECTS_ADIA = vec4(0.0, 1.0, 0.0, 1.0);

vec4 floorRho(vec4 f, vec4 q){ return f.x < smallr ? vec4(q.x, f.yzw) : f; }

// The predictor, minus whatever the parabola average already carries. hx is the
// parabola's effective half slope along the normal, hy the transverse one.
vec4 srcH(vec4 h, vec4 f, vec4 hx, vec4 hy, vec4 gx, vec4 gy){
  float r = max(h.x, smallr), v = h.z, ri = 1.0 / r;
  float jz = gx.y - gy.x;
  vec4 s;
  s.x = -v * hy.x - (hx.y + hy.z) * r;
#ifdef ISOTHERMAL_ONLY
  s.y = -v * hy.y + (-cs2 * hx.x - f.y * jz) * ri;
  s.z = -v * hy.z + (-cs2 * hy.x + f.x * jz) * ri;
#else
  s.y = -v * hy.y + (-gx.w - f.y * jz) * ri;
  s.z = -v * hy.z + (-gy.w + f.x * jz) * ri;
#endif
  s.w = -v * hy.w;
  return s;
}
vec4 srcF(vec4 h, vec4 f, vec4 hx, vec4 hy, vec4 gx, vec4 gy, float ch2){
  float u = h.y, v = h.z;
  vec4 s;
  s.x =  hy.y * f.y + u * gy.y - hy.z * f.x - v * gy.x - gx.z;
  s.y = -hx.y * f.y             + hx.z * f.x + v * gx.x - gy.z;   // -u*gx.y dropped
  s.z = -ch2 * (gx.x + gy.y);
  // The pressure, on the same convention: -gamma P div u and the transverse advection
  // stay, and -u dP/dx is dropped because the slab average now carries it -- the
  // pressure having joined the advection mask. Isothermal mode never reads this.
#ifdef ISOTHERMAL_ONLY
  s.w = 0.0;
#else
  s.w = -v * gy.w - gam * f.w * (hx.y + hy.z);
#endif
  return s;
}



vec4 swapv(vec4 a){ return vec4(a.x, a.z, a.y, a.w); }
Prim rot(Prim q){ return Prim(swapv(q.h), vec4(q.f.y, q.f.x, q.f.z, q.f.w)); }

// One sweep, in the frame where the normal velocity is in .y and the normal field
// in .x. a3..b3 run along the normal; the six t* are the transverse neighbours of
// the three cells that own the two faces. An* come back so the caller can build the
// same face-averaged divergence the psi flux sees.
void sweep(Prim am3, Prim am2, Prim am1, Prim a0, Prim ap1, Prim ap2, Prim ap3,
           Prim tmm, Prim tpm, Prim tm0, Prim tp0, Prim tmp, Prim tpp,
           float dtdx, float ch, out Flux Flo, out Flux Fhi,
           out float Anlo, out float Anhi){
  vec4 h2m = uslope(am3.h, am2.h, am1.h), f2m = uslope(am3.f, am2.f, am1.f);
  vec4 h1m = uslope(am2.h, am1.h, a0.h),  f1m = uslope(am2.f, am1.f, a0.f);
  vec4 h00 = uslope(am1.h, a0.h,  ap1.h), f00 = uslope(am1.f, a0.f,  ap1.f);
  vec4 h1p = uslope(a0.h,  ap1.h, ap2.h), f1p = uslope(a0.f,  ap1.f, ap2.f);
  vec4 h2p = uslope(ap1.h, ap2.h, ap3.h), f2p = uslope(ap1.f, ap2.f, ap3.f);

  vec4 hlL, hlR, hcL, hcR, hrL, hrR, flL, flR, fcL, fcR, frL, frR;
  ppmEdges(am2.h, am1.h, a0.h,  h2m, h1m, h00, hlL, hlR);
  ppmEdges(am1.h, a0.h,  ap1.h, h1m, h00, h1p, hcL, hcR);
  ppmEdges(a0.h,  ap1.h, ap2.h, h00, h1p, h2p, hrL, hrR);
  ppmEdges(am2.f, am1.f, a0.f,  f2m, f1m, f00, flL, flR);
  ppmEdges(am1.f, a0.f,  ap1.f, f1m, f00, f1p, fcL, fcR);
  ppmEdges(a0.f,  ap1.f, ap2.f, f00, f1p, f2p, frL, frR);

  vec4 hyL = 0.5 * uslope(tmm.h, am1.h, tpm.h), gyL = 0.5 * uslope(tmm.f, am1.f, tpm.f);
  vec4 hy0 = 0.5 * uslope(tm0.h, a0.h,  tp0.h), gy0 = 0.5 * uslope(tm0.f, a0.f,  tp0.f);
  vec4 hyR = 0.5 * uslope(tmp.h, ap1.h, tpp.h), gyR = 0.5 * uslope(tmp.f, ap1.f, tpp.f);

  float sL = clamp(abs(am1.h.y) * dtdx, 0.0, 1.0);
  float s0 = clamp(abs(a0.h.y)  * dtdx, 0.0, 1.0);
  float sR = clamp(abs(ap1.h.y) * dtdx, 0.0, 1.0);

  float ch2 = ch * ch;
  vec4 hxL = 0.5 * (hlR - hlL), gxL = 0.5 * (flR - flL);
  vec4 hx0 = 0.5 * (hcR - hcL), gx0 = 0.5 * (fcR - fcL);
  vec4 hxR = 0.5 * (hrR - hrL), gxR = 0.5 * (frR - frL);

  vec4 ehL = srcH(am1.h, am1.f, hxL, hyL, gxL, gyL) * dtdx;
  vec4 ehC = srcH(a0.h,  a0.f,  hx0, hy0, gx0, gy0) * dtdx;
  vec4 ehR = srcH(ap1.h, ap1.f, hxR, hyR, gxR, gyR) * dtdx;
  vec4 efL = srcF(am1.h, am1.f, hxL, hyL, gxL, gyL, ch2) * dtdx;
  vec4 efC = srcF(a0.h,  a0.f,  hx0, hy0, gx0, gy0, ch2) * dtdx;
  vec4 efR = srcF(ap1.h, ap1.f, hxR, hyR, gxR, gyR, ch2) * dtdx;

  // The reference floors its reconstructed states, and a parabola through a fresh
  // discontinuity can undershoot: a face density at the floor beside a large
  // velocity makes rcl tiny, ustar enormous and the flux garbage.
#ifdef ISOTHERMAL_ONLY
  vec4 adv = ADVECTS;
#else
  vec4 adv = ADVECTS_ADIA;
#endif
  Prim lo0 = Prim(floorRho(ppmFace(am1.h, hlL, hlR, sL, 1.0) + ehL, am1.h),
                  mix(flR, ppmFace(am1.f, flL, flR, sL, 1.0), adv) + efL);
  Prim lo1 = Prim(floorRho(ppmFace(a0.h,  hcL, hcR, s0, 0.0) + ehC, a0.h),
                  mix(fcL, ppmFace(a0.f,  fcL, fcR, s0, 0.0), adv) + efC);
  Prim hi0 = Prim(floorRho(ppmFace(a0.h,  hcL, hcR, s0, 1.0) + ehC, a0.h),
                  mix(fcR, ppmFace(a0.f,  fcL, fcR, s0, 1.0), adv) + efC);
  Prim hi1 = Prim(floorRho(ppmFace(ap1.h, hrL, hrR, sR, 0.0) + ehR, ap1.h),
                  mix(frL, ppmFace(ap1.f, frL, frR, sR, 0.0), adv) + efR);

  Anlo = 0.5 * (lo0.f.x + lo1.f.x);
  Anhi = 0.5 * (hi0.f.x + hi1.f.x);
  Flo = riemann(lo0, lo1, ch);
  Fhi = riemann(hi0, hi1, ch);
}

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  const ivec2 ex = ivec2(1, 0), ey = ivec2(0, 1);

  Prim q0 = prim(c);
  Prim x1p = prim(c + ex),     x1m = prim(c - ex);
  Prim x2p = prim(c + 2 * ex), x2m = prim(c - 2 * ex);
  Prim x3p = prim(c + 3 * ex), x3m = prim(c - 3 * ex);
  Prim y1p = prim(c + ey),     y1m = prim(c - ey);
  Prim y2p = prim(c + 2 * ey), y2m = prim(c - 2 * ey);
  Prim y3p = prim(c + 3 * ey), y3m = prim(c - 3 * ey);
  Prim dNE = prim(c + ex + ey), dSE = prim(c + ex - ey);
  Prim dNW = prim(c - ex + ey), dSW = prim(c - ex - ey);

  float dtdx = stepDtDx();
  float ch = chSpeed();

  Flux Fw, Fe; float Aw, Ae;
  sweep(x3m, x2m, x1m, q0, x1p, x2p, x3p,
        dSW, dNW, y1m, y1p, dSE, dNE, dtdx, ch, Fw, Fe, Aw, Ae);

  Flux Fs, Fn; float As, An;
  sweep(rot(y3m), rot(y2m), rot(y1m), rot(q0), rot(y1p), rot(y2p), rot(y3p),
        rot(dSW), rot(dSE), rot(x1m), rot(x1p), rot(dNW), rot(dNE),
        dtdx, ch, Fs, Fn, As, An);

  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);

  vec4 dU = (Fw.h - Fe.h) + swapv(Fs.h - Fn.h);
  vec2 dB = (vec2(Fw.bn, Fw.bt) - vec2(Fe.bn, Fe.bt))
          + (vec2(Fs.bt, Fs.bn) - vec2(Fn.bt, Fn.bn));
  float dpsi = (Fw.ps - Fe.ps) + (Fs.ps - Fn.ps);
  vec4 Un = U + dU * dtdx;
#ifdef ISOTHERMAL_ONLY
  Un.yz -= powell * ((Ae - Aw) + (An - As)) * B.xy * dtdx;
  Un.yz = Un.yz / (1.0 + fric * dtdx * dxCell);
  Un.x = max(Un.x, smallr);
  Un.w = Un.w / (1.0 + dyeDiss * dtdx * dxCell);
  float psi = (B.z + dpsi * dtdx) * exp(-psiDamp * ch * dtdx);
  outColor  = Un;
  outColor1 = vec4(B.xy + dB * dtdx, psi, 0.0);
#else
  float dEn = (Fw.en - Fe.en) + (Fs.en - Fn.en);
  vec2 Bn = B.xy + dB * dtdx;
  // the thermal energy, held out of the way of the momentum sinks: the note in
  // F_GODUNOV's update says why, and it applies word for word here
  float rf = max(Un.x, smallr);
  float rfi = 1.0 / rf;
  float pth = pfromE(B.w + dEn * dtdx, rf, Un.y * rfi, Un.z * rfi, dot(Bn, Bn));
  Un.yz -= powell * ((Ae - Aw) + (An - As)) * B.xy * dtdx;
  Un.yz = Un.yz / (1.0 + fric * dtdx * dxCell);
  Un.x = max(Un.x, smallr);
  Un.w = Un.w / (1.0 + dyeDiss * dtdx * dxCell);
  float psi = (B.z + dpsi * dtdx) * exp(-psiDamp * ch * dtdx);
  float En = etot(Un.x, Un.y / Un.x, Un.z / Un.x, pth, dot(Bn, Bn));

  outColor  = Un;
  outColor1 = vec4(Bn, psi, En);
#endif
}`;

  // ---- forcing ------------------------------------------------------------
  //
  // Every pass that writes the state writes both targets, so the field rides
  // through untouched rather than going stale behind a ping-pong swap. It costs one
  // extra texel read and write on a grid of fifty thousand, which is nothing beside
  // the solver.
  // adia rather than the whole EOS block, because these four need one uniform out of
  // it and three of them declare smallr themselves -- and a duplicate declaration is
  // a link failure, which on this page is a black screen.
  const PASS = `
uniform sampler2D uB;
uniform float adia;

// Kinetic energy of a momentum density, on a density the caller has already floored.
float kinE(vec2 m, float r){ return 0.5 * dot(m, m) / max(r, 1e-30); }

// The field, and the energy channel corrected for the work this pass just did.
//
// Passing a total energy through a pass that changed the velocity is not leaving the
// gas alone: it holds the sum fixed while the kinetic part moves, so the difference
// comes silently out of the thermal part. Stirring the box would chill it, and stopping
// the box would set it on fire. What a stirring force does is work -- it adds kinetic
// energy at fixed temperature -- so the channel is moved by exactly the change in
// kinetic energy and the pressure is what survives the pass. That is also the right
// answer for a blast, which changes the density too: at fixed pressure adding mass
// cools the gas, and if the cooling is on that pushes the compressed region towards the
// cold branch, which is what a real compression does.
//
// With adia = 0 the slot carries nothing and this is the passthrough it always was.
vec4 passField(float dKin){
  vec4 F = texelFetch(uB, ivec2(gl_FragCoord.xy), 0);
  F.w += adia > 0.5 ? dKin : 0.0;
  return F;
}
`;

  const F_STIR = HEAD2 + PASS + DTDX + `
uniform sampler2D uU;
uniform vec2  c0, c1;
uniform vec3  dye;
uniform float amp, radius, dyeRadius, zeta, aspect, dxCell, smallr;

vec2 kick(vec2 p, vec2 c, float dir){
  vec2 d = p - c;
  d.x *= aspect;
  float r2 = dot(d, d);
  float w  = exp(-r2 / radius);
  float r  = sqrt(r2) + 1e-5;
  vec2 t   = vec2(-d.y, d.x) / r;          // tangential: divergence-free
  vec2 g   = d / r;                        // radial: curl-free
  return dir * w * mix(t, g, zeta);
}

void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 a = kick(vUv, c0, 1.0) + kick(vUv, c1, -1.0);
  float dt = stepDtDx() * dxCell;
  float r  = max(U.x, smallr);
  vec2 m = U.yz + r * a * amp * dt;
  float w = exp(-dot((vUv - c0) * vec2(aspect, 1.0), (vUv - c0) * vec2(aspect, 1.0)) / dyeRadius)
          + exp(-dot((vUv - c1) * vec2(aspect, 1.0), (vUv - c1) * vec2(aspect, 1.0)) / dyeRadius);
  outColor  = vec4(U.x, m, U.w + r * dye.x * w * dt);
  outColor1 = passField(kinE(m, r) - kinE(U.yz, r));
}`;

  const F_SPLAT = HEAD2 + PASS + `
uniform sampler2D uU;
uniform vec2  point, delta;
uniform float aspect, radius, ink, smallr;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float w = exp(-dot(p, p) / radius);
  float r = max(U.x, smallr);
  vec2 m = U.yz + delta * w * r;
  outColor  = vec4(U.x, m, U.w + r * ink * w);
  outColor1 = passField(kinE(m, r) - kinE(U.yz, r));
}`;

  // A tanh shear layer with a sinusoidal seed, as on the live site. Whether it
  // rolls up is now a race with the field as well as with the ambient strain --
  // see the note at the top of this file.
  const F_SHEAR = HEAD2 + PASS + `
uniform sampler2D uU;
uniform vec3  band;
uniform float amp, width, seed, phase, cycles, smallr;
void main(){
  vec4  U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  float r = max(U.x, smallr);
  float y = vUv.y - 0.5;
  float env = exp(-pow(y / (2.6 * width), 2.0));
  vec2  du = vec2(amp * tanh(y / width) * env,
                  seed * sin(6.2831853 * (vUv.x * cycles + phase)) * env);
  float bw = exp(-pow(y / (1.5 * width), 2.0));
  vec2 m = U.yz + r * du;
  outColor  = vec4(U.x, m, U.w + r * band.x * bw);
  outColor1 = passField(kinE(m, r) - kinE(U.yz, r));
}`;

  const F_BLAST = HEAD2 + PASS + `
uniform sampler2D uU;
uniform vec2  point;
uniform float aspect, radius, dens, kick, smallr;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  vec2 p = vUv - point;
  p.x *= aspect;
  float rr = length(p) + 1e-5;
  float w  = exp(-dot(p, p) / radius);
  float r0 = max(U.x, smallr);
  float r1 = U.x + dens * w;
  vec2  m  = U.yz * (r1 / r0) + kick * w * (p / rr) * r0;
  outColor  = vec4(r1, m, U.w * (r1 / r0));
  outColor1 = passField(kinE(m, max(r1, smallr)) - kinE(U.yz, r0));
}`;

  // The initial field is uniform and vertical, which is exactly divergence-free:
  // the cleaning diagnostic then measures what the scheme does to a clean field
  // rather than what it inherited. |B| = sqrt(2 P / beta), and with rho = cs = 1
  // that is sqrt(2/beta) -- 1.414 at beta = 1, so vA = 1.41 against cs = 1.
  const F_INIT_ISO = HEAD2 + `
uniform float b0;
uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float h1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}
const float TAU = 6.28318530718;
void main(){
  vec2  uv  = vUv;
  float rho = 1.0 + 0.02 * (h1(uvec2(gl_FragCoord.xy), 17u) - 0.5);
  float a = sin(TAU * uv.x * 2.0), b = cos(TAU * uv.y);
  float c = sin(TAU * (uv.x + uv.y)), d = cos(TAU * (uv.x * 2.0 - uv.y));
  vec2  u = 0.8 * vec2(a * b + 0.5 * d, -c + 0.4 * a);
  float Y = 0.04 + 0.03 * sin(TAU * (uv.x * 3.0 + uv.y)) * cos(TAU * (uv.y * 2.0 - uv.x));
  outColor  = vec4(rho, rho * u, rho * Y);
  outColor1 = vec4(0.0, b0, 0.0, 0.0);
}`;

  const F_INIT = HEAD2 + EOS + `
uniform float b0;
uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float h1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}
const float TAU = 6.28318530718;
void main(){
  vec2  uv  = vUv;
  float rho = 1.0 + 0.02 * (h1(uvec2(gl_FragCoord.xy), 17u) - 0.5);
  float a = sin(TAU * uv.x * 2.0), b = cos(TAU * uv.y);
  float c = sin(TAU * (uv.x + uv.y)), d = cos(TAU * (uv.x * 2.0 - uv.y));
  vec2  u = 0.8 * vec2(a * b + 0.5 * d, -c + 0.4 * a);
  float Y = 0.04 + 0.03 * sin(TAU * (uv.x * 3.0 + uv.y)) * cos(TAU * (uv.y * 2.0 - uv.x));
  outColor  = vec4(rho, rho * u, rho * Y);
  // The energy channel, consistent with the pressure the isothermal page has always
  // started from: P = rho cs^2, so the internal energy is rho cs^2/(gamma-1) and is
  // uniform up to the density's own two per cent of noise. This is the seam the unit
  // calibration exists to hide -- see the note on the reference state further down --
  // and it is why the cooling can be switched on without the pressure moving.
  outColor1 = vec4(0.0, b0, 0.0,
                   adia > 0.5 ? etot(rho, u.x, u.y, rho * cs2, b0 * b0) : 0.0);
}`;

  // Re-seed the field alone. Changing the plasma beta, or turning the field off
  // entirely, does not need the gas re-run: each half of the state pair has its own
  // single-attachment framebuffer as well as the shared MRT one, so this writes the
  // field target and leaves the flow going. What it does discard is psi and whatever
  // structure the field had, which is the right thing for an answer to "what would this
  // look like at a different beta".
  //
  // Under a gamma law it discards the thermal state as well, for the same reason and by
  // the same choice: the energy channel lives in the target being rewritten, and a pass
  // cannot read the texture it is writing. So the gas comes back at the isothermal
  // pressure, rho cs^2, which is the state the whole calibration is built around. Move
  // the beta on a two-phase box and it starts its thermal history again -- which is the
  // honest reading of a question that has just changed how much of the pressure was
  // magnetic. The density structure is left alone.
  const F_BSEED_ISO = HEAD + `
uniform float b0;
void main(){ outColor = vec4(0.0, b0, 0.0, 0.0); }`;

  const F_BSEED = HEAD + EOS + `
uniform sampler2D uU;
uniform float b0;
void main(){
  vec4 U = texelFetch(uU, ivec2(gl_FragCoord.xy), 0);
  float r = max(U.x, smallr), ri = 1.0 / r;
  outColor = vec4(0.0, b0, 0.0,
                  adia > 0.5 ? etot(r, U.y * ri, U.z * ri, r * cs2, b0 * b0) : 0.0);
}`;

  // Turning the gamma law on under a running box, which is a change of state and not of
  // a uniform. The energy channel holds whatever the isothermal path left in it -- zero
  // -- and zero read as a total energy is a negative pressure on the very next frame,
  // so it has to be filled in before the first adiabatic step rather than by it. Two
  // ways in, on one switch:
  //
  //   fromE = 0   the pressure is the isothermal one, rho cs^2. This is the seam the
  //               calibration was chosen to make invisible, so switching the cooling
  //               on mid-flight changes nothing about the pressure field at all.
  //   fromE = 1   the pressure comes out of the energy channel under the *previous*
  //               gamma, which is what moving gamma on a live box needs: the gas keeps
  //               the pressure it had and only its heat capacity changes.
  //
  // Both attachments, and the caller swaps. A pass that touched the field target alone
  // would have to read the channel it was rewriting, and a texture cannot be both.
  const F_ENERGIZE = HEAD2 + EOS + `
uniform sampler2D uU, uB;
uniform float fromE, gamPrev;
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);
  float r = max(U.x, smallr), ri = 1.0 / r;
  float u = U.y * ri, v = U.z * ri;
  float b2 = dot(B.xy, B.xy);
  float pg = fromE > 0.5
    ? max((gamPrev - 1.0) * (B.w - 0.5 * r * (u * u + v * v) - 0.5 * b2), smallr * 1e-4)
    : r * cs2;
  outColor  = U;
  outColor1 = vec4(B.xyz, adia > 0.5 ? etot(r, u, v, pg, b2) : 0.0);
}`;

  // The cooling function itself, in cgs, and therefore the one place on the page where
  // the code units have to be pinned to physical ones. They are, once, in UNITS below;
  // the shader receives scaleT2, scaleT and nRef and asks no questions.
  const ISM = `
// ---- the two-phase ISM: a cooling function and the ODE that integrates it -
//
// Ported from ism_net_cooling_h2 and solve_cooling_ism_h2_gpu in
// gpu/gpu_cooling.cuf, reduced to the one case this page runs. The reduction is
// not an approximation of the reference: each clause below makes whole terms
// identically zero or identically one, and they are gone rather than multiplied
// out on every texel of every frame.
//
//   fH2 = 0 exactly. Nothing here forms molecular hydrogen and there is no
//     chemistry carried to evolve it, so h2_line_cooling, the nH2 channel of
//     the C+ 158 um line and the cosmic-ray heating of H2 are identically zero
//     and are dropped. In mu the (1 - fH2/2) factor is 1, and the neutral
//     fraction is simply 1 - xHII.
//   No shielding: fDust = fLW = 1, so G0 is the bare 1.7 Habing everywhere and
//     the Lyman-Werner rate never enters.
//   Dust traces the gas: dust_to_gas = 0.01, so h2_dust_norm is exactly 1 and
//     the photoelectric and grain-recombination dust factors are both 1, folded
//     away. dust_chem_coupling is off.
//   dust_ne_grain_recombination is off, which collapses
//     dust_stage8_composition to h2_composition -- the xHII/mu fixed point
//     below, with no grain-assisted channel in the ionization balance and no
//     27-pass secondary loop.
//
// Everything else in the inventory survives: photoelectric heating with the
// grain-charge efficiency fit, cosmic-ray heating of atomic H, C+ fine
// structure off H and off electrons, O fine structure off H, Lyman-alpha,
// grain-assisted recombination cooling, and the CIE metal term blended in over
// the octave above 10^4 K.
//
// That inventory is the whole point of having it. Lambda changes sign twice
// along an isobar, so the equilibrium pressure P_eq(n) turns over at
// nH = 3.7305 and again at nH = 34.011 and the branch between them is
// isobarically unstable -- dlnP/dlnn = -0.98 at the centre. The cold phase this
// page shows is that instability doing its own work; nothing switches it on.
//
// On precision, because the reference is real(kind=8) throughout and this is
// not. Five things are written differently here. All five are exact algebra --
// no term is clamped, dropped or floored to make float32 reach -- and each one
// is marked where it happens:
//
//   1  Lambda is assembled per nH^2, term by term, instead of forming
//      (cool - heat) and dividing by nH*nH at the end. Every cooling channel is
//      a two-body rate and every heating channel a one-body one, so this is a
//      rearrangement and nothing else. It lifts the smallest live intermediate
//      from ~1e-34 to ~1e-27, and the final division by nH*nH -- a product that
//      is itself zero in float32 below nH = 1.1e-19, so the division is 0/0 and
//      the frame is black -- becomes a single division by nH.
//   2  The ionization balance is divided through by zeta_cr. Written the
//      reference's way the discriminant is built from b*b, which bottoms out at
//      zeta_cr^2 = 4e-32, and 4*a*zeta_cr, which reaches 7e-32 at the thin hot
//      corner and 7e-36 at the density guard: six decades over the smallest
//      normal float32 where the page runs and two at the guard, for quantities
//      that are squares of numbers somebody might change. Scaled, the same
//      quadratic runs 1 .. 1e9 and the question does not come up.
//   3  The ODE carries its integration variable in units of time_max. time_max
//      is ~1e26 for a page-sized step and ~1e35 for a Gyr, and the reference
//      accumulates it by adding 1/wcool, so the increment is lost to rounding
//      as soon as it falls below 1e-7 of the total -- and then the terminal
//      interpolation subtracts two numbers that agree to seven digits and keeps
//      the noise. Normalised, wmax = 1/time_max is the literal 1.0, Lambda and
//      Lambda' appear only through the O(1) products Lambda*time_max, and the
//      time still remaining is carried directly instead of being differenced.
//   4  The CIE exponential is factored, 10^(-21.7 - lc*lc) written as
//      1.995e-22 * 2^(-3.3219 lc*lc). The reference's exponent reaches -57 at
//      the CMB floor, and exp2 of -190 is not a float32 number; factored, the
//      exponential is 1 at the peak and never below 2.5e-36 anywhere T can be.
//   5  The CIE blend and fit variable are logs of a ratio rather than a constant
//      subtracted from a log. This is the one that mattered most: the reference's
//      (log10 T - 4) throws away five leading digits where the blend switches on,
//      and it cost a factor of ten on Lambda as a whole. Details at the term.
//
// And one more, in ismMu, which is not about range at all but about there being
// one rounding fewer in the innermost helper.
//
// The one true clamp is exp(-x) above x = 80, and it is a clamp on nothing. Of
// the three exponentials, 91.2/T never reaches 80 at all (T is floored at
// 2.725 K, so it stops at 33); 228/T reaches it only in the sliver below
// T = 2.85 K; and 1.184e5/T reaches it below T = 1480 K, where the term it
// multiplies is 1e-56 against a Lambda of 1e-26. e^-80 is 1.8e-35 and e^-88 is
// already subnormal, so returning zero past 80 is the value to thirty digits and
// nothing depends on a driver keeping subnormals.
//
// Called as: T2 = T/mu in kelvin is the reference's thermal variable and it is
// exactly what f.w carries once the pressure is a primitive, since
// p = rho kB T/(mu mH) = rho kB T2/mH with no ionization state in it -- which is
// why the reference integrates T2 and not T. So T2 = p mH/(rho kB) and
// nH = X_H rho/mH in cgs, and the step is one call:
//   T2 += solveCoolingIsm(nH, T2, dtSec, gam);
// Isothermal mode has no thermal energy to integrate, so this is called only
// when adia > 0.5.

const float ISM_A_C   = 1.6e-4;          // gas-phase carbon, amr_commons default
const float ISM_A_O   = 3.2e-4;          // gas-phase oxygen
const float ISM_ZCR   = 2.0e-16;         // cosmic-ray ionization rate, 1/s
const float ISM_IZCR  = 5.0e15;          // 1/zeta_cr, folded (note 2 above)
const float ISM_TCMB  = 2.725;           // the floor everything sits on
const float ISM_G0    = 1.7;             // Habing, unattenuated: fDust = 1
const float ISM_Z     = 1.0;             // metallicity, solar
const float ISM_XH    = 0.76;            // hydrogen mass fraction
const float ISM_KB    = 1.380649e-16;    // erg/K
const float ISM_VARMAX = 4.0;            // the reference's step limiter
const int   ISM_ITER_CAP = 40;           // see the note on solveCoolingIsm

// Case-B recombination. pow of a strictly positive base: T is floored at the
// CMB temperature before it gets here, and 1e4 is exact in float32 so the
// division is one correctly-rounded op rather than a rounded reciprocal.
float ismAlphaB(float T){
  return 2.59e-13 * pow(max(T, ISM_TCMB) / 1.0e4, -0.7);
}

// The HII fraction from ionization balance, scaled by zeta_cr. The reference is
//   xHII = 2 z rem / (b + sqrt(b*b + 4 a z rem)),  a = alphaB nH,  b = a A_C + z
// with z = zeta_cr = 2e-16 and rem = 1 - fH2 = 1 here. Dividing numerator and
// denominator by z is exact and moves the whole quadratic up thirty decades:
//   as = a/z,  bs = as A_C + 1,  xHII = 2/(bs + sqrt(bs*bs + 4 as))
// The reference's rationalised form is kept as it stands -- it is the branch of
// the quadratic that has no cancellation in it, which matters far more here than
// it does in double. xHII <= 1 follows from bs >= 1, so the neutral fraction
// comes out non-negative without a guard.
float ismXhii(float nH, float T){
  float as = ismAlphaB(T) * max(nH, 1.0e-30) * ISM_IZCR;
  float bs = as * ISM_A_C + 1.0;
  return 2.0 / (bs + sqrt(bs * bs + 4.0 * as));
}

// Mean molecular weight. fH2 = 0, so the reference's (1 - fH2/2) is 1, and then
// its X (1 + yHe + xe) with yHe = (1 - X)/(4X) is X + (1 - X)/4 + X xe, and for
// X = 0.76 the constant part is exactly 0.82. Same number, one rounding fewer in
// the function that every other one of these calls five times over.
float ismMu(float xe){
  return 1.0 / (0.82 + ISM_XH * xe);
}

// The composition closure of h2_composition: xHII, xe, mu and T from nH and
// T2 = T/mu, by fixed point. mu depends on the free electrons and the ionization
// balance depends on T = mu T2, so the two are solved together; four passes of
// ismXhii, which is the reference's count and worth keeping to the letter,
// because the equilibrium curve the page is initialised on was tabulated with
// it. It is converged long before that -- xe only ever moves mu by a part in
// 1e2 -- but "long before" is not a reason to return a different number.
//
// mu is what a caller needs to turn T2 back into a temperature to colour by.
struct IsmComp { float xHII; float xe; float mu; float T; };

IsmComp ismComposition(float nH, float T2){
  float xHII = ismXhii(nH, max(ismMu(ISM_A_C) * T2, ISM_TCMB));
  float xe = ISM_A_C + xHII;
  float mu = ismMu(xe);
  float T  = max(mu * T2, ISM_TCMB);
  for (int i = 0; i < 3; ++i) {
    xHII = ismXhii(nH, T);
    xe   = ISM_A_C + xHII;
    mu   = ismMu(xe);
    T    = max(mu * T2, ISM_TCMB);
  }
  return IsmComp(xHII, xe, mu, T);
}

// exp(-x) for x >= 0, with the tail cut where float32 stops having normals.
float ismExpNeg(float x){
  return x > 80.0 ? 0.0 : exp(-x);
}

// Lambda(nH, T2) in erg cm^3 / s. Positive is net cooling and the sign change is
// the equilibrium curve, so this function's zero set is the physics: everything
// the page does thermally is where cool - heat crosses.
//
// Per note 1, the cooling terms have shed their two powers of nH and the heating
// terms carry one in the denominator, which is why nH appears exactly once,
// under a single division. That division is the one thing here that has to be
// guarded, and 1e-6 is the value that makes the guard do two jobs: xe is never
// below A_C = 1.6e-4, so at nH >= 1e-6 the electron density is never below
// 1.6e-10 and the reference's own ne = max(xe nH, 1e-30) floor is unreachable --
// which is what lets ne be carried as the fraction xe and not as a density.
// 1e-6 cm^-3 is eight decades under anything this page can make.
float ismNetCooling(float nH, float T2){
  float n = max(nH, 1.0e-6);
  IsmComp c = ismComposition(n, T2);
  float T   = c.T;
  float xe  = c.xe;
  float xHI = 1.0 - c.xHII;            // fH2 = 0: all the rest of the H is neutral
  float t100 = T / 100.0;

  // the grain charge parameter G0 sqrt(T)/ne and the two-branch photoelectric
  // efficiency it feeds. The 1e8 ceiling is the reference's; over the range this
  // page runs -- nH >= 1e-2, T <= 1.3e6 -- it never binds, the largest charge
  // reached being 4e6 at the hot, thin corner.
  float gc = min(ISM_G0 * sqrt(T) / (xe * n), 1.0e8);
  // t100^0.7 serves twice: (T/1e4)^0.7 is the same power up to 100^-0.7, and
  // 3.65e-2 * 100^-0.7 = 1.4530911725e-3 is folded. One pow, not two.
  float p07 = pow(t100, 0.7);
  float epspe = 4.87e-2 / (1.0 + 4.0e-3 * pow(gc, 0.73))
              + 1.4530911725e-3 * p07 / (1.0 + 2.0e-2 * gc);

  // heating: grains, then cosmic rays on atomic H. Both are one-body.
  float heat = (1.3e-24 * epspe * ISM_G0 + ISM_ZCR * 1.6e-11 * xHI) / n;

  // C+ 158 um, excited by neutral H and by electrons. t100^(-0.5) is
  // inversesqrt, which is the same number computed without a log and an exp.
  float kHcii = 7.6e-10 * pow(t100, 0.14);
  float kecii = 3.0e-7 * inversesqrt(t100);
  float cool = ISM_A_C * ISM_Z * 1.26e-14 * 2.0 * ismExpNeg(91.2 / T)
             * (xHI * kHcii + xe * kecii);

  // O 63 um off neutral H, and Lyman-alpha off electrons. The Lyman-alpha
  // exponential is the one that falls off the bottom of float32, at T = 1480 K,
  // by which point the term is 1e-56 against a Lambda of 1e-26.
  float kHoi = 4.0e-11 * p07;
  cool += ISM_A_O * ISM_Z * 3.14e-14 * 0.6 * ismExpNeg(228.0 / T) * xHI * kHoi;
  cool += 7.3e-19 * xe * xHI * ismExpNeg(1.184e5 / T);

  // recombination of electrons onto grains: the charge parameter again, this
  // time with a temperature-dependent index.
  float beta = 0.74 / pow(T, 0.068);
  cool += 4.65e-30 * pow(T, 0.94) * pow(gc, beta) * xe;

  // the CIE metal term, blended in across the octave above 10^4 K so that the
  // fit is not asked to cover gas it was never meant to. Note 4 is here, and so
  // is note 5, which is the largest float32 effect in this whole function:
  //
  // The reference's blend is (log10 T - 4)/log10 2 and its fit variable is
  // (log10 T - 5.2)/0.8. Both subtract a constant from a logarithm, and where the
  // blend is switching on log10 T is 4.08 -- so (log10 T - 4) throws away the
  // leading five digits and hands the rest of the term a relative error of 3e-6
  // where float32 had 6e-8. The CIE term is 86% of the cooling rate there, so
  // Lambda inherits it: 1.1e-5, two orders worse than anything else here.
  //
  // Both are logs of the same ratio, so take the log of the ratio instead.
  // q = log2(T/1e4) is the blend outright, and
  // (log10 T - 5.2)/0.8 = q log10(2)/0.8 - 1.2/0.8, exactly. One log2 as before,
  // no cancellation in either, and the error goes back to being T's own.
  float q     = log2(T / 1.0e4);
  float blend = clamp(q, 0.0, 1.0);
  float lc    = q * 0.37628749457997650 - 1.5;   // (log10 T - 5.2)/0.8
  cool += 1.9952623149688828e-22 * ISM_Z * blend * exp2(-3.3219280948873622 * lc * lc);

  return cool - heat;
}

// The CMB floor in the solver's own variable: h2_t2_floor, which is T_CMB/mu
// evaluated at the electron fraction the floor itself implies. Not T_CMB -- mu
// is 1.22, so the floor on T2 is 2.24 K and a solver that used 2.725 would sit
// half a kelvin above the radiation field forever.
float ismT2Floor(float nH){
  return ISM_TCMB / ismMu(ISM_A_C + ismXhii(nH, ISM_TCMB));
}

// The thermal ODE: returns the change in T2 over dtSec, which is what the caller
// adds. solve_cooling_ism_h2_gpu with fH2 = 0, algorithm for algorithm.
//
// The reference integrates dtau/dtime = -Lambda(tau) in a time variable that
// carries the density and the units, time_max = dt (gamma-1) X_H nH/kB, and
// takes linearised-implicit steps
//
//   tau <- tau (1 + L'/w - L/(tau w)) / (1 + L'/w)
//
// which is one Newton step on the backward-Euler equation with the step size
// itself, 1/w, chosen by a three-way limiter: w = max(|L|/tau varmax, wmax,
// -L' varmax) with varmax = 4 and wmax = 1/time_max. Each of the three earns its
// place. The first bounds the relative change per step -- |L|/tau <= w/4 makes
// the correction at most a quarter, so tau moves by no more than a third and
// cannot go negative or overflow. The second stops the step from overshooting dt.
// The third is what makes it safe where -L' is large, and incidentally what
// makes 1 + L'/w >= 3/4 rather than something that could vanish: the branch
// where L' < 0 is exactly the branch that limiter catches, so the denominator
// needs no guard.
//
// Per note 3 the time runs in units of time_max, so wmax is the literal 1.0,
// w >= 1 means no step can be longer than what is left, and Lambda enters only
// as Lh = Lambda time_max -- a product of 1e-26 and 1e26, which is the point.
//
// The cap: the Fortran's is 500, which is about 35 transcendentals x 2 rate
// evaluations x 500 in a fragment shader and not a serious proposition. Over the
// range this page sees -- nH in [1e-2, 1e3], T2 in [10, 1e6], dt anywhere from
// 1e-3 to 1e2 cooling times -- the measured trip count is 1 at the median, 15 at
// the 90th percentile, 25 at the 99th and 31 at the worst; 32 is the smallest cap
// that reproduces the uncapped answer bit for bit over that whole sweep, and 40
// is what is set here, for margin. It is a parameter because that measurement is
// a property of the range and not of the algorithm: the trip count is set by how
// many factors of 4/3 in tau the step has to cross, so a page that let the gas
// swing further would need a larger one.
//
// One deliberate divergence, on the last line before the return. The reference
// interpolates back to time_max unconditionally, and when it leaves by the
// iteration cap rather than by finishing, time < time_max and that
// "interpolation" is an extrapolation with a weight of 1/(1 - t/time_max) -- it
// can be 90, and tau = 90 tau_new - 89 tau_old is not an answer, it is a black
// frame. min(w, 1.0) makes a cap-out stop early instead, which is wrong by a
// bounded amount rather than unbounded. It cannot fire unless the cap fires, so
// in every converged case this is the reference to the bit.
float solveCoolingIsm(float nH, float T2, float dtSec, float gamma, int iterCap){
  if (dtSec <= 0.0) return 0.0;
  float n = max(nH, 1.0e-6);
  float floorT2 = ismT2Floor(n);
  float tauIni = T2;
  float tau = max(tauIni, floorT2);
  float tauOld = tau;

  // dt (gamma-1) X_H nH / kB. 1/kB is 7.2e15, so this is 1e26 for a page-sized
  // step -- and it overflows float32 at dt nH = 9e22, which is a step of 3e14
  // years at nH = 1 and nothing this page can ask for. It is bounded anyway,
  // because the overflow is not the harm: inf reaches W, ds becomes 0, and the
  // update evaluates inf*0, which is the one NaN in here that could get out. min
  // takes the finite side of an inf, so one min covers both this product and the
  // dtSec*n inside it, and at the bound the integration is already 1e11 cooling
  // times long -- a decade either way changes nothing about where it lands.
  float timeMax = min(dtSec * n * ((gamma - 1.0) * ISM_XH / ISM_KB), 3.0e37);
  float rem = 1.0;          // time left, in units of time_max
  float remOld = 1.0;
  float ds = 0.0;           // length of the last step, same units

  for (int i = 0; i < iterCap; ++i) {
    if (rem <= 0.0) break;
    float L  = ismNetCooling(n, tau);
    // the reference's one-sided difference, at a 1e-3 relative offset. This is
    // the one place the float32 cancellation is not fixable: the two rates agree
    // to 1e-3 by construction and each is good to ~1e-6, so L' carries 2e-4 at
    // the median and 2% at the worst. It is affordable because L' appears only
    // as a Newton correction and inside a max() -- degrade it and the step gets
    // smaller and the answer does not move, which is exactly what the wcool
    // limiter is for. Widening the offset trades this for truncation error and
    // was not worth it: at 1e-2 the limiter starts to see the curvature.
    float Lp = (ismNetCooling(n, tau * (1.0 + 1.0e-3)) - L) / (tau * 1.0e-3);
    float Lh  = L  * timeMax;
    float Lph = Lp * timeMax;
    float W = max(max(abs(Lh) / tau * ISM_VARMAX, 1.0), -Lph * ISM_VARMAX);
    ds = 1.0 / W;
    float g = Lph * ds;                                  // L'/wcool
    tauOld = tau;
    tau = tau * (1.0 + g - Lh * ds / tau) / (1.0 + g);
    tau = max(tau, floorT2);
    remOld = rem;
    rem -= ds;
  }

  // land on time_max exactly: linear interpolation back across the last step,
  // its weight being (time_max - time_old)/(time - time_old) with the two
  // subtractions already done for us.
  if (ds > 0.0) {
    float w = min(remOld / ds, 1.0);
    tau = tau * w + tauOld * (1.0 - w);
  }
  return max(tau, floorT2) - tauIni;
}

float solveCoolingIsm(float nH, float T2, float dtSec, float gamma){
  return solveCoolingIsm(nH, T2, dtSec, gamma, ISM_ITER_CAP);
}
`;

  // ---- cooling -------------------------------------------------------------
  //
  // One operator-split pass, after the hydrodynamic update, at the full step, once per
  // step -- which is where amr_step puts it: godunov_fine, then cooling_fine, and no
  // subcycling between them. The gas does not know it is being cooled while it is being
  // advected and it does not have to: the two operators commute to first order in dt and
  // dt here is a few thousandths of a cooling time.
  //
  // Everything the pass has to do is a consequence of the energy being conserved rather
  // than the temperature being carried. cooling_fine.f90 takes the kinetic energy out of
  // the conserved total to get the thermal part, and takes the magnetic energy out as
  // well when there is a field; both come out here, in that order, and go back
  // afterwards untouched, because cooling changes the pressure and nothing else. The
  // density, the momentum and the field are written back bit for bit.
  //
  // Both attachments and a swap, for the reason F_ENERGIZE has: the channel being
  // rewritten is in the texture being read.
  const F_COOL = HEAD2 + EOS + DTDX + ISM + `
uniform sampler2D uU, uB;
uniform float dxCell, scaleT2, scaleT, nRef;

void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);
  outColor = U;
  if (adia < 0.5) { outColor1 = B; return; }

  float r = max(U.x, smallr), ri = 1.0 / r;
  float u = U.y * ri, v = U.z * ri;
  float b2 = dot(B.xy, B.xy);
  // the thermal energy, which is what is left of the conserved total once the kinetic
  // and the magnetic have been taken out of it, as a pressure
  float pg = pfromE(B.w, r, u, v, b2);

  // and the two physical numbers the cooling function wants. T2 = T/mu is the
  // reference's own thermal variable and it comes straight out of the code pressure,
  // because p = rho kB T/(mu mH) = rho kB T2/mH has no ionization state in it.
  float T2 = pg * ri * scaleT2;
  float nH = r * nRef;
  // dt in seconds. The same stepDtDx() the hydrodynamic pass used -- the Courant
  // texture has not been rewritten in between -- so the split is over one step and not
  // over one and a bit.
  float dtSec = stepDtDx() * dxCell * scaleT;

  float T2n = T2 + solveCoolingIsm(nH, T2, dtSec, gam);
  // A non-finite excursion in one texel would poison the whole frame through the next
  // step's stencil, so it is refused here rather than diagnosed later. Nothing measured
  // reaches this, which is the point of it being a single comparison.
  if (!(T2n > 0.0)) T2n = T2;

  float pn = T2n * r / scaleT2;
  outColor1 = vec4(B.xyz, etot(r, u, v, pn, b2));
}`;

  // What the dust samples: velocity and field in one filterable texture, so a
  // grain gathers both in four bilinear taps instead of eight. RGBA16F, because
  // asking for LINEAR on RGBA32F makes the texture incomplete and every read --
  // texelFetch included -- silently returns zero.
  const F_FLOW = HEAD + `
uniform sampler2D uU, uB;
uniform float smallr;
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);
  outColor = vec4(U.yz / max(U.x, smallr), B.xy);
}`;

  // Everything the picture needs, at grid resolution, so the composite can afford
  // a proper reconstruction: the dye, the convergence of the velocity field, and -- for
  // the multiphase view -- the temperature and the density.
  //
  // The temperature is stored as log2 T in kelvin, because that is what a phase map is a
  // function of and because a half float holds three decades of it with a resolution of
  // a thousandth. mu is taken at its reference value rather than solved for: across the
  // whole two-phase range it moves between 1.193 and 1.219, two per cent, which is a
  // hundredth of a decade in log T and nothing at all to a colour. Solving for it would
  // mean the composition fixed point, four pow calls, in a pass that exists to feed the
  // display. Isothermal mode returns T_ref in every cell, which is the honest answer: a
  // gas that cannot change its temperature has one.
  const F_DISP = HEAD + EOS + `
uniform sampler2D uU, uB;
uniform ivec2 size;
uniform float tScale;
vec4 cell(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uU, p, 0);
}
vec3 q(ivec2 c){ vec4 U = cell(c); float ri = 1.0 / max(U.x, smallr); return vec3(U.w * ri, U.yz * ri); }
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec3 q0 = q(c);
  vec3 qe = q(c + ivec2(1, 0)), qw = q(c - ivec2(1, 0));
  vec3 qn = q(c + ivec2(0, 1)), qs = q(c - ivec2(0, 1));
  float div = 0.5 * ((qe.y - qw.y) + (qn.z - qs.z));

  vec4 U = cell(c);
  vec4 B = texelFetch(uB, c, 0);
  float r = max(U.x, smallr), ri = 1.0 / r;
  float pg = adia > 0.5 ? pfromE(B.w, r, U.y * ri, U.z * ri, dot(B.xy, B.xy)) : r * cs2;
  float T = max(pg * ri * tScale, 1.0);
  outColor = vec4(q0.x, max(-div, 0.0), log2(T), r);
}`;


  // ---- the display field, smoothed ----------------------------------------
  //
  // Why this pass exists. The convergence channel is a centred difference of the
  // velocity, so a captured front is one cell wide and effectively all of its power
  // sits at the grid's Nyquist frequency. No interpolant can fix that: Catmull-Rom
  // is C1 and passes through the samples, which is exactly right for the dye, but
  // fed a one-cell spike it has no sub-cell information to work with and reproduces
  // the lattice -- lozenges along cell boundaries, and ringing from its negative
  // lobes. Stretched over eight screen pixels per cell that reads as grid artifacts
  // on every shock, which is what it was.
  //
  // So the smoothing happens before the reconstruction, not in it: a 5x5 binomial on
  // the convergence channel at grid resolution, which turns the spike into a
  // resolved two-and-a-half-cell feature the interpolant can actually represent. It
  // costs one pass over seventeen thousand pixels, against doing anything at all at
  // display resolution, where there are seventy times as many. The dye is passed
  // through untouched -- it is a smooth field already and blurring it would cost the
  // wisps their contrast.
  const F_DSMOOTH = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 size;
uniform float sigma;
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  float s2 = 2.0 * max(sigma * sigma, 1e-6);
  float acc = 0.0, wsum = 0.0;
  for (int j = -2; j <= 2; j++) {
    for (int i = -2; i <= 2; i++) {
      ivec2 p = c + ivec2(i, j);
      p = ivec2((p.x + size.x) % size.x, (p.y + size.y) % size.y);
      float w = exp(-float(i * i + j * j) / s2);
      acc += w * texelFetch(uSrc, p, 0).y;
      wsum += w;
    }
  }
  outColor = vec4(texelFetch(uSrc, c, 0).x, acc / wsum, 0.0, 1.0);
}`;

  // ---- reconstruction -----------------------------------------------------

  // Catmull-Rom for the picture: it interpolates, so a front stays where the
  // solver put it, and it is C1 across cell boundaries, which is what kills the
  // facets bilinear leaves when a 128-cell grid is stretched over 800 pixels.
  // Sixteen point samples in nine hardware taps.
  const CATROM = `
vec4 texCR(sampler2D tex, vec2 uv, vec2 size){
  vec2 p  = uv * size;
  vec2 c1 = floor(p - 0.5) + 0.5;
  vec2 f  = p - c1;
  vec2 f2 = f * f, f3 = f2 * f;
  vec2 w0 = -0.5 * f3 + f2 - 0.5 * f;
  vec2 w1 =  1.5 * f3 - 2.5 * f2 + 1.0;
  vec2 w2 = -1.5 * f3 + 2.0 * f2 + 0.5 * f;
  vec2 w3 =  0.5 * f3 - 0.5 * f2;
  vec2 w12 = w1 + w2;
  vec2 t0 = (c1 - 1.0) / size;
  vec2 t3 = (c1 + 2.0) / size;
  vec2 t12 = (c1 + w2 / w12) / size;
  return (texture(tex, vec2(t0.x,  t0.y))  * w0.x +
          texture(tex, vec2(t12.x, t0.y))  * w12.x +
          texture(tex, vec2(t3.x,  t0.y))  * w3.x) * w0.y
       + (texture(tex, vec2(t0.x,  t12.y)) * w0.x +
          texture(tex, vec2(t12.x, t12.y)) * w12.x +
          texture(tex, vec2(t3.x,  t12.y)) * w3.x) * w12.y
       + (texture(tex, vec2(t0.x,  t3.y))  * w0.x +
          texture(tex, vec2(t12.x, t3.y))  * w12.x +
          texture(tex, vec2(t3.x,  t3.y))  * w3.x) * w3.y;
}
`;

  // A cubic B-spline for the display field, alongside the Catmull-Rom above --
  // because the two channels want opposite things.
  //
  // Catmull-Rom interpolates, which is what the dye wants: it passes through the
  // samples, so a wisp keeps its amplitude. But it is built from a kernel with
  // negative lobes, and negative lobes on a one-cell feature ring: an undershoot
  // ridge on each side, aligned with the axes because the kernel is separable.
  // That is what the facets on the fronts actually were -- not the interpolation
  // being too coarse, but it overshooting a feature it could not represent.
  //
  // The B-spline is C2 and strictly positive, so it cannot ring at all. It
  // approximates instead of interpolating, which costs a front about half a cell of
  // width, and that is exactly the trade to take here: nobody is reading a
  // convergence value off the screen. Four bilinear taps against nine.
  const BSPLINE4 = `
vec4 texBS(sampler2D tex, vec2 uv, vec2 size){
  vec2 p = uv * size - 0.5;
  vec2 i = floor(p), f = p - i;
  vec2 f2 = f * f, f3 = f2 * f;
  vec2 w0 = (1.0 - 3.0 * f + 3.0 * f2 - f3) * (1.0 / 6.0);
  vec2 w1 = (4.0 - 6.0 * f2 + 3.0 * f3) * (1.0 / 6.0);
  vec2 w2 = (1.0 + 3.0 * f + 3.0 * f2 - 3.0 * f3) * (1.0 / 6.0);
  vec2 w3 = f3 * (1.0 / 6.0);
  vec2 g0 = w0 + w1, g1 = w2 + w3;
  vec2 h0 = (i - 0.5 + w1 / g0) / size;
  vec2 h1 = (i + 1.5 + w3 / g1) / size;
  return g0.y * (g0.x * texture(tex, vec2(h0.x, h0.y)) +
                 g1.x * texture(tex, vec2(h1.x, h0.y)))
       + g1.y * (g0.x * texture(tex, vec2(h0.x, h1.y)) +
                 g1.x * texture(tex, vec2(h1.x, h1.y)));
}
`;

  // A cubic B-spline for the grains: the shape function a particle-in-cell code
  // gathers with, and for the same reason -- it is C2, so a grain crossing a cell
  // boundary feels no kink in its acceleration and the population cannot print the
  // grid onto its own clustering. Four bilinear taps, and here they carry the
  // field as well as the velocity.
  const BSPLINE = `
vec4 gather4(sampler2D tex, vec2 uv, vec2 size){
  vec2 p = uv * size - 0.5;
  vec2 i = floor(p), f = p - i;
  vec2 f2 = f * f, f3 = f2 * f;
  vec2 w0 = (1.0 - 3.0 * f + 3.0 * f2 - f3) * (1.0 / 6.0);
  vec2 w1 = (4.0 - 6.0 * f2 + 3.0 * f3) * (1.0 / 6.0);
  vec2 w2 = (1.0 + 3.0 * f + 3.0 * f2 - 3.0 * f3) * (1.0 / 6.0);
  vec2 w3 = f3 * (1.0 / 6.0);
  vec2 g0 = w0 + w1, g1 = w2 + w3;
  vec2 h0 = (i - 0.5 + w1 / g0) / size;
  vec2 h1 = (i + 1.5 + w3 / g1) / size;
  return g0.y * (g0.x * texture(tex, vec2(h0.x, h0.y)) +
                 g1.x * texture(tex, vec2(h1.x, h0.y)))
       + g1.y * (g0.x * texture(tex, vec2(h0.x, h1.y)) +
                 g1.x * texture(tex, vec2(h1.x, h1.y)));
}
`;

  // ---- dust ---------------------------------------------------------------

  const SPECTRUM = `
uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float hash1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}
float stRank(ivec2 tc){ return hash1(uvec2(tc), 7u); }
float tauMul(float s){ return exp((s * 2.0 - 1.0) * 1.6118); }
float aRel(float s){ return exp((s - 1.0) * 2.0 * 1.6118); }
`;

  // compute_lorentz, from mini-ramses pm/move_fine.f90, with Bz = 0 substituted
  // and the sign of dteff as the reference carries it. The matrix is the Cayley
  // transform of the cross-product matrix,
  //
  //   w' = Inverse[e - (dteff/2) m] . [e + (dteff/2) m] . w
  //
  // which is an exact rotation about B through 2 atan(|B| dteff / 2) -- so |w| is
  // preserved for any step and the drift along B is untouched, both of which are
  // visible in the algebra below: the .y row divides by det where the diagonal is
  // det. That is what makes this safe at charge 100, where an explicit push would
  // need the gyration substepped.
  const LORENTZ = `
vec3 lorentz(vec3 w, vec2 B, float dt, float q){
  float d  = -dt * q;                       // dteff
  float d2 = d * d;
  float bx = B.x, by = B.y;
  float b2 = bx * bx + by * by;
  float det = 1.0 + 0.25 * b2 * d2;
  float k = 0.25 * (bx * bx - by * by) * d2;
  float s = 0.5 * d2 * bx * by;
  vec3 o;
  o.x = (1.0 + k) * w.x + s * w.y           - d * by * w.z;
  o.y = s * w.x         + (1.0 - k) * w.y   + d * bx * w.z;
  o.z = d * by * w.x    - d * bx * w.y      + (1.0 - 0.25 * b2 * d2) * w.z;
  return o / det;
}
`;

  // The grain push: gather, rotate, drag, move. Two targets, because a charged
  // grain in a planar field does not stay planar.
  const F_PART = HEAD2 + SPECTRUM + BSPLINE + LORENTZ + `
uniform sampler2D uPart, uPartZ, uFlow;
uniform vec2  gridSize;
uniform float dragA, dt, aspect, brown, charge, stream;
uniform uint  uFrame;

void main(){
  ivec2 id  = ivec2(gl_FragCoord.xy);
  vec4  P   = texture(uPart, vUv);
  vec4  Z   = texture(uPartZ, vUv);
  vec2  pos = P.xy;
  vec3  v   = vec3(P.zw, Z.x);

  vec4 g = gather4(uFlow, pos, gridSize);
  vec3 u = vec3(g.xy, 0.0);
  vec2 B = g.zw;

  // In ideal MHD the lab-frame electric field is -u x B, so the whole Lorentz
  // force on a grain is q (v - u) x B: it acts on the drift and nothing else.
  vec3 w = v - u;
  w = lorentz(w, B, dt, charge);

  // then one backward Euler drag stage, on a per-grain stopping time. tau varies
  // across the population but not in time, so the median grain's dt/tau is the
  // only scalar the CPU supplies and a grain rescales it by its place in the size
  // spectrum.
  float a = dragA / tauMul(stRank(id));
  w = w / (1.0 + a);

  if (brown > 0.0) {
    w.xy += brown * (vec2(hash1(uvec2(id), uFrame + 11u),
                          hash1(uvec2(id), uFrame + 523u)) - 0.5);
  }
  v = u + w;

  // bulk streaming, along the local field now that there is one to follow
  vec2 sdir = stream * B / max(length(B), 1e-4);
  vec3 drift = w + vec3(sdir, 0.0);

  pos = fract(pos + vec2((v.x + sdir.x) / aspect, v.y + sdir.y) * dt);
  outColor  = vec4(pos, v.xy);
  outColor1 = vec4(v.z, length(drift), 0.0, 1.0);
}`;

  // The drift scale, measured rather than predicted. On the live site the
  // lighting was normalised by tau * M^2 -- the drift a grain of stopping time tau
  // reaches in a flow of rms Mach M -- and that estimate is simply wrong once the
  // grains are magnetised: the cross-field drift is suppressed by 1/(1 + (omega
  // tau)^2), which at charge 100 is a factor of two hundred, so every grain would
  // sit at the floor of the ramp. The population's own mean drift is reduced on the
  // GPU and read straight out of the vertex shader, so nothing has to be guessed
  // and nothing has to come back to the CPU.
  const F_PDRIFT = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
uniform float first;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 8;
  float sum = 0.0, mx = 0.0;
  for (int j = 0; j < 8; j++) {
    for (int i = 0; i < 8; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      vec4 v = texelFetch(uSrc, p, 0);
      // The first pass reads the particle state, where the drift is .y; later
      // passes read their own output, where the running sum is .x. The maximum
      // lives in .y in both layouts, so it needs no branch.
      sum += first > 0.5 ? v.y : v.x;
      mx = max(mx, v.y);
    }
  }
  outColor = vec4(sum, mx, 0.0, 1.0);
}`;

  const V_PARTICLE = `#version 300 es
precision highp float;
uniform sampler2D uPart, uPartZ, uDrift;
uniform int   uPW;
uniform float uPointSize, uDriftRef, uInvN, uAlpha, uVignette, uAspect;
out float vBright;
out float vDrift;
` + SPECTRUM + `
void main(){
  int   id  = gl_VertexID;
  ivec2 tc  = ivec2(id % uPW, id / uPW);
  vec2  pos = texelFetch(uPart, tc, 0).xy;
  float w   = texelFetch(uPartZ, tc, 0).y;
  float st  = stRank(tc);

  // Drift, not speed: a grain swept along with the gas is doing nothing
  // interesting however fast it is going. The drift is the one the grain actually
  // felt, computed in the push from the gather the drag used, and it includes the
  // out-of-plane component the field put there.
  float mean = texelFetch(uDrift, ivec2(0, 0), 0).x * uInvN;
  vDrift = clamp(w / max(uDriftRef * mean, 1e-5), 0.0, 1.0);

  gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);

  vec2  q  = (pos - 0.5) * vec2(uAspect, 1.0);
  float vg = 1.0 - uVignette * smoothstep(0.30, 1.05, length(q));
  vBright  = uAlpha * vg * (0.30 + 0.70 * vDrift);
  gl_PointSize = max(1.0, uPointSize * aRel(st));
}`;

  const F_PARTICLE = `#version 300 es
precision highp float;
in float vBright;
in float vDrift;
out vec4 outColor;
uniform vec3 uColor;
void main(){
  vec2  d  = gl_PointCoord * 2.0 - 1.0;
  float r2 = dot(d, d);
  if (r2 > 1.0) discard;
  float a = exp(-r2 * 2.7) * (1.0 - r2);
  vec3  c = mix(uColor, mix(uColor, vec3(1.0), 0.42), vDrift * vDrift);
  float b = a * vBright;
  outColor = vec4(c * b, b);
}`;

  // ---- the field, drawn ---------------------------------------------------
  //
  // Line integral convolution: walk a few cells along the field line through both
  // ends of each pixel, average a static noise field over the walk, and what comes
  // back is a texture whose grain lies along B. It is the standard way to draw a
  // vector field without drawing arrows on it, and it is the right one here for a
  // reason that is not only aesthetic: the dust at charge 100 is already stuck to
  // field lines, so the streaks it produces run the same way the grains do. The
  // field reinforces the thing the dust is telling you instead of annotating it.
  //
  // Three choices that keep it inside the palette rather than on top of it:
  //
  //   - the contribution is *signed* around zero, so the field grains the medium
  //     rather than lighting it, and the mean brightness of the frame does not move
  //   - it is drawn in the chapter's own accent, like the dye and the fronts. No
  //     hue appears here that appears nowhere else
  //   - it is scaled by |B| / B0, so a stretched or compressed field shows itself
  //     and an undisturbed one stays quiet
  //
  // The kernel travels along the walk, which makes bright packets run along the
  // field lines at a fixed rate. That is not decoration either: it is what an Alfven
  // wave does, and it is the only motion in the picture that is not advective.
  //
  // Run at half the display resolution into one RG16F, because the pattern is smooth
  // by construction and eleven taps per pixel at full resolution is not free.
  const F_LIC = HEAD + `
uniform sampler2D uFlow, uNoise;
uniform float aspect, time, b0, stepUv, noiseFreq, waveLen, waveSpeed;

void main(){
  vec2 uv = vUv;
  vec4 g0 = texture(uFlow, uv);
  float bmag = length(g0.zw);
  if (bmag < 1e-5) { outColor = vec4(0.0, 0.0, 0.0, 1.0); return; }

  // Both the noise grain and the walk step are set in pixels of *this* target, not
  // in box units, and that is not a detail -- it is the whole thing working or not.
  // Written in box units first, the noise repeated a hundred and ten times across
  // the box, which put its texels thirty-five times below the pixel Nyquist: every
  // sample along the walk drew an independent random number, the average went to a
  // flat half everywhere, and there was no streak to see. The step has the same
  // requirement from the other side -- step further than the noise is correlated and
  // consecutive samples are again independent. So one noise texel is a few target
  // pixels, and the step is the same few pixels.
  const int N = 12;                      // 2N+1 samples along the field line
  float h = stepUv;                      // step, in uv, sized from the target
  float acc = 0.0, wsum = 0.0, w2sum = 0.0;

  // Both directions from the pixel, following the local field as it turns. The
  // walk re-reads B at every step, so a curved field line stays followed rather
  // than being approximated by its tangent at the centre.
  for (int side = 0; side < 2; side++) {
    vec2 p = uv;
    float sgn = side == 0 ? 1.0 : -1.0;
    vec2 b = g0.zw;
    for (int i = 0; i <= N; i++) {
      if (i > 0) {
        b = texture(uFlow, p).zw;
        float m = length(b);
        if (m < 1e-5) break;
        b /= m;
        // physical -> uv: x spans aspect box widths, y spans one
        p += sgn * h * vec2(b.x / aspect, b.y);
      }
      // arc length from the centre, and a kernel that travels along it
      float s = float(i) * h;
      float w = 0.55 + 0.45 * cos(6.2831853 * (s / waveLen - time * waveSpeed));
      if (side == 1 && i == 0) continue;  // the centre sample belongs to one side only
      float n = texture(uNoise, vec2(p.x * aspect, p.y) * noiseFreq).x;
      acc += n * w;
      wsum += w;
      w2sum += w * w;
    }
  }
  // Normalised to unit variance, computed rather than tuned. Averaging n samples of
  // white noise leaves a standard deviation of sigma/sqrt(n_eff), where sigma is
  // 1/sqrt(12) for the uniform bytes in the noise texture and n_eff = wsum^2/w2sum is
  // the kernel's effective sample count. Dividing that out makes the gain in the
  // composite a perceptual knob, and keeps it one when the walk length or the kernel
  // shape changes -- which the hand-set constant that used to be here could not do.
  float lic = wsum > 0.0 ? acc / wsum : 0.5;
  float neff = w2sum > 0.0 ? wsum * wsum / w2sum : 1.0;
  outColor = vec4((lic - 0.5) * sqrt(neff * 12.0), bmag / b0, 0.0, 1.0);
}`;

  // ---- composite ----------------------------------------------------------
  //
  // Otherwise identical to the live site's, deliberately: this page is a change of
  // physics, not a change of picture.
  const F_COMPOSITE = HEAD + CATROM + BSPLINE4 + `
uniform sampler2D uDisp, uDispS, uLic;
uniform float uFieldGain, uFieldFloor;
uniform vec2  gridSize;
uniform vec2  uRes;
uniform vec3  uBg, uTint, uAccent, uCold, uWarm;
uniform float uTime, uGrain, uDyeGain, uVignette, uShockGain, uShockLift, uShockPow;
uniform float uPhase, uPhaseMode;

uint uhash(uint x){
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float h1(uvec2 p, uint s){
  return float(uhash(p.x * 1973u + p.y * 9277u + s * 26699u)) * (1.0 / 4294967296.0);
}

void main(){
  // The dye, interpolated with Catmull-Rom off the raw field so a wisp keeps its
  // amplitude. The convergence comes from a different texture through a different
  // interpolant, below, and only if it is being drawn at all.
  vec4  d0   = texCR(uDisp, vUv, gridSize);
  float dye0 = d0.x;

  float dens = clamp(dye0 * uDyeGain, 0.0, 1.0);
  vec3  dye  = uAccent * dens * dens;

  float g = smoothstep(1.15, -0.15, vUv.y + vUv.x * 0.22);
  vec3  c = uBg + uTint * g * 0.5;
  if (uPhaseMode < 1.5) c += dye;

  // The multiphase view, and the two channels the display pass carries for it come out
  // of the nine taps the dye already paid for.
  //
  // Keyed on temperature, because that is what the cooling actually decides. The two
  // stable phases at the reference pressure sit at 49 K and 7800 K, and between them is
  // the branch where dP/dn < 0 -- the gas that cannot stay where it is. So the
  // transition is placed on that branch and nowhere else: 142.7 K to 4211.5 K, the two
  // turning points of P_eq(n), which in log2 are the two numbers below. Cold gas comes
  // out one side of it and warm gas the other, and the only thing in between is the gas
  // that is genuinely in between.
  //
  // Two tones, and both of them are the chapter's own accent -- the cold phase in it,
  // the warm phase in a pale wash of it. No hue appears here that does not appear
  // everywhere else on the page, which is the difference between one more palette and a
  // colourbar laid over the picture.
  //
  // Brightness follows the density and not the temperature, on a third-power root so
  // that a factor of 160 between the phases becomes a factor of five on screen. The cold
  // phase is where the mass is and should read as substance; the warm phase fills most of
  // the box and would drown it if the two were lit alike.
  if (uPhaseMode > 0.5 && uPhaseMode < 1.5) {
    float ph = smoothstep(7.1564, 12.0401, d0.z);
    float amp = clamp(0.42 * pow(max(d0.w, 1e-3), 0.30), 0.0, 1.0);
    c += mix(uCold, uWarm, ph) * amp * uPhase;
  } else if (uPhaseMode > 1.5) {
    // Temperature alone, diverging about the calibrated 762 K equilibrium. The
    // two sides are normalised separately so the 49 K and 7800 K stable phases
    // reach equal strength even though they are not equally far away in log T.
    float cold = smoothstep(0.0, 1.0, clamp((9.57404 - d0.z) / 3.95933, 0.0, 1.0));
    float hot  = smoothstep(0.0, 1.0, clamp((d0.z - 9.57404) / 3.35522, 0.0, 1.0));
    c += (uCold * cold + uWarm * hot) * uPhase;
  }

  // A soft shoulder rather than a clamp. Where the old expression saturated, sh was
  // exactly 1 over a whole region and the edge of that region was a level set of the
  // interpolant -- an angular, grid-following contour that read as a facet on the
  // brightest fronts, which are the ones you look at. 1 - exp(-x) is linear at the
  // onset, bounded by one, and has no contour anywhere.
  // The fronts, if they are being drawn -- and on this page they are not. At beta = 1
  // they are long-lived, and with the field itself drawn they were two overlays
  // competing for the same picture. uShockLift = 0 skips the smoothed texture and its
  // four taps here, and the pass that fills it does not run either; see smoothDisp.
  // The exponent is a display knob and nothing here touches the solver: a front
  // weakens as it propagates, so its convergence falls, so a steeper exponent fades an
  // old front faster while leaving a fresh one where it was.
  if (uShockLift > 0.0) {
    float x = max(texBS(uDispS, vUv, gridSize).y, 0.0) * uShockGain;
    float sh = 1.0 - exp(-x);
    c += mix(uAccent, vec3(1.0), 0.30) * pow(sh, uShockPow) * uShockLift;
  }

  // The field. Signed about zero so it grains rather than lights, in the chapter's
  // own accent, and scaled by its own strength so only a field that has been done
  // something to shows up. uFieldGain = 0 removes it entirely, at the cost of one
  // texture fetch.
  if (uFieldGain > 0.0) {
    vec2 L = texture(uLic, vUv).xy;
    // |B| still modulates it, but from a floor rather than from zero: with the fronts
    // gone this is the structural element of the picture, and a region the field has
    // not been done anything to should still read as having a field in it.
    c += uAccent * L.x * mix(uFieldFloor, 1.0, clamp(L.y, 0.0, 1.8)) * uFieldGain;
  }

  float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
  c *= 1.0 / (1.0 + 0.62 * l);

  float d = length((vUv - 0.5) * vec2(uRes.x / uRes.y, 1.0));
  c *= 1.0 - uVignette * smoothstep(0.30, 1.05, d);

  c += (h1(uvec2(gl_FragCoord.xy), uint(uTime * 60.0)) - 0.5) * uGrain;
  outColor = vec4(max(c, 0.0), 1.0);
}`;

  // Diagnostics, two maxima and two sums through one reduction so the whole
  // measurement is a single readback:
  //   x = ctot, max      y = |u|^2, sum      z = rho, sum      w = |u|, max
  const F_METRICS = HEAD + EOS + `
uniform sampler2D uU, uB;
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);
  float r  = max(U.x, smallr);
  float ri = 1.0 / r;
  vec2  u  = U.yz * ri;
  float b2 = dot(B.xy, B.xy);
  float pg = adia > 0.5 ? pfromE(B.w, r, u.x, u.y, b2) : r * cs2;
  float c2 = sound2(r, pg);
  float d2 = 0.5 * (b2 * ri + c2);
  float cx = sqrt(d2 + sqrt(max(d2 * d2 - c2 * B.x * B.x * ri, 0.0)));
  float cy = sqrt(d2 + sqrt(max(d2 * d2 - c2 * B.y * B.y * ri, 0.0)));
  outColor = vec4(abs(u.x) + cx + abs(u.y) + cy, dot(u, u), U.x, length(u));
}`;

  // The temperature, cell by cell, for the readout's min/median/max. A median is not a
  // reduction -- it needs the whole distribution -- so this one goes back to the CPU
  // whole, and it goes back on the divergence diagnostic's own sixth-of-a-hertz cadence
  // for the same reason: a readback drains a queue a couple of hundred steps a second
  // deep. Twelve thousand texels every six seconds is a stall nobody can see; the same
  // readback every frame would be the frame budget.
  //   x = T in kelvin      y = rho in code units
  const F_TEMP = HEAD + EOS + `
uniform sampler2D uU, uB;
uniform float tScale;
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec4 U = texelFetch(uU, c, 0);
  vec4 B = texelFetch(uB, c, 0);
  float r = max(U.x, smallr), ri = 1.0 / r;
  float pg = adia > 0.5 ? pfromE(B.w, r, U.y * ri, U.z * ri, dot(B.xy, B.xy)) : r * cs2;
  outColor = vec4(max(pg * ri * tScale, 0.0), r, 0.0, 1.0);
}`;

  // On demand only, so it costs nothing per frame: how well the cleaning is working,
  // and how strong the field and its current sheets are.
  //   x = |div B| dx, max      y = (div B dx)^2, sum
  //   z = |B|, sum             w = |J| dx, max
  //
  // Absolute, and normalised by the *mean* field strength on the way out, because
  // normalising per cell by the local |B| does not survive a weak mean field: at
  // beta = 2 the turbulence reverses the field locally, |B| passes through zero at the
  // nulls, and the ratio there went to 4 and the current to 37 while the absolute
  // divergence had barely moved. That was the diagnostic breaking, not the cleaning.
  const F_DIVB = HEAD + `
uniform sampler2D uB;
uniform ivec2 size;
vec2 bf(ivec2 c){
  ivec2 p = ivec2((c.x + size.x) % size.x, (c.y + size.y) % size.y);
  return texelFetch(uB, p, 0).xy;
}
void main(){
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec2 b0 = bf(c);
  vec2 be = bf(c + ivec2(1, 0)), bw = bf(c - ivec2(1, 0));
  vec2 bn = bf(c + ivec2(0, 1)), bs = bf(c - ivec2(0, 1));
  float div = 0.5 * ((be.x - bw.x) + (bn.y - bs.y));
  float jz  = 0.5 * ((be.y - bw.y) - (bn.x - bs.x));
  outColor = vec4(abs(div), div * div, length(b0), abs(jz));
}`;

  const F_REDUCE = HEAD + `
uniform sampler2D uSrc;
uniform ivec2 srcSize;
void main(){
  ivec2 o = ivec2(gl_FragCoord.xy) * 4;
  float mx = 0.0, mw = 0.0;
  vec2  s  = vec2(0.0);
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
      ivec2 p = o + ivec2(i, j);
      if (p.x >= srcSize.x || p.y >= srcSize.y) continue;
      vec4 v = texelFetch(uSrc, p, 0);
      mx = max(mx, v.x);
      mw = max(mw, v.w);
      s += v.yz;
    }
  }
  outColor = vec4(mx, s, mw);
}`;

  // ------------------------------------------------------------- gl helpers

  function compile(gl, type, src) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
      throw new Error('shader: ' + gl.getShaderInfoLog(s));
    }
    return s;
  }

  function program(gl, vsrc, fsrc) {
    const p = gl.createProgram();
    gl.attachShader(p, compile(gl, gl.VERTEX_SHADER, vsrc));
    gl.attachShader(p, compile(gl, gl.FRAGMENT_SHADER, fsrc));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
      throw new Error('link: ' + gl.getProgramInfoLog(p));
    }
    const u = {};
    const n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
    for (let i = 0; i < n; i++) {
      const nm = gl.getActiveUniform(p, i).name.replace(/\[0\]$/, '');
      u[nm] = gl.getUniformLocation(p, nm);
    }
    return { p, u };
  }

  function makeFBO(gl, w, h, internal, format, type, filter, wrap) {
    const tex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, wrap);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, wrap);
    gl.texImage2D(gl.TEXTURE_2D, 0, internal, w, h, 0, format, type, null);
    const fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) {
      throw new Error('incomplete framebuffer ' + internal);
    }
    return {
      tex, fbo, w, h,
      bind(unit) { gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, tex); return unit; }
    };
  }

  // Static value noise for the LIC to average along. A texture rather than a hash
  // in the shader: eleven samples a pixel is eleven bilinear taps this way and
  // forty-four integer hashes the other.
  function makeNoise(gl, n) {
    const data = new Uint8Array(n * n);
    for (let i = 0; i < data.length; i++) data[i] = Math.random() * 256;
    const tex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, n, n, 0, gl.RED, gl.UNSIGNED_BYTE, data);
    return {
      tex,
      bind(unit) { gl.activeTexture(gl.TEXTURE0 + unit); gl.bindTexture(gl.TEXTURE_2D, tex); return unit; }
    };
  }

  // Two attachments on one framebuffer, so the state and the field are written by
  // the same pass. The draw-buffer list is per-framebuffer state in ES 3.0, so it
  // is set once here and every single-target framebuffer keeps its default.
  function makeMRT(gl, a, b) {
    const fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, a.tex, 0);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, b.tex, 0);
    gl.drawBuffers([gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1]);
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) {
      throw new Error('incomplete MRT framebuffer');
    }
    return { fbo, w: a.w, h: a.h, u: a, f: b };
  }

  // ------------------------------------------------------------- presets

  const BASE = {
    accent: [0.62, 0.72, 0.98],
    tint: [0.020, 0.024, 0.040],
    tau: 0.10,
    mach: 0.55,
    dustGain: 1.00,
    dyeGain: 0.75,
    dyeDiss: 0.75,
    pointSize: 1.5,
    // Saturation point of the dust brightness ramp, as a multiple of the
    // population's measured mean drift. A relative number now, not the absolute
    // tau*M^2 coefficient the unmagnetised build used -- see F_PDRIFT.
    driftRef: 2.0,
    brown: 0.0,
    bhat: [1.0, 0.0],
    stream: 0.0,
    stirGain: 1.00,
    zeta: 0.14,
    grain: 0.030,
    vignette: 0.62
  };

  const PRESETS = {
    hero:      { accent: [0.50, 0.64, 1.00], tau: 0.14, dustGain: 1.48, stirGain: 1.25 },
    bio:       { accent: [0.96, 0.80, 0.50], tau: 0.10, dustGain: 1.05 },
    picdust:   { accent: [1.00, 0.58, 0.18], tau: 0.030, dustGain: 1.45, pointSize: 1.2, stirGain: 1.15 },
    dfmm:      { accent: [0.26, 0.95, 0.72], tau: 0.85,  dustGain: 2.24, pointSize: 1.75, brown: 0.08 },
    mhd:       { accent: [0.60, 0.46, 1.00], tau: 0.12,  dustGain: 1.81, bhat: [0.94, 0.34], zeta: 0.08 },
    cosmicray: { accent: [0.22, 0.86, 1.00], tau: 0.020, dustGain: 2.40, stream: 1.5, bhat: [0.94, 0.34], pointSize: 1.15 },
    phrike:    { accent: [1.00, 0.24, 0.48], tau: 0.06,  dustGain: 1.38, stirGain: 1.35, zeta: 0.26 },
    papers:    { accent: [0.94, 0.82, 0.52], tau: 0.09,  dustGain: 0.95, stirGain: 0.85 },
    blog:      { accent: [0.56, 0.72, 1.00], tau: 0.11,  dustGain: 0.95, stirGain: 0.8 },
    iron:      { accent: [1.00, 0.13, 0.10], tau: 0.55,  dustGain: 1.00, pointSize: 1.6, stirGain: 0.65, vignette: 0.74, zeta: 0.28 },
    cactus:    { accent: [0.54, 0.84, 0.32], tau: 0.95,  dustGain: 2.06, pointSize: 1.7, brown: 0.16, stirGain: 0.42, grain: 0.040 },
    skate:     { accent: [1.00, 0.48, 0.06], tau: 0.20,  dustGain: 1.19, pointSize: 1.3, stirGain: 1.1 },
    contact:   { accent: [0.92, 0.86, 0.58], tau: 0.10,  dustGain: 1.05, stirGain: 0.9 },
    // The play chapter needs its own entry: resolvePreset falls back to hero for an
    // unknown key rather than to BASE, silently, so without this the panel's chapter
    // would quietly retune the servo to hero's driving while you were reading numbers
    // off it.
    play:      { accent: [0.62, 0.78, 1.00], tau: 0.12,  dustGain: 1.30, stirGain: 1.05 }
  };

  function resolvePreset(key) {
    return Object.assign({}, BASE, PRESETS[key] || PRESETS.hero);
  }

  // Palettes for the play panel. 'chapter' means the accent each chapter already
  // carries, which is the default and the only one the rest of the site knows about;
  // the others override it everywhere at once, because every coloured thing on the page
  // -- dye, fronts, field, grains -- is drawn in the accent by construction.
  const PALETTES = {
    chapter:  null,
    ice:      [0.62, 0.80, 1.00],
    amber:    [1.00, 0.74, 0.34],
    ember:    [1.00, 0.42, 0.26],
    viridian: [0.32, 0.94, 0.70],
    violet:   [0.66, 0.50, 1.00],
    mono:     [0.86, 0.88, 0.94]
  };

  // ---------------------------------------------------------------- engine

  function Field(canvas) {
    const gl = canvas.getContext('webgl2', {
      alpha: false, depth: false, stencil: false, antialias: false,
      premultipliedAlpha: false, preserveDrawingBuffer: false,
      powerPreference: 'default'
    });
    if (!gl) throw new Error('no webgl2');
    if (!gl.getExtension('EXT_color_buffer_float')) throw new Error('no float render targets');
    if (gl.getParameter(gl.MAX_DRAW_BUFFERS) < 2) throw new Error('no MRT');

    const reduced = global.matchMedia && global.matchMedia('(prefers-reduced-motion: reduce)').matches;
    const mobile = Math.min(global.innerWidth, global.innerHeight) < 720;

    const PART_SIDE = 320;
    const CBLOCK = 8;
    const SMALLR = 1e-3;
    const SMALLC = 1e-3;
    const CS = 1.0;

    // ------------------------------------------------------------------- UNITS
    //
    // The whole calibration, and the reason the cooling can be switched on without the
    // picture moving. Every number here comes out of the equilibrium curve of the
    // cooling function in ISM above -- not out of a fit to it, out of it.
    //
    // Lambda changes sign twice along an isobar, so the equilibrium pressure P_eq(n)
    // turns over twice, at nH = 3.7305 and nH = 34.011, and the branch between them has
    // dlnP/dlnn = -0.98: gas there cannot stay where it is. rho_code = 1 is put at the
    // geometric centre of that branch,
    //
    //   n_ref  = 11.263964 cm^-3      T_ref = 762.209 K      mu_ref = 1.215638
    //   T2_ref = 627.0030 K           (= T/mu, which is what the solver integrates)
    //
    // and the velocity unit is that state's own isothermal sound speed,
    // scale_v = 2.283e5 cm/s. Then
    //
    //   scale_T2 = (mH/kB) scale_v^2 = 627.003042 K
    //
    // which is T2_ref exactly -- not approximately, and not by tuning. It is the same
    // statement twice: if the velocity unit is the sound speed at the reference state
    // then (P/rho) = 1 there in code units, and (P/rho) = 1 is what T2 = scale_T2 means.
    // So P_code = rho * 1 = 1 at rho = 1, which is precisely the isothermal page's own
    // pressure at its own reference density. Switching the cooling on therefore changes
    // nothing about the pressure field at all. It changes only how the gas answers being
    // compressed -- which is the entire point, and the only thing that should change.
    //
    // Everything the shaders need follows: T2 = (P/rho)_code * scale_T2 in kelvin,
    // T = mu T2, and nH = rho_code * n_ref in cm^-3.
    //
    // The one thing left free is the length unit, and it is a real knob rather than a
    // leftover: it sets scale_t = scale_l/scale_v and so how many cooling times fit into
    // a box crossing. See the boxPc control.
    const UNITS = {
      nRef: 11.263964,          // cm^-3 at rho_code = 1
      tRef: 762.2087,           // K
      t2Ref: 627.003042,        // K, and equal to scaleT2 by construction
      muRef: 1.2156380,
      scaleT2: 627.003042,      // K per unit of (P/rho) in code units
      scaleV: 2.28324185e5,     // cm/s -- the isothermal sound speed is 1 code unit
      scaleD: 2.46108581e-23,   // g/cm^3 at rho_code = 1
      tauCool: 7.15769e12,      // s, the cooling time at the reference state
      // the two turning points of P_eq(n), which are what the phase map keys on
      nCold: 34.0108, nWarm: 3.7305,
      tCold: 142.656, tWarm: 4211.494,
      PC: 3.0856776e18,         // cm
      MYR: 3.15576e13           // s
    };
    // seconds per code time unit, from the one free parameter
    function scaleT() { return Math.max(cfg.boxPc, 0.05) * UNITS.PC / UNITS.scaleV; }
    // kelvin per unit of (P/rho) in code units, mu folded in -- the display's T, not T2
    function tScale() { return UNITS.muRef * UNITS.scaleT2; }
    // cooling time / box-crossing time, and therefore the cooling length in cells
    function coolPerCross() { return UNITS.tauCool / scaleT(); }

    // Reconstruction: 'ppm' or 'plm'. Backend only, no control. PLM is the
    // reference's own configuration and is what runs while the box is being driven
    // either way -- a parabola fitted across a freshly injected discontinuity is
    // monotonised variable by variable, so a state that is monotone in every
    // variable separately can still be unphysical taken together, and with a
    // magnetic field among those variables there is more to get wrong.



    // Plasma beta = P / (B^2/2). At beta = 5 with rho = cs = 1 the seed field is
    // sqrt(2/5) = 0.632, so the Alfven speed is 0.63 against a sound speed of 1 and the
    // rms flow at Mach 0.69 is mildly *super*-Alfvenic. That is a different regime from
    // where this page started: the turbulence wins against the field rather than the
    // other way round, so the mean field gets tangled and amplified instead of
    // organising the flow, and a scroll can roll up again -- the Kelvin-Helmholtz
    // threshold is 2 vA = 1.26 against the 1.4 a gesture deposits, so the layer is now
    // unstable rather than marginal.
    //
    // Measured across beta: the Courant sum barely moves (11 at beta = 1, 11.5 at
    // beta = 0.5, 10-12 at beta = 2), because ctot is set by the Alfven speed in the
    // most rarefied cells and a stronger field resists being rarefied -- the two effects
    // very nearly cancel. A weaker field is the dirtier one for the cleaning, not the
    // stronger: div B rms went 0.0039 -> 0.0052 -> 0.0066 as beta went 0.5 -> 1 -> 2.
    // (plasma beta now lives in cfg.beta; see b0Now)
    
    // Dedner's parabolic term, as the factor exp(-psiDamp * ch * dt/dx) applied
    // once per step. Larger damps the divergence error faster and closer to where
    // it was made; too large and psi cannot carry it out of the box at all.


    // Large-scale drag, as an inverse timescale: see the note in the solver. 1/0.625
    // sim time units, a little under a box-scale turnover at the Mach number this
    // page holds.


    // Grain charge-to-mass ratio: charge_parameter in the reference, one value for
    // the whole population. At 100 the gyrofrequency is q|B| = 141, the gyration
    // is resolved thirty times over per step, and the gyroradius is a twelfth of a
    // cell -- the grains are stuck to field lines.


    // The field visualisation. FIELD_VIS = false removes the pass and the fetch;
    // setFieldVis(false) does the same at runtime, and the page binds it to a key.

    // A standard deviation in linear colour before the accent multiplies it, against a
    // background of about 0.02 -- the LIC is normalised to unit variance, so this is a
    // real level and not an arbitrary coefficient.
    //
    // 0.037 is the level the hand-set version worked out to: it had a normalisation
    // constant of 3.3 against a true standard deviation of 0.289/sqrt(12.7) = 0.081,
    // so 3.3 * 0.081 * 0.14 = 0.037. Stronger than this read as noise rather than as a
    // field, so it is back where it was. The floor is back at zero with it, so |B|
    // modulates the grain from nothing exactly as before.

    const FIELD_FLOOR = 0.0;
    // ------------------------------------------------------------------ the knobs
    //
    // Everything below this line used to be a const, and the play panel is the reason it
    // is not. One object, one manifest, one setter: js/play.js renders whatever
    // controls() returns, so a knob added here appears on the page with no change to the
    // panel, and the panel can never disagree with the engine because set() returns the
    // value it actually took.
    //
    // Two of them cannot be live and the manifest says so. `grid` reallocates every
    // texture; `beta` and `mhd` re-seed the field, which is cheap (one pass into the
    // field half of the state pair) but throws away psi and whatever structure the field
    // had. Everything else is a uniform read per frame or a JS branch, so it takes effect
    // on the next frame.
    const cfg = {
      mhd: true, beta: 5.0,
      solver: 3,               // 0 LLF, 1 HLL, 2 HLLC, 3 HLLD
      recon: 'ppm', slopeType: 2, plmWindow: 1.6,
      // The equation of state. Isothermal is the page and stays the default; the energy
      // equation exists for the cooling option, which cannot work without one. gamma is
      // 5/3 rather than the reference namelist's 1.4 because the cooling implemented
      // here is the atomic branch with no H2 -- a monatomic gas, whose heat capacity 1.4
      // would misstate by a fifth.
      adia: false, gamma: 5 / 3,
      // The cooling, off by default, and the one switch that turns it on turns the
      // energy equation on with it -- a gas that cannot change its temperature cannot
      // have two phases, so the two are not independent choices and are not offered as
      // two. boxPc is the length unit in parsecs, which is the only free parameter left
      // in UNITS and the one that decides how fast the gas is allowed to cool.
      cool: false, boxPc: 2,
      psiDamp: 0.4, powell: 1.0, fric: 1.6, dtdxMax: 0.060, cfl: 0.8,
      tier: mobile ? 1 : 0,
      charge: 100, charged: true,
      fieldVis: true, fieldGain: 0.037,
      shockVis: false, shockLift: 0.36,
      palette: 'chapter',
      phase: 'off', phaseGain: 0.40
    };
    const cfg0 = Object.assign({}, cfg);
    let phaseBeforeCooling = cfg.phase;

    // Overrides the panel lays on top of whatever chapter preset is current. Only keys
    // the user has actually touched appear here, so an untouched knob still follows the
    // chapter as it always did.
    const ovr = {};
    let ovrCount = 0;

    function b0Now() { return cfg.mhd ? Math.sqrt(2 * CS * CS / Math.max(cfg.beta, 1e-3)) : 0.0; }

    const SOLVERS = ['LLF', 'HLL', 'HLLC', 'HLLD'];
    function solverName() {
      const r = (cfg.recon === 'ppm' && agitated <= 0) ? 'PPM (CW84) · MonCen slopes'
              : cfg.slopeType > 0 ? 'MonCen slopes (slope_type 2)' : 'piecewise constant';
      return SOLVERS[cfg.solver] + ' · ' + r + (cfg.mhd ? ' · Dedner GLM' : ' · unmagnetised') + ' · unsplit';
    }

    // Grain and step in pixels of the LIC target, and they have to match each other:
    // step further than the noise is correlated and consecutive samples along the walk
    // are independent again. Coarser and longer than the first pass at this, which read
    // as noise; at five pixels with a twelve-step walk the streaks are long enough to
    // be silk instead, and the sheared regions of the field become legible.
    const NOISE_PX = 5.0;
    const STEP_PX = 5.0;
    const WAVE_LEN = 0.055;      // travelling kernel wavelength, in box heights
    const WAVE_SPEED = 0.35;     // and its rate, wavelengths per wall second
    const LIC_SCALE = 0.5;       // resolution, relative to the canvas

    // The interactions are scaled up on this page, and here is why. At beta = 1 the
    // field carries as much energy as the gas and resists being pushed around, so an
    // impulse that reads as a perturbation on the live site barely registers here --
    // measured, a scroll flick moved the rms Mach by two per cent.
    //
    // Whether it rolls up is a race with the field. The mean field is vertical and the
    // layer a scroll deposits is horizontal, so the field threads it at right angles and
    // tension holds it against Kelvin-Helmholtz until delta_u > 2 vA. At beta = 1 that
    // was 2.83 sound speeds against the 1.4 one gesture deposits -- comfortably stable,
    // Alfven waves and current sheets instead of vortices. At beta = 2 it is 2.0 against
    // the same 1.4, which is marginal: the layer is on the edge of rolling, and beta is
    // the knob either way.
    const DRIVE_GAIN = 2.0;

    // Lower than the hydrodynamic build's 0.075, because ctot is about twice as
    // large at the same Mach number: this keeps the fixed branch of
    // min(dt_default, dt_Courant) the binding one in the developed state, which is
    // what makes the pace of the picture independent of the flow's peaks.


    // Front rendering. The 5x5 binomial spreads a one-cell spike over about two and
    // a half cells, which takes roughly a factor of two off its peak, and the soft
    // shoulder takes about another two off the mid-range where most of the visible
    // front is -- so the gain carries the first and the lift carries the second.
    // Grid-level smoothing of the convergence, in cells. The B-spline reconstruction
    // is what removes the facets; this is only here to take the front from one cell
    // wide to something a display interpolant has any information about at all, so it
    // is deliberately small -- at 1.0 the fronts went soft and stopped reading as
    // fronts.
    // The fronts are not drawn on this page. shockVis = true puts them back, and the
    // page binds it to a key; everything they need is still here and still measured.

    const SHOCK_BLUR = 0.40;
    const SHOCK_GAIN = 2.2;

    const SHOCK_POW = 2.6;   // 2.0 is the live page; higher fades a propagating front faster

    // The servo's gains. The plant here is slower and less dissipative than the
    // unmagnetised one: at beta = 1 the field resists the stirring, so it takes a
    // large amplitude to get the flow moving, and then the Alfvenic cascade carries
    // that energy to the grid more slowly than a hydrodynamic one does, so the
    // amplitude that got the box going is far above the one that holds it.
    // Proportional control alone loses that race -- measured, the rms was still
    // climbing through 1.4 while the drive was being cut at its maximum rate -- so
    // there is a derivative term, which stops the servo pushing while an earlier
    // push is still arriving.
    const SERVO_D = 2.0;
    const SERVO_CLAMP = 0.35;
    const SERVO_GAIN = 0.30;
    // and a bounded output, which is all the anti-windup a multiplicative servo
    // needs: the drive is what makes the vortices and the dye plumes, so it is not
    // allowed to reach zero however far off the rms is.
    const AMP_MIN = 0.30, AMP_MAX = 60;

    // [gridH, gasStepsPerFrame, grainFrac, dprCap]. n/gridH is held constant, so
    // every rung moves the picture at the same apparent speed:
    //   rate = M * fps * n * cfg.dtdxMax / gridH  ~  0.058 box heights per wall
    // second at rms Mach 0.69, which is what the live site measures.
    //
    // The top rung is coarser than the live site's 176 because an MHD step costs
    // more: seven primitives through the limiter instead of four, two textures on a
    // thirteen-cell stencil instead of one, and HLLD's five waves and four square
    // roots instead of HLLC's three waves and none. At 96 cells and 2.24 steps a
    // frame the cost per frame comes out level with the live site's, which is the
    // budget this page has to fit in.
    const TIERS = [
      [ 96, 2.24, 0.90, 1.30],
      [ 80, 1.86, 0.75, 1.20],
      [ 64, 1.49, 0.55, 1.05],
      [ 56, 1.30, 0.40, 1.00]
    ];
    let tier = mobile ? 1 : 0;   // mirrors cfg.tier
    let ceiling = 0;

    const P = {
      god:    program(gl, VERT, isothermalShader(F_GODUNOV)),
      god3:   program(gl, VERT, isothermalShader(F_GODUNOV3)),
      godA:   program(gl, VERT, F_GODUNOV),
      god3A:  program(gl, VERT, F_GODUNOV3),
      cmax0:  program(gl, VERT, F_CMAX_STATE),
      cmaxN:  program(gl, VERT, F_CMAX_DOWN),
      stir:   program(gl, VERT, F_STIR),
      splat:  program(gl, VERT, F_SPLAT),
      shear:  program(gl, VERT, F_SHEAR),
      blast:  program(gl, VERT, F_BLAST),
      init:   program(gl, VERT, F_INIT_ISO),
      initA:  program(gl, VERT, F_INIT),
      bseed:  program(gl, VERT, F_BSEED_ISO),
      bseedA: program(gl, VERT, F_BSEED),
      energ:  program(gl, VERT, F_ENERGIZE),
      cool:   program(gl, VERT, F_COOL),
      flow:   program(gl, VERT, F_FLOW),
      disp:   program(gl, VERT, F_DISP),
      part:   program(gl, VERT, F_PART),
      pdrift: program(gl, VERT, F_PDRIFT),
      pdraw:  program(gl, V_PARTICLE, F_PARTICLE),
      lic:    program(gl, VERT, F_LIC),
      dsm:    program(gl, VERT, F_DSMOOTH),
      comp:   program(gl, VERT, F_COMPOSITE),
      metric: program(gl, VERT, F_METRICS),
      temp:   program(gl, VERT, F_TEMP),
      divb:   program(gl, VERT, F_DIVB),
      reduce: program(gl, VERT, F_REDUCE)
    };

    const quad = gl.createVertexArray();
    gl.bindVertexArray(quad);
    const vb = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vb);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    const emptyVAO = gl.createVertexArray();

    let S = null, flow = null, disp = null, disp2 = null, part = null, lic = null, noise = null;
    let cmax = [], met = null, red = [], pdr = [];
    let grid = { w: 0, h: 0 }, pSide = 0, nPart = 0, nDraw = 0, dpr = 1, aspect = 1;
    let owned = [];

    const state = {
      preset: resolvePreset('hero'),
      target: resolvePreset('hero'),
      live: null,
      blend: 1,
      wall: 0, time: 0, steps: 0,
      running: true, fps: 60,
      amp: 4.0,
      ctot: 2 * CS, maxSig: 1, machRms: 0, machMax: 0, dens: 1, div: 0,
      divb: 0, divbRms: 0, bmean: 1, jmax: 0,
      // the temperature distribution, refreshed on the divergence diagnostic's cadence
      tmin: 0, tmed: 0, tmax: 0
    };
    const pending = { splats: [], shear: 0, blasts: [] };

    // The chapter's preset, with the panel's overrides on top of it. Only touched keys
    // are in ovr, so an untouched knob still follows the chapter from one to the next.
    function activePreset() {
      const p = state.live || state.preset;
      if (!ovrCount && !PALETTES[cfg.palette]) return p;
      const out = Object.assign({}, p, ovr);
      const pal = PALETTES[cfg.palette];
      if (pal) out.accent = pal;
      return out;
    }
    function machNow() { return Math.max(state.machRms, 0.5); }

    function drawQuad(f) {
      if (f) { gl.bindFramebuffer(gl.FRAMEBUFFER, f.fbo); gl.viewport(0, 0, f.w, f.h); }
      else { gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.viewport(0, 0, canvas.width, canvas.height); }
      gl.bindVertexArray(quad);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function mk(w, h, i, f, t, fil, wr, keep) {
      const o = makeFBO(gl, w, h, i, f, t, fil, wr);
      if (!keep) owned.push(o);
      return o;
    }
    // A double-buffered pair of pairs: (state, field) written together through MRT,
    // swapped together, so the two can never fall out of step.
    function mkDoubleState(w, h, keep) {
      const N = gl.NEAREST, R = gl.REPEAT;
      const a0 = mk(w, h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R, keep);
      const a1 = mk(w, h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R, keep);
      const b0 = mk(w, h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R, keep);
      const b1 = mk(w, h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, R, keep);
      let A = makeMRT(gl, a0, b0), B = makeMRT(gl, a1, b1);
      return {
        w, h,
        get read() { return A; }, get write() { return B; },
        swap() { const s = A; A = B; B = s; }
      };
    }

    function eos(pr) {
      gl.uniform1f(pr.u.cs, CS);
      gl.uniform1f(pr.u.cs2, CS * CS);
      gl.uniform1f(pr.u.smallr, SMALLR);
      // Programs that never mention these get null locations, and uniform1f(null, x) is
      // a defined no-op, so one setter serves every program that includes EOS.
      gl.uniform1f(pr.u.adia, cfg.adia ? 1 : 0);
      gl.uniform1f(pr.u.gam, cfg.gamma);
    }
    let warming = false;
    function dtdxCap() { return warming ? 1e9 : cfg.dtdxMax; }
    function dtUniforms(pr) {
      gl.uniform1f(pr.u.cfl, cfg.cfl);
      gl.uniform1f(pr.u.smallc, SMALLC);
      gl.uniform1f(pr.u.dtdxMax, dtdxCap());
    }

    let allocW = -1, allocH = -1, allocTier = -1;

    function allocate(force) {
      const t = TIERS[tier];
      dpr = Math.min(global.devicePixelRatio || 1, t[3]);
      const cw = Math.max(1, Math.round(canvas.clientWidth * dpr));
      const ch = Math.max(1, Math.round(canvas.clientHeight * dpr));
      if (!force && cw === allocW && ch === allocH && tier === allocTier) return;
      allocW = cw; allocH = ch; allocTier = tier;

      const prev = owned;
      owned = [];
      canvas.width = cw; canvas.height = ch;
      aspect = cw / ch;

      grid.h = t[0];
      grid.w = Math.max(32, Math.round(t[0] * aspect / 4) * 4);

      const R = gl.REPEAT, C = gl.CLAMP_TO_EDGE, L = gl.LINEAR, N = gl.NEAREST;
      S = mkDoubleState(grid.w, grid.h);
      // velocity and field together, filterable: the grains gather both at once
      flow = mk(grid.w, grid.h, gl.RGBA16F, gl.RGBA, gl.HALF_FLOAT, L, R);
      // Four channels rather than two: the dye and the convergence as before, plus log T
      // and the density for the multiphase view. They ride in the same texture because
      // the composite's Catmull-Rom already fetches all four components of it, so the
      // phase map costs nine taps that were being thrown away rather than nine more.
      disp = mk(grid.w, grid.h, gl.RGBA16F, gl.RGBA, gl.HALF_FLOAT, L, R);
      disp2 = mk(grid.w, grid.h, gl.RG16F, gl.RG, gl.HALF_FLOAT, L, R);
      lic = mk(Math.max(2, Math.round(cw * LIC_SCALE)), Math.max(2, Math.round(ch * LIC_SCALE)),
               gl.RG16F, gl.RG, gl.HALF_FLOAT, L, C);
      if (!noise) noise = makeNoise(gl, 256);
      met = mk(grid.w, grid.h, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C);

      cmax = [];
      let a = grid.w, b = grid.h;
      do {
        a = Math.max(1, Math.ceil(a / CBLOCK));
        b = Math.max(1, Math.ceil(b / CBLOCK));
        cmax.push(mk(a, b, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C));
      } while (a > 1 || b > 1);

      red = [];
      let rw = grid.w, rh = grid.h;
      while (rw > 4 || rh > 4) {
        rw = Math.max(1, Math.ceil(rw / 4));
        rh = Math.max(1, Math.ceil(rh / 4));
        red.push(mk(rw, rh, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C));
      }

      // Grains live for the life of the context: a tier change draws a shorter
      // prefix of the same population, which is a uniform random subsample.
      if (!part) {
        pSide = PART_SIDE;
        nPart = pSide * pSide;
        part = mkDoubleState(pSide, pSide, true);
        let pw = pSide, ph = pSide;
        do {
          pw = Math.max(1, Math.ceil(pw / CBLOCK));
          ph = Math.max(1, Math.ceil(ph / CBLOCK));
          pdr.push(mk(pw, ph, gl.RGBA32F, gl.RGBA, gl.FLOAT, N, C, true));
        } while (pw > 1 || ph > 1);
        seedParticles();
      }
      nDraw = Math.max(4096, Math.round(nPart * t[2]));

      for (let i = 0; i < prev.length; i++) { gl.deleteFramebuffer(prev[i].fbo); gl.deleteTexture(prev[i].tex); }

      reset();
      warmUp();
      writeFlow();
      driftReduce();
      paint(activePreset());
    }

    // Spin the initial condition up into developed turbulence before anything is
    // shown. The step count scales with the grid, because a box crossing costs
    // ctot/(C M) * gridH steps: 448 at gridH = 112 is the same half-crossing the
    // live site buys with 700 at 176.
    function warmUp() {
      warming = true;
      const n = Math.round(4 * grid.h);
      for (let i = 0; i < n; i++) {
        stepGas();
        if (i % 48 === 47) { measure(); servo(0.8); }
      }
      warming = false;
      cflReduce();
    }

    function reset() {
      const I = cfg.adia ? P.initA : P.init;
      gl.useProgram(I.p);
      eos(I);
      gl.uniform1f(I.u.b0, b0Now());
      drawQuad(S.read);
      drawQuad(S.write);
      state.time = 0; state.steps = 0;
      cflReduce();
      const probe = new Float32Array(4);
      gl.bindFramebuffer(gl.FRAMEBUFFER, S.read.fbo);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, probe);
      if (!(probe[0] > 0.1)) throw new Error('state reads back as ' + probe[0]);
    }

    function seedParticles() {
      const data = new Float32Array(nPart * 4);
      for (let i = 0; i < nPart; i++) {
        data[i * 4] = Math.random();
        data[i * 4 + 1] = Math.random();
      }
      const zero = new Float32Array(nPart * 4);
      for (const f of [part.read.u, part.write.u]) {
        gl.bindTexture(gl.TEXTURE_2D, f.tex);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, pSide, pSide, 0, gl.RGBA, gl.FLOAT, data);
      }
      for (const f of [part.read.f, part.write.f]) {
        gl.bindTexture(gl.TEXTURE_2D, f.tex);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, pSide, pSide, 0, gl.RGBA, gl.FLOAT, zero);
      }
    }

    // ------------------------------------------------------------- stepping

    function cflReduce() {
      gl.useProgram(P.cmax0.p);
      eos(P.cmax0);
      gl.uniform1i(P.cmax0.u.uSrc, S.read.u.bind(0));
      gl.uniform1i(P.cmax0.u.uSrcB, S.read.f.bind(1));
      gl.uniform2i(P.cmax0.u.srcSize, grid.w, grid.h);
      drawQuad(cmax[0]);
      for (let i = 1; i < cmax.length; i++) {
        gl.useProgram(P.cmaxN.p);
        gl.uniform1i(P.cmaxN.u.uSrc, cmax[i - 1].bind(0));
        gl.uniform2i(P.cmaxN.u.srcSize, cmax[i - 1].w, cmax[i - 1].h);
        drawQuad(cmax[i]);
      }
    }
    function cmaxTex() { return cmax[cmax.length - 1]; }
    function dtdxNow() {
      return Math.min(dtdxCap(), Math.min(cfg.cfl / SMALLC, cfg.cfl / Math.max(state.ctot, 1e-20)));
    }
    function dtNow() { return dtdxNow() / grid.h; }

    function applyForcing(pr) {
      const dx = 1 / grid.h;
      const t = state.time;
      const c0 = [0.5 + 0.31 * Math.cos(t * 1.7 * 0.11 + 0), 0.5 + 0.27 * Math.sin(t * 1.31 * 0.11)];
      const c1 = [0.5 + 0.31 * Math.cos(t * 1.7 * 0.162 + 1), 0.5 + 0.27 * Math.sin(t * 1.31 * 0.162 + 2)];
      const a = pr.accent;
      gl.useProgram(P.stir.p);
      gl.uniform1f(P.stir.u.adia, cfg.adia ? 1 : 0);
      dtUniforms(P.stir);
      gl.uniform1i(P.stir.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.stir.u.uB, S.read.f.bind(1));
      gl.uniform1i(P.stir.u.uCmax, cmaxTex().bind(2));
      gl.uniform2f(P.stir.u.c0, c0[0], c0[1]);
      gl.uniform2f(P.stir.u.c1, c1[0], c1[1]);
      gl.uniform1f(P.stir.u.amp, state.amp);
      gl.uniform1f(P.stir.u.radius, 0.055);
      gl.uniform1f(P.stir.u.dyeRadius, 0.020);
      gl.uniform1f(P.stir.u.zeta, pr.zeta);
      gl.uniform1f(P.stir.u.aspect, aspect);
      gl.uniform1f(P.stir.u.dxCell, dx);
      gl.uniform1f(P.stir.u.smallr, SMALLR);
      gl.uniform3f(P.stir.u.dye, 2.1, 0.0, 0.0);
      drawQuad(S.write); S.swap();

      for (const s of pending.splats) {
        gl.useProgram(P.splat.p);
        gl.uniform1f(P.splat.u.adia, cfg.adia ? 1 : 0);
        gl.uniform1i(P.splat.u.uU, S.read.u.bind(0));
        gl.uniform1i(P.splat.u.uB, S.read.f.bind(1));
        gl.uniform2f(P.splat.u.point, s.x, s.y);
        gl.uniform2f(P.splat.u.delta, s.dx, s.dy);
        gl.uniform1f(P.splat.u.aspect, aspect);
        gl.uniform1f(P.splat.u.radius, s.radius);
        gl.uniform1f(P.splat.u.ink, s.ink);
        gl.uniform1f(P.splat.u.smallr, SMALLR);
        drawQuad(S.write); S.swap();
      }
      pending.splats.length = 0;

      if (Math.abs(pending.shear) > 1e-4) {
        const amp = Math.max(-6, Math.min(6, pending.shear));   // entry-point bounded; backstop
        const full = 1.5 * machNow();
        gl.useProgram(P.shear.p);
        gl.uniform1f(P.shear.u.adia, cfg.adia ? 1 : 0);
        gl.uniform1i(P.shear.u.uU, S.read.u.bind(0));
        gl.uniform1i(P.shear.u.uB, S.read.f.bind(1));
        gl.uniform1f(P.shear.u.amp, amp);
        const w = Math.max(3 / grid.h, 0.012);
        gl.uniform1f(P.shear.u.width, w);
        gl.uniform1f(P.shear.u.cycles, Math.round(0.44 * aspect / (6.2831853 * w)));
        gl.uniform1f(P.shear.u.seed, amp * 0.50);
        gl.uniform1f(P.shear.u.phase, state.time * 0.21);
        gl.uniform1f(P.shear.u.smallr, SMALLR);
        const s = Math.min(0.9, Math.abs(amp) / full) * 0.60;
        gl.uniform3f(P.shear.u.band, a[0] * s, a[1] * s, a[2] * s);
        drawQuad(S.write); S.swap();
        pending.shear *= 0.55;
        if (Math.abs(pending.shear) < 1e-3) pending.shear = 0;
      }

      for (const b of pending.blasts) {
        gl.useProgram(P.blast.p);
        gl.uniform1f(P.blast.u.adia, cfg.adia ? 1 : 0);
        gl.uniform1i(P.blast.u.uU, S.read.u.bind(0));
        gl.uniform1i(P.blast.u.uB, S.read.f.bind(1));
        gl.uniform2f(P.blast.u.point, b.x, b.y);
        gl.uniform1f(P.blast.u.aspect, aspect);
        gl.uniform1f(P.blast.u.radius, b.radius);
        gl.uniform1f(P.blast.u.dens, b.dens);
        gl.uniform1f(P.blast.u.kick, b.kick);
        gl.uniform1f(P.blast.u.smallr, SMALLR);
        drawQuad(S.write); S.swap();
      }
      pending.blasts.length = 0;
    }

    function stepGas() {
      const pr = activePreset();
      cflReduce();
      applyForcing(pr);
      // Not during the warm-up: nothing is on screen for it to resolve, and the
      // parabola nearly doubles the cost of a step, which is the whole boot stall.
      const ppm = cfg.recon === 'ppm' && agitated <= 0 && !warming;
      const G = cfg.adia ? (ppm ? P.god3A : P.godA) : (ppm ? P.god3 : P.god);
      gl.useProgram(G.p);
      gl.uniform1f(G.u.solver, cfg.solver);
      eos(G);
      dtUniforms(G);
      gl.uniform1i(G.u.uU, S.read.u.bind(0));
      gl.uniform1i(G.u.uB, S.read.f.bind(1));
      gl.uniform1i(G.u.uCmax, cmaxTex().bind(2));
      gl.uniform2i(G.u.size, grid.w, grid.h);
      gl.uniform1f(G.u.slopeType, cfg.slopeType);
      gl.uniform1f(G.u.dyeDiss, pr.dyeDiss);
      gl.uniform1f(G.u.dxCell, 1 / grid.h);
      gl.uniform1f(G.u.psiDamp, cfg.psiDamp);
      gl.uniform1f(G.u.powell, cfg.powell);
      gl.uniform1f(G.u.fric, cfg.fric);
      drawQuad(S.write); S.swap();
      // and then the cooling, operator-split, at the full step, once per step -- the
      // order amr_step uses. Nothing between the two passes touches the Courant texture,
      // so the dt the cooling integrates over is the dt the fluxes were built with.
      if (cfg.adia) coolStep();
      state.time += dtNow();
      state.steps++;
    }

    function coolStep() {
      gl.useProgram(P.cool.p);
      eos(P.cool);
      dtUniforms(P.cool);
      gl.uniform1i(P.cool.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.cool.u.uB, S.read.f.bind(1));
      gl.uniform1i(P.cool.u.uCmax, cmaxTex().bind(2));
      gl.uniform1f(P.cool.u.dxCell, 1 / grid.h);
      gl.uniform1f(P.cool.u.scaleT2, UNITS.scaleT2);
      gl.uniform1f(P.cool.u.scaleT, scaleT());
      gl.uniform1f(P.cool.u.nRef, UNITS.nRef);
      drawQuad(S.write); S.swap();
    }

    function writeFlow() {
      gl.useProgram(P.flow.p);
      gl.uniform1f(P.flow.u.smallr, SMALLR);
      gl.uniform1i(P.flow.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.flow.u.uB, S.read.f.bind(1));
      drawQuad(flow);

      gl.useProgram(P.disp.p);
      eos(P.disp);
      gl.uniform1i(P.disp.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.disp.u.uB, S.read.f.bind(1));
      gl.uniform2i(P.disp.u.size, grid.w, grid.h);
      gl.uniform1f(P.disp.u.tScale, tScale());
      drawQuad(disp);
      smoothDisp();
    }

    function smoothDisp() {
      if (!cfg.shockVis) return;
      gl.useProgram(P.dsm.p);
      gl.uniform1i(P.dsm.u.uSrc, disp.bind(0));
      gl.uniform2i(P.dsm.u.size, grid.w, grid.h);
      gl.uniform1f(P.dsm.u.sigma, SHOCK_BLUR);
      drawQuad(disp2);
    }

    function stepDust(pr, dt) {
      const a = dt / Math.max(pr.tau, 1e-4);
      gl.useProgram(P.part.p);
      gl.uniform1i(P.part.u.uPart, part.read.u.bind(0));
      gl.uniform1i(P.part.u.uPartZ, part.read.f.bind(1));
      gl.uniform1i(P.part.u.uFlow, flow.bind(2));
      gl.uniform1f(P.part.u.dragA, a);
      gl.uniform1f(P.part.u.dt, dt);
      gl.uniform1f(P.part.u.aspect, aspect);
      gl.uniform2f(P.part.u.gridSize, grid.w, grid.h);
      gl.uniform1f(P.part.u.brown, pr.brown * (dt * 60));
      gl.uniform1f(P.part.u.charge, cfg.charged ? cfg.charge : 0.0);
      gl.uniform1f(P.part.u.stream, pr.stream);
      gl.uniform1ui(P.part.u.uFrame, state.steps >>> 0);
      drawQuad(part.write); part.swap();
    }

    // The population's mean drift, reduced on the GPU and never read back: the
    // vertex shader fetches the tip directly, so the lighting scale is measured
    // every time it is used and nothing has to be predicted from tau and M.
    function driftReduce() {
      let src = part.read.f, first = 1;
      for (const d of pdr) {
        gl.useProgram(P.pdrift.p);
        gl.uniform1i(P.pdrift.u.uSrc, src.bind(0));
        gl.uniform2i(P.pdrift.u.srcSize, src.w, src.h);
        gl.uniform1f(P.pdrift.u.first, first);
        drawQuad(d);
        src = d; first = 0;
      }
    }
    function driftTex() { return pdr[pdr.length - 1]; }

    // -------------------------------------------------------------- painting

    // One pass, before the composite, at half resolution.
    function drawLic() {
      gl.useProgram(P.lic.p);
      gl.uniform1i(P.lic.u.uFlow, flow.bind(0));
      gl.uniform1i(P.lic.u.uNoise, noise.bind(1));
      gl.uniform1f(P.lic.u.aspect, aspect);
      // wall time, not sim time: the packets travel at a rate you can see rather
      // than at the clock the gas happens to be running at
      gl.uniform1f(P.lic.u.time, state.wall);
      gl.uniform1f(P.lic.u.b0, Math.max(b0Now(), 0.05));
      // one noise texel = NOISE_PX target pixels, and the walk steps by the same
      gl.uniform1f(P.lic.u.stepUv, STEP_PX / lic.h);
      gl.uniform1f(P.lic.u.noiseFreq, lic.h / (256.0 * NOISE_PX));
      gl.uniform1f(P.lic.u.waveLen, WAVE_LEN);
      gl.uniform1f(P.lic.u.waveSpeed, WAVE_SPEED);
      drawQuad(lic);
    }

    function composite(pr) {
      gl.useProgram(P.comp.p);
      gl.uniform1i(P.comp.u.uLic, lic.bind(3));
      gl.uniform1f(P.comp.u.uFieldGain, cfg.fieldVis ? cfg.fieldGain : 0.0);
      gl.uniform1f(P.comp.u.uFieldFloor, FIELD_FLOOR);
      gl.uniform1i(P.comp.u.uDisp, disp.bind(0));
      gl.uniform1i(P.comp.u.uDispS, disp2.bind(4));
      gl.uniform2f(P.comp.u.gridSize, grid.w, grid.h);
      gl.uniform2f(P.comp.u.uRes, canvas.width, canvas.height);
      gl.uniform3f(P.comp.u.uBg, 0.014, 0.015, 0.021);
      gl.uniform3f(P.comp.u.uTint, pr.tint[0], pr.tint[1], pr.tint[2]);
      gl.uniform3f(P.comp.u.uAccent, pr.accent[0], pr.accent[1], pr.accent[2]);
      gl.uniform1f(P.comp.u.uTime, state.wall);
      gl.uniform1f(P.comp.u.uGrain, pr.grain);
      gl.uniform1f(P.comp.u.uDyeGain, pr.dyeGain);
      gl.uniform1f(P.comp.u.uVignette, pr.vignette);
      gl.uniform1f(P.comp.u.uShockGain, SHOCK_GAIN / Math.max(0.22 * machNow(), 0.08));
      gl.uniform1f(P.comp.u.uShockLift, cfg.shockVis ? cfg.shockLift : 0.0);
      gl.uniform1f(P.comp.u.uShockPow, SHOCK_POW);
      // The original phase overlay follows the chapter accent. The temperature-only
      // view is fixed blue-to-red because its colour carries a physical sign.
      const a = pr.accent;
      const thermal = cfg.phase === 'thermal';
      gl.uniform1f(P.comp.u.uPhase, cfg.phase === 'off' ? 0.0 : cfg.phaseGain);
      gl.uniform1f(P.comp.u.uPhaseMode, thermal ? 2.0 : (cfg.phase === 'temperature' ? 1.0 : 0.0));
      gl.uniform3f(P.comp.u.uCold, thermal ? 0.08 : a[0], thermal ? 0.30 : a[1], thermal ? 1.35 : a[2]);
      gl.uniform3f(P.comp.u.uWarm, thermal ? 1.20 : 1 - 0.58 * (1 - a[0]),
        thermal ? 0.12 : 1 - 0.58 * (1 - a[1]), thermal ? 0.04 : 1 - 0.58 * (1 - a[2]));
      drawQuad(null);
    }

    function drawDust(pr) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, canvas.width, canvas.height);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE);
      gl.useProgram(P.pdraw.p);
      gl.uniform1i(P.pdraw.u.uPart, part.read.u.bind(0));
      gl.uniform1i(P.pdraw.u.uPartZ, part.read.f.bind(1));
      gl.uniform1i(P.pdraw.u.uDrift, driftTex().bind(2));
      gl.uniform1i(P.pdraw.u.uPW, pSide);
      gl.uniform1f(P.pdraw.u.uPointSize, pr.pointSize * dpr * 3.9);
      gl.uniform1f(P.pdraw.u.uDriftRef, pr.driftRef);
      gl.uniform1f(P.pdraw.u.uInvN, 1 / nPart);
      gl.uniform1f(P.pdraw.u.uAlpha, pr.dustGain * 0.30);
      gl.uniform1f(P.pdraw.u.uVignette, pr.vignette);
      gl.uniform1f(P.pdraw.u.uAspect, aspect);
      gl.uniform3f(P.pdraw.u.uColor, pr.accent[0], pr.accent[1], pr.accent[2]);
      gl.bindVertexArray(emptyVAO);
      gl.drawArrays(gl.POINTS, 0, nDraw);
      gl.disable(gl.BLEND);
    }

    function paint(pr) { if (cfg.fieldVis) drawLic(); composite(pr); drawDust(pr); }


    // ------------------------------------------------------- the control surface
    //
    // The panel renders whatever this returns, so this is the single place a knob is
    // declared. `live` is false for the two that cannot be: `grid` reallocates every
    // texture, and `beta`/`mhd` re-seed the field.
    function seedField() {
      const B = cfg.adia ? P.bseedA : P.bseed;
      gl.useProgram(B.p);
      eos(B);
      gl.uniform1f(B.u.b0, b0Now());
      // each field target paired with its own state target, so the energy each one
      // carries is consistent with the density and momentum beside it
      gl.uniform1i(B.u.uU, S.read.u.bind(0));
      drawQuad(S.read.f);
      gl.uniform1i(B.u.uU, S.write.u.bind(0));
      drawQuad(S.write.f);
      cflReduce();
      writeFlow();
    }

    // Fill in the energy channel of a state that does not have one yet, or reinterpret
    // one that was written under a different gamma. See F_ENERGIZE: both attachments and
    // a swap, because the pass has to read the channel it is rewriting.
    function energize(fromE, gamPrev) {
      gl.useProgram(P.energ.p);
      eos(P.energ);
      gl.uniform1i(P.energ.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.energ.u.uB, S.read.f.bind(1));
      gl.uniform1f(P.energ.u.fromE, fromE ? 1 : 0);
      gl.uniform1f(P.energ.u.gamPrev, gamPrev);
      drawQuad(S.write); S.swap();
      cflReduce();
      writeFlow();
    }

    // chrome.js owns the page furniture and holds the only copy of that state -- the
    // corner control and the switch in the panel have to be two views of one thing, not
    // two booleans that can disagree, so this reaches for its handle rather than keeping
    // one. Looked up on every call because chrome.js may well be parsed after this file.
    function pageChrome() {
      const c = global.__pageChrome;
      return (c && typeof c.setBare === 'function') ? c : null;
    }

    const OVR_RANGE = {
      tau:       [0.005, 1.5, 0.005],
      mach:      [0.05, 2.5, 0.05],
      zeta:      [0, 1, 0.02],
      dustGain:  [0, 4, 0.05],
      dyeGain:   [0, 3, 0.05],
      pointSize: [0.5, 4, 0.05],
      grain:     [0, 0.12, 0.005],
      vignette:  [0, 1, 0.02],
      stream:    [0, 4, 0.05],
      brown:     [0, 0.5, 0.01]
    };

    function controls() {
      const pr = activePreset();
      const sel = (v) => ({ value: String(v), label: String(v) });
      const list = [
        { key: 'mhd', label: 'Magnetised', group: 'gas', kind: 'toggle', value: cfg.mhd,
          note: 'off sets B = 0 everywhere, which reduces the same solver to isothermal hydrodynamics' },
        { key: 'beta', label: 'Plasma beta', group: 'gas', kind: 'range', value: cfg.beta,
          min: 0.25, max: 40, step: 0.25,
          note: '|B| = sqrt(2/beta), so vA = ' + b0Now().toFixed(2) + ' against cs = 1' },
        { key: 'mach', label: 'Target rms Mach', group: 'gas', kind: 'range',
          value: ovr.mach !== undefined ? ovr.mach : pr.mach, min: 0.05, max: 2.5, step: 0.05,
          note: 'what the driving servos to; the pace of the picture follows it' },
        { key: 'zeta', label: 'Compressive drive', group: 'gas', kind: 'range',
          value: ovr.zeta !== undefined ? ovr.zeta : pr.zeta, min: 0, max: 1, step: 0.02,
          note: '0 is purely solenoidal, which makes eddies; 1 is purely radial, which makes sound' },
        { key: 'fric', label: 'Large-scale drag', group: 'gas', kind: 'range', value: cfg.fric,
          min: 0, max: 8, step: 0.1,
          note: 'two dimensions cascade energy to large scales, where nothing dissipates it; at 0 the servo has no plant and the rms wanders' },
        { key: 'grid', label: 'Grid (cells tall)', group: 'gas', kind: 'select', value: String(tier),
          options: TIERS.map((t, i) => ({ value: String(i), label: t[0] + ' cells' })), live: false },

        { key: 'cool', label: 'ISM cooling', group: 'thermal', kind: 'toggle', value: cfg.cool,
          note: 'atomic heating and cooling -- grain photoelectric, cosmic rays, C+ and O fine '
              + 'structure, Lyman-alpha, grain recombination, CIE metals -- split after each '
              + 'hydrodynamic step, which needs the energy equation and turns it on too. '
              + 'rho = 1 is nH = ' + UNITS.nRef.toFixed(2) + ' cm-3 and T = ' + UNITS.tRef.toFixed(0)
              + ' K, where the equilibrium pressure and the isothermal one are the same number by '
              + 'construction, so switching this on does not move the pressure' },
        { key: 'boxPc', label: 'Box size (parsecs)', group: 'thermal', kind: 'range', value: cfg.boxPc,
          min: 1, max: 32, step: 1,
          note: 'the only free unit left, and it sets the thermal clock: '
              + 'the cooling time is ' + coolPerCross().toFixed(3) + ' box crossings, so its length '
              + 'is ' + (coolPerCross() * grid.h).toFixed(1) + ' cells. Below about 4 pc the cold '
              + 'phase is resolved; by 32 it fragments at the grid scale and the cold clouds are '
              + 'single cells, which is a picture of the mesh and not of the gas' },
        { key: 'gamma', label: 'Adiabatic index', group: 'thermal', kind: 'range', value: cfg.gamma,
          min: 1.05, max: 2, step: 0.05,
          note: '5/3 for a monatomic gas, which is what this cooling assumes: it is the atomic '
              + 'branch with no H2, and 1.4 would misstate the heat capacity by a fifth. Moving it '
              + 'on a running box keeps the pressure and changes only the heat capacity' },

        { key: 'solver', label: 'Riemann solver', group: 'solver', kind: 'select', value: String(cfg.solver),
          options: [{ value: '0', label: 'LLF — one wave' },
                    { value: '1', label: 'HLL — two waves' },
                    { value: '2', label: 'HLLC — three, contact' },
                    { value: '3', label: 'HLLD — five, Alfven' }],
          note: 'HLL has no contact, so it averages the transverse velocity and the dye across every interface, which is exactly what a shear layer is made of' },
        { key: 'recon', label: 'Reconstruction', group: 'solver', kind: 'select', value: cfg.recon,
          options: [{ value: 'plm', label: 'PLM — linear' }, { value: 'ppm', label: 'PPM — parabolic' }],
          note: 'PLM always runs for a moment after you disturb the box either way' },
        { key: 'slopeType', label: 'Slope limiter', group: 'solver', kind: 'select', value: String(cfg.slopeType),
          options: [{ value: '0', label: '0 — piecewise constant' }, { value: '2', label: '2 — MonCen' }] },
        { key: 'dtdxMax', label: 'Fixed timestep dt/dx', group: 'solver', kind: 'range', value: cfg.dtdxMax,
          min: 0.005, max: 0.16, step: 0.005,
          note: 'dt = min(this, Courant); above about 0.06 the flow starts setting the pace instead' },
        { key: 'cfl', label: 'Courant factor', group: 'solver', kind: 'range', value: cfg.cfl,
          min: 0.1, max: 1, step: 0.05 },

        { key: 'psiDamp', label: 'psi damping', group: 'field', kind: 'range', value: cfg.psiDamp,
          min: 0, max: 2, step: 0.02,
          note: "Dedner's parabolic term: larger damps the divergence error nearer where it was made" },
        { key: 'powell', label: 'Powell source', group: 'field', kind: 'range', value: cfg.powell,
          min: 0, max: 1, step: 0.05,
          note: 'removes the field-aligned force a residual divergence would otherwise exert on the gas' },
        { key: 'fieldVis', label: 'Draw the field (LIC)', group: 'field', kind: 'toggle', value: cfg.fieldVis,
          note: 'line integral convolution along B, in the chapter accent. Also on the B key' },
        { key: 'fieldGain', label: 'Field strength', group: 'field', kind: 'range', value: cfg.fieldGain,
          min: 0, max: 0.16, step: 0.002,
          note: 'a standard deviation in linear colour, against a background of about 0.02' },

        { key: 'charged', label: 'Charged grains', group: 'dust', kind: 'toggle', value: cfg.charged,
          note: 'the Lorentz force acts on the drift v - u, since the lab electric field is -u x B' },
        { key: 'charge', label: 'Charge-to-mass', group: 'dust', kind: 'range', value: cfg.charge,
          min: 0, max: 400, step: 5,
          note: 'at 100 the gyroradius is a twelfth of a cell and the grains are stuck to field lines' },
        { key: 'tau', label: 'Stopping time', group: 'dust', kind: 'range',
          value: ovr.tau !== undefined ? ovr.tau : pr.tau, min: 0.005, max: 1.5, step: 0.005,
          note: 'small glues a grain to the gas, large lets it slip and light up' },
        { key: 'dustGain', label: 'Grain brightness', group: 'dust', kind: 'range',
          value: ovr.dustGain !== undefined ? ovr.dustGain : pr.dustGain, min: 0, max: 4, step: 0.05 },
        { key: 'pointSize', label: 'Grain size', group: 'dust', kind: 'range',
          value: ovr.pointSize !== undefined ? ovr.pointSize : pr.pointSize, min: 0.5, max: 4, step: 0.05 },
        { key: 'brown', label: 'Brownian kick', group: 'dust', kind: 'range',
          value: ovr.brown !== undefined ? ovr.brown : pr.brown, min: 0, max: 0.5, step: 0.01 },
        { key: 'stream', label: 'Field-aligned streaming', group: 'dust', kind: 'range',
          value: ovr.stream !== undefined ? ovr.stream : pr.stream, min: 0, max: 4, step: 0.05 },

        { key: 'shockVis', label: 'Draw the fronts', group: 'display', kind: 'toggle', value: cfg.shockVis,
          note: 'the convergence of the velocity field, which is where a shock is. Also on the S key' },
        { key: 'shockLift', label: 'Front strength', group: 'display', kind: 'range', value: cfg.shockLift,
          min: 0, max: 1.2, step: 0.02 },
        { key: 'dyeGain', label: 'Dye brightness', group: 'display', kind: 'range',
          value: ovr.dyeGain !== undefined ? ovr.dyeGain : pr.dyeGain, min: 0, max: 3, step: 0.05 },
        { key: 'palette', label: 'Palette', group: 'display', kind: 'select', value: cfg.palette,
          options: Object.keys(PALETTES).map((k) => ({ value: k, label: k })),
          note: 'every coloured thing on the page is drawn in the accent, so this moves all of it at once' },
        { key: 'phase', label: 'Colour by', group: 'display', kind: 'select', value: cfg.phase,
          options: [{ value: 'off', label: 'the dye' },
                    { value: 'temperature', label: 'the dye and the phase' },
                    { value: 'thermal', label: 'temperature only — blue / red' }],
          note: cfg.phase === 'thermal'
            ? 'no dye: neutral at ' + UNITS.tRef.toFixed(0) + ' K, increasingly blue toward the '
              + 'cold phase and red toward the hot phase'
            : 'the original phase overlay, with dye preserved; brightness follows density and both '
              + 'tones follow the selected palette' },
        { key: 'phaseGain', label: 'Phase strength', group: 'display', kind: 'range', value: cfg.phaseGain,
          min: 0, max: 1, step: 0.02 },
        { key: 'grain', label: 'Film grain', group: 'display', kind: 'range',
          value: ovr.grain !== undefined ? ovr.grain : pr.grain, min: 0, max: 0.12, step: 0.005 },
        { key: 'vignette', label: 'Vignette', group: 'display', kind: 'range',
          value: ovr.vignette !== undefined ? ovr.vignette : pr.vignette, min: 0, max: 1, step: 0.02 },
        { key: 'bare', label: 'Take the page away', group: 'display', kind: 'toggle',
          value: pageChrome() ? pageChrome().bare() : false,
          note: 'everything but the gas and the switches in the corner. This is the same state the '
              + 'corner control holds -- there is only one -- so it also takes this panel away; the '
              + 'corner, or Escape, brings it back' }
      ];
      void sel;
      return list.map((c) => (c.live === undefined ? Object.assign(c, { live: true }) : c));
    }

    function set(key, value) {
      // Not the engine's state at all: chrome.js holds it, and this is a second view of
      // the same switch rather than a copy of it. Returns what chrome actually took, like
      // every other knob here.
      if (key === 'bare') {
        const c = pageChrome();
        return c ? c.setBare(!!value) : undefined;
      }
      if (key === 'grid') {
        const t = Math.max(0, Math.min(TIERS.length - 1, parseInt(value, 10) || 0));
        if (t !== tier) { tier = t; cfg.tier = t; ceiling = t; try { allocate(true); } catch (e) { dead = true; } }
        return String(tier);
      }
      if (key in OVR_RANGE) {
        const r = OVR_RANGE[key];
        const v = Math.max(r[0], Math.min(r[1], Number(value)));
        if (ovr[key] === undefined) ovrCount++;
        ovr[key] = v;
        // the chapter blend would otherwise overwrite an override mid-transition
        state.preset = Object.assign({}, state.preset, ovr);
        state.target = Object.assign({}, state.target, ovr);
        return v;
      }
      if (!(key in cfg)) return undefined;
      const was = cfg[key];
      if (typeof cfg[key] === 'boolean') cfg[key] = !!value;
      else if (typeof cfg[key] === 'number') cfg[key] = Number(value);
      else cfg[key] = String(value);
      if (key === 'beta') cfg.beta = Math.max(0.25, Math.min(40, cfg.beta));
      if (key === 'charge') cfg.charge = Math.max(0, Math.min(400, cfg.charge));
      if (key === 'solver') cfg.solver = Math.max(0, Math.min(3, Math.round(cfg.solver)));
      if (key === 'boxPc') cfg.boxPc = Math.max(0.25, Math.min(64, cfg.boxPc));
      if (key === 'gamma') cfg.gamma = Math.max(1.05, Math.min(2, cfg.gamma));
      if (key === 'phaseGain') cfg.phaseGain = Math.max(0, Math.min(1, cfg.phaseGain));
      if (key === 'phase' && !['off', 'temperature', 'thermal'].includes(cfg.phase)) cfg.phase = 'off';
      if ((key === 'beta' || key === 'mhd') && cfg[key] !== was) seedField();
      // The cooling implies the energy equation, and the energy channel of a box that has
      // been running isothermally holds nothing -- so it is filled in from the pressure
      // the gas already has, which is rho cs^2, which is the pressure the calibration was
      // built to make continuous. Nothing about the picture moves at the moment of the
      // switch; what changes is what happens next.
      if (key === 'cool' && cfg.cool !== was) {
        cfg.adia = cfg.cool;
        if (cfg.cool) {
          if (cfg.phase !== 'thermal') phaseBeforeCooling = cfg.phase;
          cfg.phase = 'thermal';
        } else if (cfg.phase === 'thermal') {
          cfg.phase = phaseBeforeCooling;
        }
        energize(0, cfg.gamma);
      }
      // and moving gamma under a live gamma-law box reinterprets the channel rather than
      // rewriting the state: the gas keeps its pressure and changes only its heat capacity
      if (key === 'gamma' && cfg.adia && cfg.gamma !== was) energize(1, was);
      frozen = false;
      return cfg[key];
    }

    // The divergence is an on-demand diagnostic -- it costs a readback, which drains a
    // queue a couple of hundred steps a second deep -- so it is refreshed on every tenth
    // call rather than every one. The panel asks about three times a second, which puts
    // this at about a sixth of a hertz, and nothing else calls it.
    let roCount = 0;
    function readout() {
      const pr = activePreset();
      if (roCount-- <= 0) { roCount = 19; divbNow(); tempNow(); }
      return {
        fps: state.fps, grid: [grid.w, grid.h], steps: state.steps,
        rms: state.machRms, machMax: state.machMax, ctot: state.ctot, dtdx: dtdxNow(),
        dens: state.dens, solver: solverName(),
        recon: (cfg.recon === 'ppm' && agitated <= 0) ? 'ppm' : 'plm',
        beta: cfg.beta, b0: b0Now(),
        alfvenMach: state.machRms / Math.max(b0Now() / Math.sqrt(Math.max(state.dens, 1e-6)), 1e-6),
        charge: cfg.charged ? cfg.charge : 0,
        divbRms: state.divbRms, divbMax: state.divb, bmean: state.bmean,
        tau: pr.tau, amp: state.amp, nDraw: nDraw, sub: TIERS[tier][1],
        // the thermal state, in units a reader can check against the literature rather
        // than against this file
        cool: cfg.cool, gamma: cfg.gamma, boxPc: cfg.boxPc,
        tmin: state.tmin, tmed: state.tmed, tmax: state.tmax,
        nRef: UNITS.nRef, nH: state.dens * UNITS.nRef,
        tauCoolMyr: UNITS.tauCool / UNITS.MYR,
        crossMyr: scaleT() / UNITS.MYR,
        coolPerCross: coolPerCross(),
        coolCells: coolPerCross() * grid.h,
        // the fixed step and the Courant one, side by side, because once the gas is
        // two-phase the warm phase's sound speed decides which of the two is binding
        dtdxFixed: cfg.dtdxMax, dtdxCourant: cfg.cfl / Math.max(state.ctot, 1e-20)
      };
    }

    function resetConfig() {
      const wasBeta = cfg.beta, wasMhd = cfg.mhd, wasTier = tier, wasAdia = cfg.adia;
      Object.keys(cfg0).forEach((k) => { cfg[k] = cfg0[k]; });
      phaseBeforeCooling = cfg.phase;
      Object.keys(ovr).forEach((k) => { delete ovr[k]; });
      ovrCount = 0;
      state.preset = resolvePreset('play');
      state.target = resolvePreset('play');
      if (cfg.beta !== wasBeta || cfg.mhd !== wasMhd) seedField();
      // back to isothermal, which means the energy channel goes back to carrying nothing;
      // the gas keeps its density structure and its pressure becomes rho cs^2 again
      if (cfg.adia !== wasAdia) energize(0, cfg.gamma);
      tempNow();
      if (cfg.tier !== wasTier) set('grid', cfg.tier);
      const c = pageChrome();
      if (c) c.setBare(false);
      return controls();
    }

    // ----------------------------------------------------------- diagnostics

    const buf = new Float32Array(4 * 64);
    let measureCountdown = 2;

    // The reduction chain and the single readback both metrics share: two maxima
    // in .x/.w and two sums in .y/.z, whatever the metric pass put there.
    function reduceAndRead() {
      let src = met;
      for (const d of red) {
        gl.useProgram(P.reduce.p);
        gl.uniform1i(P.reduce.u.uSrc, src.bind(0));
        gl.uniform2i(P.reduce.u.srcSize, src.w, src.h);
        drawQuad(d);
        src = d;
      }
      gl.bindFramebuffer(gl.FRAMEBUFFER, src.fbo);
      const n = Math.min(64, src.w * src.h);
      try { gl.readPixels(0, 0, src.w, src.h, gl.RGBA, gl.FLOAT, buf); }
      catch (e) { return null; }
      let mx = 0, sy = 0, sz = 0, mw = 0;
      for (let i = 0; i < n; i++) {
        mx = Math.max(mx, buf[i * 4]);
        sy += buf[i * 4 + 1];
        sz += buf[i * 4 + 2];
        mw = Math.max(mw, buf[i * 4 + 3]);
      }
      return { mx, sy, sz, mw };
    }

    function measure() {
      gl.useProgram(P.metric.p);
      eos(P.metric);
      gl.uniform1i(P.metric.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.metric.u.uB, S.read.f.bind(1));
      drawQuad(met);
      const r = reduceAndRead();
      if (!r) return;
      const cells = grid.w * grid.h;
      if (isFinite(r.mw)) state.machMax = r.mw / CS;
      if (isFinite(r.sy)) state.machRms = Math.sqrt(Math.max(0, r.sy / cells)) / CS;
      if (isFinite(r.sz) && r.sz > 0) state.dens = r.sz / cells;
      if (isFinite(r.mx) && r.mx > 0) state.ctot = r.mx;
      state.maxSig = state.machMax * CS + 2 * CS;

      if (!isFinite(state.machRms) || state.machMax > 400) { state.amp = 2.0; reset(); }
    }

    // How well the cleaning is working. Off the render path: a readback drains a
    // queue that is a couple of hundred steps a second deep, so this is called by
    // hand, not per frame.
    function divbNow() {
      gl.useProgram(P.divb.p);
      gl.uniform1i(P.divb.u.uB, S.read.f.bind(0));
      gl.uniform2i(P.divb.u.size, grid.w, grid.h);
      drawQuad(met);
      const r = reduceAndRead();
      if (!r) return null;
      const cells = grid.w * grid.h;
      state.bmean = r.sz / cells;
      const bref = Math.max(state.bmean, 1e-6);
      state.divb = r.mx / bref;                       // max |div B| dx / <|B|>
      state.divbRms = Math.sqrt(Math.max(0, r.sy / cells)) / bref;
      state.jmax = r.mw / bref;
      return { maxRel: state.divb, rms: state.divbRms, bmean: state.bmean, jmax: state.jmax };
    }

    // The temperature distribution: min, median, max, in kelvin. Off the render path,
    // on divbNow's cadence, and for the same reason. A median is not a reduction, so this
    // is the one diagnostic that brings a whole grid back rather than a reduced tip -- and
    // with the cooling off it brings nothing back at all, because a gas at one
    // temperature has one temperature and it is known in closed form.
    let tbuf = null, tsort = null;
    function tempNow() {
      if (!cfg.adia) {
        const T = tScale();          // (P/rho) = 1 at the reference state, so T = T_ref
        state.tmin = state.tmed = state.tmax = T;
        return { min: T, med: T, max: T };
      }
      gl.useProgram(P.temp.p);
      eos(P.temp);
      gl.uniform1i(P.temp.u.uU, S.read.u.bind(0));
      gl.uniform1i(P.temp.u.uB, S.read.f.bind(1));
      gl.uniform1f(P.temp.u.tScale, tScale());
      drawQuad(met);
      const n = grid.w * grid.h;
      if (!tbuf || tbuf.length < n * 4) { tbuf = new Float32Array(n * 4); tsort = new Float32Array(n); }
      gl.bindFramebuffer(gl.FRAMEBUFFER, met.fbo);
      try { gl.readPixels(0, 0, grid.w, grid.h, gl.RGBA, gl.FLOAT, tbuf); }
      catch (e) { return null; }
      const t = tsort.subarray(0, n);
      for (let i = 0; i < n; i++) t[i] = tbuf[i * 4];
      t.sort();
      state.tmin = t[0]; state.tmed = t[n >> 1]; state.tmax = t[n - 1];
      return { min: state.tmin, med: state.tmed, max: state.tmax };
    }

    // Hold the target Mach. Proportional-derivative, and only ever called once per
    // measurement -- per step it is a runaway. The derivative is per measurement
    // interval, which is the natural unit here: it is the rate the error is
    // changing at, and subtracting it means the servo eases off while the flow is
    // still answering rather than pushing until the error changes sign.
    let servoPrev = -1;
    function servo(gain) {
      const pr = activePreset();
      const want = pr.mach * pr.stirGain;
      const s = Math.max(want, 0.4);
      const err = (want - state.machRms) / s;
      const d = servoPrev < 0 ? 0 : (state.machRms - servoPrev) / s;
      servoPrev = state.machRms;
      const u = err - SERVO_D * d;
      const g = gain * Math.max(-SERVO_CLAMP, Math.min(SERVO_CLAMP, u));
      state.amp = Math.max(AMP_MIN, Math.min(AMP_MAX, state.amp * (1 + g)));
    }

    // ------------------------------------------------------------------ loop

    function lerpPreset(dt) {
      if (state.blend >= 1) return;
      state.blend = Math.min(1, state.blend + dt * 1.35);
      const k = state.blend, e = k * k * (3 - 2 * k);
      const cur = state.preset, tgt = state.target, out = {};
      for (const key in tgt) {
        const a = cur[key], b = tgt[key];
        if (Array.isArray(b)) out[key] = b.map((v, i) => a[i] + (v - a[i]) * e);
        else if (typeof b === 'number') out[key] = a + (b - a) * e;
        else out[key] = b;
      }
      state.live = out;
    }

    const TARGET_FPS = 60;
    const FRAME_MIN = 1 / TARGET_FPS - 0.002;
    const FPS_FLOOR = 30;
    let last = 0, frames = 0, fpsT = 0, slow = 0, fast = 0, tierChanges = 0;
    let frozen = false, lost = false, dead = false;
    // Two different stillnesses, and they are not the same thing. `frozen` is the
    // reduced-motion freeze and the allocation-failure freeze: it is the page's, and any
    // interaction lifts it. `paused` is the reader's, from the control in the corner, and
    // nothing lifts it but the reader.
    let paused = false;
    let frameRequest = 0;

    function queueFrame() {
      if (!frameRequest && !dead && !paused && state.running) {
        frameRequest = requestAnimationFrame(frame);
      }
    }

    function stopFrames() {
      if (!frameRequest) return;
      cancelAnimationFrame(frameRequest);
      frameRequest = 0;
    }

    function retune(d) {
      if (tierChanges >= 4) return;
      const next = tier + d;
      if (next < 0 || next > TIERS.length - 1) return;
      tier = next;
      if (d > 0) ceiling = tier;
      tierChanges++;
      credit = 0;
      try { allocate(true); } catch (e) { frozen = true; dead = true; }
    }

    const IMP_CAP = 1.0, IMP_LEAK = 0.7;
    // Wall seconds of PLM left to run. Any interaction refreshes it.
    let agitated = 0;
    function agitate(){ agitated = cfg.plmWindow; }
    // The pointer's own impulse budget. A splat is a velocity increment, not a
    // force, and site.js can emit one every 28 ms -- so a ceiling on each one does
    // not bound a stroke, which is the same mistake the scroll shear had. Measured
    // before this: thirty-six splats a second at the old ceiling of 1.0 M is 86
    // sound speeds a second deposited into a spot a twentieth of the box across. It
    // took the peak Mach to 8 here and 22 on the magnetised page, drove ctot past
    // the point where the Courant condition takes over from the default timestep,
    // and stayed there for seconds afterwards because nothing had removed the
    // energy: the picture visibly slowed down and did not speed back up.
    //
    // So the momentum is metered and the ink is not. The ink is a passive scalar
    // that costs the solver nothing and it is what makes a stroke visible, so it is
    // scaled by the gesture; the momentum draws on a budget that refills at a fixed
    // rate. A flick lands at a bit under half of what it used to, and a held scribble is bounded
    // by the leak rather than by how long it is held: a three-second scribble now
    // delivers a fifth of what it did, and a flick about half, with the timestep
    // measured unthrottled through both -- including a scribble held in one spot,
    // which is the case that concentrates it worst.
    const SPLAT_CAP = 0.45;      // per event, in units of the rms Mach
    const SPLAT_BUDGET = 1.60;   // how much may be outstanding
    const SPLAT_LEAK = 7.0;      // and how fast that comes back, per second
    let pimp = 0;
    let impulse = 0;

    let credit = 0, driftCountdown = 0;
    function tick(pr) {
      const n = TIERS[tier][1];
      credit += n;
      let stepped = false;
      while (credit >= 1) { stepGas(); credit -= 1; stepped = true; }
      if (stepped) writeFlow();
      stepDust(pr, Math.max(n * dtNow(), 1e-9));
      // the lighting scale changes on the timescale the whole population changes
      // on, so it does not need refreshing every frame
      if (--driftCountdown <= 0) { driftCountdown = 8; driftReduce(); }
      paint(pr);
    }

    function bench(n) {
      const pr = activePreset();
      const px = new Float32Array(4);
      const sync = function () {
        gl.bindFramebuffer(gl.FRAMEBUFFER, S.read.fbo);
        gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, px);
      };
      for (let i = 0; i < 4; i++) tick(pr);
      sync();
      const t0 = performance.now();
      for (let i = 0; i < n; i++) tick(pr);
      sync();
      return (performance.now() - t0) / n;
    }

    function frame(now) {
      frameRequest = 0;
      if (dead || paused || !state.running) return;
      queueFrame();
      if (lost || gl.isContextLost()) return;
      if (!last) { last = now; return; }
      const el = (now - last) / 1000;
      if (el < FRAME_MIN) return;
      last = now;
      const dtWall = Math.min(el, 1 / 12);
      state.wall += dtWall;

      frames++; fpsT += dtWall;
      if (state.wall < 2.5) { frames = 0; fpsT = 0; }
      else if (fpsT > 0.75) {
        state.fps = frames / fpsT; frames = 0; fpsT = 0;
        if (state.fps < FPS_FLOOR - 6 && tier < TIERS.length - 1) {
          fast = 0; if (++slow >= 3) { slow = 0; retune(+1); }
        } else if (state.fps > FPS_FLOOR + 12 && tier > ceiling) {
          slow = 0; if (++fast >= 5) { fast = 0; retune(-1); }
        } else if (state.fps > FPS_FLOOR) { slow = 0; fast = 0; }
      }

      lerpPreset(dtWall);
      impulse = Math.max(0, impulse - IMP_LEAK * machNow() * DRIVE_GAIN * dtWall);
      agitated = Math.max(0, agitated - dtWall);
      pimp = Math.max(0, pimp - SPLAT_LEAK * machNow() * dtWall);
      const pr = activePreset();

      if (!frozen) {
        tick(pr);
        if (--measureCountdown <= 0) { measureCountdown = 40; measure(); servo(SERVO_GAIN); }
      } else {
        paint(pr);
      }
      if (reduced && state.wall > 2.5) frozen = true;
    }

    canvas.addEventListener('webglcontextlost', function (e) {
      e.preventDefault(); lost = true; owned = []; part = null; S = null;
    }, false);

    let resizeTimer = null;
    global.addEventListener('resize', function () {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(function () { try { allocate(false); } catch (e) { dead = true; } }, 250);
    });
    document.addEventListener('visibilitychange', function () {
      state.running = !document.hidden; last = 0;
      if (state.running) queueFrame(); else stopFrames();
    });

    allocate(true);
    queueFrame();

    return {
      kind: 'webgl2',
      // The interaction surface site.js speaks, unchanged from the live site:
      // every impulse is a multiple of the flow's own rms Mach, drawn against a
      // budget that refills at a fixed rate, so a held scroll cannot keep
      // refilling a ceiling.
      splat(x, y, dx, dy, radius) {
        // Nothing is queued while paused. Without this the impulses would pile up in
        // `pending` and land all at once the moment the reader unfroze, which is the
        // opposite of what the control is for.
        if (paused) return;
        const M = machNow();
        const sp = Math.hypot(dx, dy);
        if (sp < 1e-6) return;
        const cap = SPLAT_CAP * M;
        // the ink follows the gesture, so the stroke stays visible even when the
        // budget has nothing left to give the gas
        const m = Math.min(1, sp / cap);
        const give = Math.min(sp, cap, Math.max(0, SPLAT_BUDGET * M - pimp));
        const k = give / sp;
        pimp += give;
        pending.splats.push({
          x, y, dx: dx * k, dy: dy * k,
          ink: 0.06 + m * 0.30,
          radius: radius || 0.0030
        });
        if (pending.splats.length > 6) pending.splats.splice(0, pending.splats.length - 6);
        frozen = false;
        agitate();
      },
      shear(amount) {
        if (paused) return;
        const M = machNow();
        const want = Math.max(-1.5, Math.min(1.5, amount / 70)) * 0.22 * M * DRIVE_GAIN;
        const room = Math.max(0, IMP_CAP * M * DRIVE_GAIN - impulse) / 2.22;
        const s = Math.sign(want) * Math.min(Math.abs(want), room);
        if (Math.abs(s) < 1e-5) return;
        impulse += Math.abs(s) * 2.22;
        pending.shear += s;
        frozen = false;
        agitate();
      },
      blast(x, y, amp, radius) {
        if (paused) return;
        const M = machNow();
        const kick = Math.min(1.8 * M, ((amp || 300) / 300) * 0.9 * M,
                              Math.max(0, SPLAT_BUDGET * M - pimp));
        pimp += kick;
        pending.blasts.push({
          x, y,
          dens: 1.0,
          kick: kick,
          radius: radius || 0.010
        });
        frozen = false;
        agitate();
      },
      setPreset(key) {
        const next = resolvePreset(key);
        state.preset = Object.assign({}, activePreset());
        state.target = next;
        state.blend = 0;
        frozen = false;
      },
      accentOf(key) { return resolvePreset(key).accent; },
      // The reader's freeze. Held in the engine rather than in the button so that the
      // corner control and the play chapter are two views of one piece of state and
      // cannot disagree about whether the box is running.
      setPaused(on) {
        paused = !!on;
        if (paused) stopFrames();
        else { last = 0; queueFrame(); }
        return paused;
      },
      paused() { return paused; },
      // The field visualisation, on a switch. Off costs one texture fetch in the
      // composite and skips the LIC pass entirely.
      setFieldVis(on) { cfg.fieldVis = !!on; return cfg.fieldVis; },
      fieldVis() { return cfg.fieldVis; },
      // and the fronts, the same way
      setShockVis(on) { cfg.shockVis = !!on; return cfg.shockVis; },
      shockVis() { return cfg.shockVis; },
      controls: controls,
      set: set,
      readout: readout,
      resetConfig: resetConfig,
      __bench: bench,
      __divb: divbNow,
      __temp: tempNow,
      // The temperature distribution as a histogram in log10 T, which is the shape the
      // two phases show up in: one peak near 50 K, one near 8000 K, and a trough on the
      // branch between them where nothing can sit.
      __thist(nbins, lo, hi) {
        const n = nbins || 40, a = lo || 1.0, b = hi || 5.0;
        if (!tempNow()) return null;
        const cells = grid.w * grid.h;
        const h = new Array(n).fill(0);
        for (let i = 0; i < cells; i++) {
          const T = cfg.adia && tbuf ? tbuf[i * 4] : state.tmed;
          const x = Math.log10(Math.max(T, 1e-3));
          const k = Math.floor((x - a) / (b - a) * n);
          if (k >= 0 && k < n) h[k]++;
        }
        return { lo: a, hi: b, bins: h, cells };
      },
      // A column of primitives -- rho, u, v, Bx, By, psi -- for measuring a
      // profile off the clock.
      __col(ix) {
        const pu = new Float32Array(grid.h * 4), pf = new Float32Array(grid.h * 4);
        gl.bindFramebuffer(gl.FRAMEBUFFER, S.read.u.fbo);
        gl.readPixels(ix | 0, 0, 1, grid.h, gl.RGBA, gl.FLOAT, pu);
        gl.bindFramebuffer(gl.FRAMEBUFFER, S.read.f.fbo);
        gl.readPixels(ix | 0, 0, 1, grid.h, gl.RGBA, gl.FLOAT, pf);
        const out = [];
        for (let j = 0; j < grid.h; j++) {
          const r = Math.max(pu[j * 4], 1e-6);
          out.push([r, pu[j * 4 + 1] / r, pu[j * 4 + 2] / r,
            pf[j * 4], pf[j * 4 + 1], pf[j * 4 + 2]]);
        }
        return out;
      },
      // A sample of the grain population: position, velocity (three components),
      // and the drift it was lit by.
      __grains(n) {
        const side = Math.max(1, Math.min(pSide, Math.round(Math.sqrt(n || 64))));
        const a = new Float32Array(side * side * 4), b = new Float32Array(side * side * 4);
        gl.bindFramebuffer(gl.FRAMEBUFFER, part.read.u.fbo);
        gl.readPixels(0, 0, side, side, gl.RGBA, gl.FLOAT, a);
        gl.bindFramebuffer(gl.FRAMEBUFFER, part.read.f.fbo);
        gl.readPixels(0, 0, side, side, gl.RGBA, gl.FLOAT, b);
        const out = [];
        for (let i = 0; i < side * side; i++) {
          out.push([a[i * 4], a[i * 4 + 1], a[i * 4 + 2], a[i * 4 + 3], b[i * 4], b[i * 4 + 1]]);
        }
        return out;
      },
      stats() {
        const pr = activePreset();
        return {
          gridW: grid.w, gridH: grid.h, nPart, nDraw, fps: state.fps, tier,
          tau: pr.tau, dpr,
          rms: state.machRms, max: state.machMax, div: state.div,
          mach: state.machRms, dens: state.dens,
          solver: solverName(),
          mgLevels: null,
          time: state.time, sub: TIERS[tier][1], steps: state.steps, amp: state.amp,
          recon: (cfg.recon === 'ppm' && agitated <= 0) ? 'ppm' : 'plm',
          timeScale: 1, targetFps: TARGET_FPS, cfl: cfg.cfl,
          drag: 'backward Euler after Cayley Lorentz', driveGain: DRIVE_GAIN,
          // MHD
          beta: cfg.beta, b0: b0Now(), charge: cfg.charged ? cfg.charge : 0, ch: state.ctot,
          psiDamp: cfg.psiDamp, powell: cfg.powell, fric: cfg.fric,
          fieldVis: cfg.fieldVis, shockVis: cfg.shockVis, phase: cfg.phase,
          // thermal
          cool: cfg.cool, adia: cfg.adia, gamma: cfg.gamma, boxPc: cfg.boxPc,
          nRef: UNITS.nRef, scaleT2: UNITS.scaleT2,
          tmin: state.tmin, tmed: state.tmed, tmax: state.tmax,
          coolPerCross: coolPerCross(), coolCells: coolPerCross() * grid.h,
          divb: state.divb, divbRms: state.divbRms, bmean: state.bmean, jmax: state.jmax,
          alfvenMach: state.machRms * CS / Math.max(state.bmean / Math.sqrt(Math.max(state.dens, 1e-6)), 1e-6)
        };
      }
    };
  }

  global.Field = {
    create(canvas) {
      try {
        // site.js holds the instance privately, and this page wants a key bound to
        // the field visualisation, so publish it. Test page; one global.
        const api = Field(canvas);
        global.__mhdField = api;
        return api;
      }
      catch (e) {
        if (global.console) console.warn('[field-mhd] ' + e.message);
        return null;
      }
    }
  };
})(window);

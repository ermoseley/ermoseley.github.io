// ---- the corner EMF: one Riemann problem in two directions at once -------
//
// riemann2d_hlld, from the reference's godunov_utils: the two-dimensional Riemann
// problem at a cell *corner*, where four states meet instead of two, reduced to the
// single number a constrained-transport stage wants from it -- the z-component of the
// EMF along the edge that pierces the corner. Miyoshi & Kusano's HLLD is a
// one-dimensional solver; this is that construction run along x and along y
// simultaneously, with the two fans allowed to overlap, and the EMF read off wherever
// the corner turns out to sit inside the resulting two-dimensional fan.
//
// The header above says why the gas itself is GLM and not CT, and none of that
// changes: this is the extra pass and the second stencil, and it is here because the
// EMF is the one quantity GLM cannot give you -- a cleaning wave keeps div B small,
// where a corner EMF keeps the *staggered* divergence identically zero because the
// same number is differenced into both of the two faces it borders, with opposite
// signs, and cancels. That is a stronger statement about the field than anything
// riemann() above can make, and it costs a corner solve per corner per step.
//
// What it returns. The reference's E is
//
//     Ez_code = u By - v Bx = (v x B)_z
//
// which is MINUS the physical Ez = -(v x B)_z of ideal MHD. The convention is kept
// exactly as the reference writes it, sign and all, because the reference's own
// update differences it with that sign; anything that consumes this has to agree with
// the reference, not with a textbook.
//
// Two dimensions, and what that removes. Bz and vz are absent from every state here.
// This is not an approximation and not a truncation: with Bz = vz = 0 the in-plane
// equations are homogeneous in them and they stay zero exactly, the same argument the
// header makes for the gas. So the reference's C, CLL, CLR, CRL, CRR are all zero,
// which drops C*C from the four total pressures and C*C from B2 in the fast speed,
// and that is the whole of it. Nothing else is dropped. All six -- nine, counting the
// leaves rather than the clauses -- branches of the final selection are here, and so
// are both Alfven-speed pairs, because the x pair and the y pair are what distinguish
// this from four independent one-dimensional solves and the reason a rotational
// discontinuity through the corner does not lose its rotation.
//
// Not ported: the switch_llf_dmin / switch_llf_pmin fallback, which replaces the
// whole solve with a Rusanov average of the four corner EMFs when the density or the
// pressure anywhere in the corner falls below a threshold. Both thresholds default to
// zero in the reference, which disables the branch, and there is no uniform for them
// here; the guards below are what stands in for it, and they are cheaper because they
// are local to the divide that would have failed rather than to the whole corner.
//
// The guards, and why each is the smallest one that cannot move a well-conditioned
// answer. The reference divides by all of these unprotected because a Fortran run can
// afford to trap and inspect; a fragment shader cannot, and one NaN texel propagates
// through the next stencil until the page is black.
//
//   the ustar / vstar denominators : the four rc are each rho*(u - SL) or
//     rho*(SR - u) with SL below every u and SR above every u by at least the fast
//     speed, so all four are strictly positive and the sum cannot vanish for a state
//     that came out of a reconstruction. Floored at 1e-20. Over 33k randomised corners
//     spanning eight decades of density the smallest sum seen was 7e-3, so this one is
//     unreachable in practice and exists for the state that has already gone wrong.
//   (SL - ustar), (SB - vstar) : negative by the same argument, floored at -1e-8;
//     (SR - ustar), (ST - vstar) positive, floored at +1e-8. riemann()'s own choice
//     above, for the same quantity, and kept for consistency with it rather than
//     because the number matters: measured against the fast speed instead of against
//     zero, the contact never came closer to the edge of the fan than 0.10 cx in the
//     same 33k corners, so the floor cannot bind unless every wave speed in the corner
//     is itself of order 1e-8 -- a state with no dynamics in it -- or unless a density
//     or a pressure has gone negative, which is where it does bind and where it is the
//     only thing standing between one bad texel and a black frame.
//   sqrt(rho*) : every star density is a product of positive ratios, so the floor is
//     smallr, the same one riemann() puts on its two star densities, reached through
//     one helper so there is one place to change it. This is the one guard here that
//     does bind on input the engine could plausibly produce -- a corner whose densities
//     span four decades can leave a star density a factor of 40 below smallr -- and
//     where it binds this deliberately stops agreeing with the reference. The
//     justification is the reference's own: a density below the floor is not a density
//     this scheme represents, prim() has already clamped the four inputs, and being
//     consistent with the sibling solver about what a representable density is matters
//     more than tracking a Fortran run into a state the rest of the file has floored.
//   (SAR - SAL), (SAT - SAB) : floored at smallc. Not a fudge factor -- it is the
//     exact lower bound. SAR >= 0 >= SAL always, and whichever sign ustar has, one of
//     the two is at least calfven away from zero while the other is on the far side
//     of it, and calfven is itself max'd with smallc; so SAR - SAL >= smallc is an
//     identity and the floor can only bind on round-off. smallc is floored in turn,
//     in case a caller passes zero.

// find_speed_fast, again. cfast() above hardwires cs2 and takes the normal field
// already squared, which is all the one-dimensional solver ever needs; here the same
// state is asked for its fast speed along x and then along y, so the two in-plane
// components come in separately and trade places between the two calls, and the sound
// speed comes from sound2 so the corner works under either EOS. B2 has lost its third
// term with Bz = 0.
//
// The inner max is round-off insurance and cfast()'s above: d2*d2 - c2*A2/rho equals
// (A2/rho - c2)^2/4 plus transverse terms, so it is non-negative exactly and only
// arithmetic can make it slightly less. The max on the sound speed is a different kind
// of thing and is honestly redundant: it fires only where a pressure floor upstream has
// been missed, and even then the outer root is safe without it, because the inner root
// returns at least |d2| and IEEE makes sqrt(x*x) exactly |x|. It stays because a
// negative c2 is not a small error to propagate into four wave speeds -- zero is the
// defensible fallback for a gas whose pressure has gone negative -- but the test cannot
// observe it and does not pretend to.
//
// The two divides by r are deliberately not folded into one reciprocal, which is the
// one place in this port where the cheaper form was rejected. The reason is that these
// eight numbers set SL, SR, SB and ST, and the reference's EMF selection is genuinely
// discontinuous where one of those crosses zero -- measured on one such state, E jumps
// from -6.7625 to -10.0151 as ST passes through zero, because the two formulas either
// side of the boundary average over different sub-fans and are not required to meet.
// So the corner that happens to sit on a boundary is decided by the last bit of this
// expression, and the only defensible choice is to write it exactly as the reference
// writes it, divide for divide, rather than one bit away from it.
float cfast2(float r, float p, float An, float At){
  float c2 = max(sound2(r, p), 0.0);
  float d2 = 0.5 * ((An * An + At * At) / r + c2);
  return sqrt(d2 + sqrt(max(d2 * d2 - c2 * An * An / r, 0.0)));
}

// 1/sqrt of a star density. They appear only as Alfven speeds |B|/sqrt(rho*), so the
// reciprocal is one operation where a root and a divide are two, and the floor lives
// here rather than at twelve call sites.
float irs(float rs){ return 1.0 / sqrt(max(rs, smallr)); }

float riemann2d(Prim LL, Prim LR, Prim RL, Prim RR, float smallc){
  // The quadrant naming is the reference's and is worth stating, because it is not the
  // (row, column) reading one would guess: the FIRST letter is the side in x, the
  // SECOND the side in y, with L = left/bottom and R = right/top. So LL = (-x,-y),
  // LR = (-x,+y), RL = (+x,-y), RR = (+x,+y). Two lines pin it down on their own --
  // rcLLx = rLL*(uLL - SL) measures LL against the leftmost x speed, so LL is left in
  // x, and rcLRy = rLR*(ST - vLR) measures LR against the topmost y speed, so LR is
  // top in y. Everything downstream agrees: AstarT, the top Alfven-averaged Bx, is
  // built from LR and RR, and BstarR, the right one, from RL and RR.
  float rLL = max(LL.h.x, smallr), rLR = max(LR.h.x, smallr);
  float rRL = max(RL.h.x, smallr), rRR = max(RR.h.x, smallr);
  float uLL = LL.h.y, uLR = LR.h.y, uRL = RL.h.y, uRR = RR.h.y;
  float vLL = LL.h.z, vLR = LR.h.z, vRL = RL.h.z, vRR = RR.h.z;
  float ALL = LL.f.x, ALR = LR.f.x, ARL = RL.f.x, ARR = RR.f.x;   // A = Bx
  float BLL = LL.f.y, BLR = LR.f.y, BRL = RL.f.y, BRR = RR.f.y;   // B = By
  float pLL = LL.f.w, pLR = LR.f.w, pRL = RL.f.w, pRR = RR.f.w;   // gas pressure
  float sc  = max(smallc, 1e-20);   // everything below leans on calfven >= sc

  // eight fast speeds, and only the largest in each direction survives. Along y the
  // normal field is B and the transverse one is A, which is the whole of what "rotate
  // the state" means for a solver that only ever looks at one normal.
  float cx = max(max(cfast2(rLL, pLL, ALL, BLL), cfast2(rLR, pLR, ALR, BLR)),
                 max(cfast2(rRL, pRL, ARL, BRL), cfast2(rRR, pRR, ARR, BRR)));
  float cy = max(max(cfast2(rLL, pLL, BLL, ALL), cfast2(rLR, pLR, BLR, ALR)),
                 max(cfast2(rRL, pRL, BRL, ARL), cfast2(rRR, pRR, BRR, ARR)));

  // the four edges of the two-dimensional fan. Deliberately not clamped through zero,
  // as in riemann(): the selection at the bottom needs the true signs to know which
  // quadrant, or which pair, the corner fell in.
  float SL = min(min(uLL, uLR), min(uRL, uRR)) - cx;
  float SR = max(max(uLL, uLR), max(uRL, uRR)) + cx;
  float SB = min(min(vLL, vLR), min(vRL, vRR)) - cy;
  float ST = max(max(vLL, vLR), max(vRL, vRR)) + cy;

  // each quadrant's own EMF, u By - v Bx, which is what the supersonic branches return
  float ELL = uLL * BLL - vLL * ALL;
  float ELR = uLR * BLR - vLR * ALR;
  float ERL = uRL * BRL - vRL * ARL;
  float ERR = uRR * BRR - vRR * ARR;

  // total pressure, gas plus magnetic, with the Bz term gone
  float PtotLL = pLL + 0.5 * (ALL * ALL + BLL * BLL);
  float PtotLR = pLR + 0.5 * (ALR * ALR + BLR * BLR);
  float PtotRL = pRL + 0.5 * (ARL * ARL + BRL * BRL);
  float PtotRR = pRR + 0.5 * (ARR * ARR + BRR * BRR);

  float rcLLx = rLL * (uLL - SL), rcRLx = rRL * (SR - uRL);
  float rcLRx = rLR * (uLR - SL), rcRRx = rRR * (SR - uRR);
  float rcLLy = rLL * (vLL - SB), rcLRy = rLR * (ST - vLR);
  float rcRLy = rRL * (vRL - SB), rcRRy = rRR * (ST - vRR);

  // ustar and vstar: mass and normal-momentum conservation across the whole corner,
  // one direction each, and the same two-wave argument riemann() uses for its contact
  // -- except that here four states contribute to each and the two transverse
  // pressure differences enter in the pairing the reference chose (LL-RL and LR-RR
  // across x, LL-LR and RL-RR across y). That pairing is what makes the whole function
  // antisymmetric under transposing the two directions -- exchange u with v, Bx with
  // By and the LR quadrant with RL, and E changes sign and nothing else, which the test
  // measures at 2e-12 relative over twenty thousand corners. It is the strongest
  // internal check there is on the quadrant naming: get one slot wrong and it fails.
  float iux = 1.0 / max(rcLLx + rcLRx + rcRLx + rcRRx, 1e-20);
  float ivy = 1.0 / max(rcLLy + rcLRy + rcRLy + rcRRy, 1e-20);
  float ustar = (rcLLx * uLL + rcLRx * uLR + rcRLx * uRL + rcRRx * uRR
                 + (PtotLL - PtotRL + PtotLR - PtotRR)) * iux;
  float vstar = (rcLLy * vLL + rcLRy * vLR + rcRLy * vRL + rcRRy * vRR
                 + (PtotLL - PtotLR + PtotRL - PtotRR)) * ivy;

  // the compression factors (S - q)/(S - qstar), four in x and four in y, each used
  // three times below: once for a star density, once for a star field component and
  // once for the doubly-starred density of the quadrant. Four reciprocals in place of
  // the reference's sixteen divides.
  float iL = 1.0 / min(SL - ustar, -1e-8), iR = 1.0 / max(SR - ustar, 1e-8);
  float iB = 1.0 / min(SB - vstar, -1e-8), iT = 1.0 / max(ST - vstar, 1e-8);
  float fLLx = (SL - uLL) * iL, fLLy = (SB - vLL) * iB;
  float fLRx = (SL - uLR) * iL, fLRy = (ST - vLR) * iT;
  float fRLx = (SR - uRL) * iR, fRLy = (SB - vRL) * iB;
  float fRRx = (SR - uRR) * iR, fRRy = (ST - vRR) * iT;

  // Per quadrant: the state behind the x fast wave (starred in x, so B is compressed
  // and the EMF uses ustar), the state behind the y fast wave (starred in y, A
  // compressed, vstar), and the state behind both. Three EMFs per quadrant because
  // the corner can sit in any of the three.
  float rstarLLx = rLL * fLLx, BstarLL = BLL * fLLx;
  float rstarLLy = rLL * fLLy, AstarLL = ALL * fLLy;
  float rstarLL  = rLL * fLLx * fLLy;
  float EstarLLx = ustar * BstarLL - vLL   * ALL;
  float EstarLLy = uLL   * BLL     - vstar * AstarLL;
  float EstarLL  = ustar * BstarLL - vstar * AstarLL;

  float rstarLRx = rLR * fLRx, BstarLR = BLR * fLRx;
  float rstarLRy = rLR * fLRy, AstarLR = ALR * fLRy;
  float rstarLR  = rLR * fLRx * fLRy;
  float EstarLRx = ustar * BstarLR - vLR   * ALR;
  float EstarLRy = uLR   * BLR     - vstar * AstarLR;
  float EstarLR  = ustar * BstarLR - vstar * AstarLR;

  float rstarRLx = rRL * fRLx, BstarRL = BRL * fRLx;
  float rstarRLy = rRL * fRLy, AstarRL = ARL * fRLy;
  float rstarRL  = rRL * fRLx * fRLy;
  float EstarRLx = ustar * BstarRL - vRL   * ARL;
  float EstarRLy = uRL   * BRL     - vstar * AstarRL;
  float EstarRL  = ustar * BstarRL - vstar * AstarRL;

  float rstarRRx = rRR * fRRx, BstarRR = BRR * fRRx;
  float rstarRRy = rRR * fRRy, AstarRR = ARR * fRRy;
  float rstarRR  = rRR * fRRx * fRRy;
  float EstarRRx = ustar * BstarRR - vRR   * ARR;
  float EstarRRy = uRR   * BRR     - vstar * AstarRR;
  float EstarRR  = ustar * BstarRR - vstar * AstarRR;

  // the twelve Alfven denominators, floored once each
  float qLLx = irs(rstarLLx), qLRx = irs(rstarLRx);
  float qRLx = irs(rstarRLx), qRRx = irs(rstarRRx);
  float qLLy = irs(rstarLLy), qLRy = irs(rstarLRy);
  float qRLy = irs(rstarRLy), qRRy = irs(rstarRRy);
  float qLL  = irs(rstarLL),  qLR  = irs(rstarLR);
  float qRL  = irs(rstarRL),  qRR  = irs(rstarRR);

  // The Alfven speed on each of the four sides of the corner, taken as the largest of
  // the singly- and doubly-starred estimates on that side. Both pairs are kept: the x
  // pair (L, R) and the y pair (B, T). smallc is the reference's floor and it is what
  // keeps the two Alfven waves from landing exactly on the contact when the normal
  // field vanishes -- with A = 0 the L/R pair collapses to +-smallc and the x part of
  // the selection degenerates gracefully to an average of the two adjacent EMFs
  // rather than to 0/0.
  float calfvenL = max(sc, max(max(abs(ALR) * qLRx, abs(AstarLR) * qLR),
                               max(abs(ALL) * qLLx, abs(AstarLL) * qLL)));
  float calfvenR = max(sc, max(max(abs(ARR) * qRRx, abs(AstarRR) * qRR),
                               max(abs(ARL) * qRLx, abs(AstarRL) * qRL)));
  float calfvenB = max(sc, max(max(abs(BLL) * qLLy, abs(BstarLL) * qLL),
                               max(abs(BRL) * qRLy, abs(BstarRL) * qRL)));
  float calfvenT = max(sc, max(max(abs(BLR) * qLRy, abs(BstarLR) * qLR),
                               max(abs(BRR) * qRRy, abs(BstarRR) * qRR)));

  // Clamped through zero, unlike the fast speeds: these four bound the region the
  // corner itself is in, so the reference wants SAL <= 0 <= SAR and SAB <= 0 <= SAT
  // whether or not the true Alfven waves straddle it.
  float SAL = min(ustar - calfvenL, 0.0), SAR = max(ustar + calfvenR, 0.0);
  float SAB = min(vstar - calfvenB, 0.0), SAT = max(vstar + calfvenT, 0.0);
  float iSAx = 1.0 / max(SAR - SAL, sc);
  float iSAy = 1.0 / max(SAT - SAB, sc);

  // the transverse field averaged across the Alfven pair of the *other* direction:
  // AstarT / AstarB are Bx on the top and bottom edges of the corner region, BstarR /
  // BstarL are By on its right and left edges. These are the two terms that carry the
  // rotation, and they exist only because both pairs were kept.
  float AstarT = (SAR * AstarRR - SAL * AstarLR) * iSAx;
  float AstarB = (SAR * AstarRL - SAL * AstarLL) * iSAx;
  float BstarR = (SAT * BstarRR - SAB * BstarRL) * iSAy;
  float BstarL = (SAT * BstarLR - SAB * BstarLL) * iSAy;

  // The selection. Read as a decision on where the corner sits: first whether the y
  // fan has left it behind entirely (SB > 0, everything moving up, so only the bottom
  // pair is upstream; ST < 0, only the top pair), then the same question in x. Where
  // one direction is supersonic and the other is not, the EMF is an HLL average over
  // that one direction's Alfven pair -- the middle three lines -- and where neither
  // is, the last expression averages all four doubly-starred EMFs over both pairs and
  // adds the two rotation terms. The sign in front of SAT*SAB differs from the one in
  // front of SAR*SAL because E is antisymmetric in the two directions: it is +u*By and
  // -v*Bx, so the jump in A enters the y average with the opposite sign to the jump in
  // B in the x average. Nine leaves, all of them the reference's.
  if (SB > 0.0) {
    if (SL > 0.0) return ELL;
    if (SR < 0.0) return ERL;
    return (SAR * EstarLLx - SAL * EstarRLx + SAR * SAL * (BRL - BLL)) * iSAx;
  }
  if (ST < 0.0) {
    if (SL > 0.0) return ELR;
    if (SR < 0.0) return ERR;
    return (SAR * EstarLRx - SAL * EstarRRx + SAR * SAL * (BRR - BLR)) * iSAx;
  }
  if (SL > 0.0) return (SAT * EstarLLy - SAB * EstarLRy - SAT * SAB * (ALR - ALL)) * iSAy;
  if (SR < 0.0) return (SAT * EstarRLy - SAB * EstarRRy - SAT * SAB * (ARR - ARL)) * iSAy;
  return (SAL * SAB * EstarRR - SAL * SAT * EstarRL
          - SAR * SAB * EstarLR + SAR * SAT * EstarLL) * iSAx * iSAy
         - SAT * SAB * iSAy * (AstarT - AstarB)
         + SAR * SAL * iSAx * (BstarR - BstarL);
}

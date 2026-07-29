// ref-riemann2d.mjs -- the ORACLE.
//
// A direct, line-by-line transcription of riemann2d_hlld from
//   mr/godunov_utils.f90  lines 1995-2212
// and of find_speed_fast from the same file, lines 1389-1425, into plain JS doubles.
// Written from the Fortran only. No guards, no floors, no reassociation: where the
// Fortran writes a/b/c this writes a/b/c, and where it divides by something that can
// be zero the result here is Infinity or NaN, exactly as the Fortran would trap.
//
// The two departures from the Fortran, both forced and both minimal:
//   * Bz (the reference's C, and CLL/CLR/CRL/CRR) is zero. It appears in only two
//     places, C*C in the four total pressures and C*C in B2 inside find_speed_fast,
//     so both terms are simply absent.
//   * find_speed_fast's c2 = gamma*P/d becomes cs2 when the EOS is isothermal, which
//     is the one substitution the isothermal mode needs anywhere in this routine.
// NENER is 0, so the radiation loops are absent. switch_llf_dmin / switch_llf_pmin
// are zero (the reference's defaults), so the LLF escape hatch never fires and is not
// transcribed.
//
// A quadrant q is {r, u, v, A, B, p}: A = Bx = qvar(6), B = By = qvar(7), p = qvar(5).
// params is {adia, gam, cs2, smallc}.

const half = 0.5;
const zero = 0.0;

// find_speed_fast(qvar) with qvar = (d, P, A, B, C=0). Returns cf.
// Exported under a second name so the test can aim a state at a chosen branch using
// the oracle's own wave speeds rather than an independent guess at them.
export { findSpeedFast as findSpeedFastRef };
function findSpeedFast(d, P, A, B, params) {
  const c2 = params.adia > 0.5 ? params.gam * P / d : params.cs2;
  const B2 = A * A + B * B;
  const d2 = half * (B2 / d + c2);
  return Math.sqrt(d2 + Math.sqrt(d2 * d2 - c2 * A * A / d));
}

// Returns { E, branch }. branch is 1..9 in the order the Fortran's leaves appear:
//   1 ELL, 2 ERL, 3 x-fan bottom, 4 ELR, 5 ERR, 6 x-fan top,
//   7 y-fan left, 8 y-fan right, 9 full two-dimensional star.
export function riemann2dRef(qLL, qLR, qRL, qRR, params) {
  const rLL = qLL.r, pLL = qLL.p, uLL = qLL.u, vLL = qLL.v, ALL = qLL.A, BLL = qLL.B;
  const rLR = qLR.r, pLR = qLR.p, uLR = qLR.u, vLR = qLR.v, ALR = qLR.A, BLR = qLR.B;
  const rRL = qRL.r, pRL = qRL.p, uRL = qRL.u, vRL = qRL.v, ARL = qRL.A, BRL = qRL.B;
  const rRR = qRR.r, pRR = qRR.p, uRR = qRR.u, vRR = qRR.v, ARR = qRR.A, BRR = qRR.B;

  const smallc = params.smallc;

  // Compute 4 fast magnetosonic velocity relative to x direction
  // qtmp(6)=A, qtmp(7)=B
  const cfastLLx = findSpeedFast(rLL, pLL, ALL, BLL, params);
  const cfastLRx = findSpeedFast(rLR, pLR, ALR, BLR, params);
  const cfastRLx = findSpeedFast(rRL, pRL, ARL, BRL, params);
  const cfastRRx = findSpeedFast(rRR, pRR, ARR, BRR, params);

  // Compute 4 fast magnetosonic velocity relative to y direction
  // qtmp(6)=B, qtmp(7)=A
  const cfastLLy = findSpeedFast(rLL, pLL, BLL, ALL, params);
  const cfastLRy = findSpeedFast(rLR, pLR, BLR, ALR, params);
  const cfastRLy = findSpeedFast(rRL, pRL, BRL, ARL, params);
  const cfastRRy = findSpeedFast(rRR, pRR, BRR, ARR, params);

  const SL = Math.min(uLL, uLR, uRL, uRR) - Math.max(cfastLLx, cfastLRx, cfastRLx, cfastRRx);
  const SR = Math.max(uLL, uLR, uRL, uRR) + Math.max(cfastLLx, cfastLRx, cfastRLx, cfastRRx);
  const SB = Math.min(vLL, vLR, vRL, vRR) - Math.max(cfastLLy, cfastLRy, cfastRLy, cfastRRy);
  const ST = Math.max(vLL, vLR, vRL, vRR) + Math.max(cfastLLy, cfastLRy, cfastRLy, cfastRRy);

  const ELL = uLL * BLL - vLL * ALL;
  const ELR = uLR * BLR - vLR * ALR;
  const ERL = uRL * BRL - vRL * ARL;
  const ERR = uRR * BRR - vRR * ARR;

  // C = 0 in every quadrant, so C*C is gone from all four
  const PtotLL = pLL + half * (ALL * ALL + BLL * BLL);
  const PtotLR = pLR + half * (ALR * ALR + BLR * BLR);
  const PtotRL = pRL + half * (ARL * ARL + BRL * BRL);
  const PtotRR = pRR + half * (ARR * ARR + BRR * BRR);

  const rcLLx = rLL * (uLL - SL), rcRLx = rRL * (SR - uRL);
  const rcLRx = rLR * (uLR - SL), rcRRx = rRR * (SR - uRR);
  const rcLLy = rLL * (vLL - SB), rcLRy = rLR * (ST - vLR);
  const rcRLy = rRL * (vRL - SB), rcRRy = rRR * (ST - vRR);

  const ustar = (rcLLx * uLL + rcLRx * uLR + rcRLx * uRL + rcRRx * uRR
    + (PtotLL - PtotRL + PtotLR - PtotRR)) / (rcLLx + rcLRx + rcRLx + rcRRx);
  const vstar = (rcLLy * vLL + rcLRy * vLR + rcRLy * vRL + rcRRy * vRR
    + (PtotLL - PtotLR + PtotRL - PtotRR)) / (rcLLy + rcLRy + rcRLy + rcRRy);

  const rstarLLx = rLL * (SL - uLL) / (SL - ustar), BstarLL = BLL * (SL - uLL) / (SL - ustar);
  const rstarLLy = rLL * (SB - vLL) / (SB - vstar), AstarLL = ALL * (SB - vLL) / (SB - vstar);
  const rstarLL = rLL * (SL - uLL) / (SL - ustar) * (SB - vLL) / (SB - vstar);
  const EstarLLx = ustar * BstarLL - vLL * ALL;
  const EstarLLy = uLL * BLL - vstar * AstarLL;
  const EstarLL = ustar * BstarLL - vstar * AstarLL;

  const rstarLRx = rLR * (SL - uLR) / (SL - ustar), BstarLR = BLR * (SL - uLR) / (SL - ustar);
  const rstarLRy = rLR * (ST - vLR) / (ST - vstar), AstarLR = ALR * (ST - vLR) / (ST - vstar);
  const rstarLR = rLR * (SL - uLR) / (SL - ustar) * (ST - vLR) / (ST - vstar);
  const EstarLRx = ustar * BstarLR - vLR * ALR;
  const EstarLRy = uLR * BLR - vstar * AstarLR;
  const EstarLR = ustar * BstarLR - vstar * AstarLR;

  const rstarRLx = rRL * (SR - uRL) / (SR - ustar), BstarRL = BRL * (SR - uRL) / (SR - ustar);
  const rstarRLy = rRL * (SB - vRL) / (SB - vstar), AstarRL = ARL * (SB - vRL) / (SB - vstar);
  const rstarRL = rRL * (SR - uRL) / (SR - ustar) * (SB - vRL) / (SB - vstar);
  const EstarRLx = ustar * BstarRL - vRL * ARL;
  const EstarRLy = uRL * BRL - vstar * AstarRL;
  const EstarRL = ustar * BstarRL - vstar * AstarRL;

  const rstarRRx = rRR * (SR - uRR) / (SR - ustar), BstarRR = BRR * (SR - uRR) / (SR - ustar);
  const rstarRRy = rRR * (ST - vRR) / (ST - vstar), AstarRR = ARR * (ST - vRR) / (ST - vstar);
  const rstarRR = rRR * (SR - uRR) / (SR - ustar) * (ST - vRR) / (ST - vstar);
  const EstarRRx = ustar * BstarRR - vRR * ARR;
  const EstarRRy = uRR * BRR - vstar * AstarRR;
  const EstarRR = ustar * BstarRR - vstar * AstarRR;

  const calfvenL = Math.max(Math.abs(ALR) / Math.sqrt(rstarLRx), Math.abs(AstarLR) / Math.sqrt(rstarLR),
    Math.abs(ALL) / Math.sqrt(rstarLLx), Math.abs(AstarLL) / Math.sqrt(rstarLL), smallc);
  const calfvenR = Math.max(Math.abs(ARR) / Math.sqrt(rstarRRx), Math.abs(AstarRR) / Math.sqrt(rstarRR),
    Math.abs(ARL) / Math.sqrt(rstarRLx), Math.abs(AstarRL) / Math.sqrt(rstarRL), smallc);
  const calfvenB = Math.max(Math.abs(BLL) / Math.sqrt(rstarLLy), Math.abs(BstarLL) / Math.sqrt(rstarLL),
    Math.abs(BRL) / Math.sqrt(rstarRLy), Math.abs(BstarRL) / Math.sqrt(rstarRL), smallc);
  const calfvenT = Math.max(Math.abs(BLR) / Math.sqrt(rstarLRy), Math.abs(BstarLR) / Math.sqrt(rstarLR),
    Math.abs(BRR) / Math.sqrt(rstarRRy), Math.abs(BstarRR) / Math.sqrt(rstarRR), smallc);

  const SAL = Math.min(ustar - calfvenL, zero), SAR = Math.max(ustar + calfvenR, zero);
  const SAB = Math.min(vstar - calfvenB, zero), SAT = Math.max(vstar + calfvenT, zero);
  const AstarT = (SAR * AstarRR - SAL * AstarLR) / (SAR - SAL);
  const AstarB = (SAR * AstarRL - SAL * AstarLL) / (SAR - SAL);
  const BstarR = (SAT * BstarRR - SAB * BstarRL) / (SAT - SAB);
  const BstarL = (SAT * BstarLR - SAB * BstarLL) / (SAT - SAB);

  let E, branch;
  if (SB > 0.0) {
    if (SL > 0.0) {
      E = ELL; branch = 1;
    } else if (SR < 0.0) {
      E = ERL; branch = 2;
    } else {
      E = (SAR * EstarLLx - SAL * EstarRLx + SAR * SAL * (BRL - BLL)) / (SAR - SAL); branch = 3;
    }
  } else if (ST < 0.0) {
    if (SL > 0.0) {
      E = ELR; branch = 4;
    } else if (SR < 0.0) {
      E = ERR; branch = 5;
    } else {
      E = (SAR * EstarLRx - SAL * EstarRRx + SAR * SAL * (BRR - BLR)) / (SAR - SAL); branch = 6;
    }
  } else if (SL > 0.0) {
    E = (SAT * EstarLLy - SAB * EstarLRy - SAT * SAB * (ALR - ALL)) / (SAT - SAB); branch = 7;
  } else if (SR < 0.0) {
    E = (SAT * EstarRLy - SAB * EstarRRy - SAT * SAB * (ARR - ARL)) / (SAT - SAB); branch = 8;
  } else {
    E = (SAL * SAB * EstarRR - SAL * SAT * EstarRL - SAR * SAB * EstarLR + SAR * SAT * EstarLL) / (SAR - SAL) / (SAT - SAB)
      - SAT * SAB / (SAT - SAB) * (AstarT - AstarB) + SAR * SAL / (SAR - SAL) * (BstarR - BstarL);
    branch = 9;
  }

  // diagnostics the test uses to classify conditioning and to see whether any of the
  // mirror's guards would have bound on this state; not part of the algorithm
  const minRstar = Math.min(rstarLLx, rstarLRx, rstarRLx, rstarRRx,
    rstarLLy, rstarLRy, rstarRLy, rstarRRy, rstarLL, rstarLR, rstarRL, rstarRR);
  return {
    E, branch, SL, SR, SB, ST, SAL, SAR, SAB, SAT, ustar, vstar, minRstar,
    sumRcx: rcLLx + rcLRx + rcRLx + rcRRx, sumRcy: rcLLy + rcLRy + rcRLy + rcRRy
  };
}

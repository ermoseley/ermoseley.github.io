// mirror-riemann2d.mjs -- a mechanical transcription of piece-riemann2d.glsl.
//
// Same expressions, same order, same guards, same names. GLSL's min/max/abs/sqrt
// become Math.min/Math.max/Math.abs/Math.sqrt; Prim's h/f become plain objects. The
// only additions are the branch tag the test uses for coverage and the diagnostics at
// the bottom. This file is what the test measures; it is not the oracle.
//
// A quadrant is {h:{x:rho,y:u,z:v,w:Y}, f:{x:Bx,y:By,z:psi,w:p}} -- the engine's Prim.
// The uniforms the shader reads (cs2, adia, gam, smallr) come in as a `U` object here
// because JS has no uniforms; in the shader they are global and the signature is
//   float riemann2d(Prim LL, Prim LR, Prim RL, Prim RR, float smallc)

function sound2(r, p, U) { return U.adia > 0.5 ? U.gam * p / r : U.cs2; }

function cfast2(r, p, An, At, U) {
  const c2 = Math.max(sound2(r, p, U), 0.0);
  const d2 = 0.5 * ((An * An + At * At) / r + c2);
  return Math.sqrt(d2 + Math.sqrt(Math.max(d2 * d2 - c2 * An * An / r, 0.0)));
}

function irs(rs, U) { return 1.0 / Math.sqrt(Math.max(rs, U.smallr)); }

export function riemann2dMirror(LL, LR, RL, RR, smallc, U) {
  const rLL = Math.max(LL.h.x, U.smallr), rLR = Math.max(LR.h.x, U.smallr);
  const rRL = Math.max(RL.h.x, U.smallr), rRR = Math.max(RR.h.x, U.smallr);
  const uLL = LL.h.y, uLR = LR.h.y, uRL = RL.h.y, uRR = RR.h.y;
  const vLL = LL.h.z, vLR = LR.h.z, vRL = RL.h.z, vRR = RR.h.z;
  const ALL = LL.f.x, ALR = LR.f.x, ARL = RL.f.x, ARR = RR.f.x;   // A = Bx
  const BLL = LL.f.y, BLR = LR.f.y, BRL = RL.f.y, BRR = RR.f.y;   // B = By
  const pLL = LL.f.w, pLR = LR.f.w, pRL = RL.f.w, pRR = RR.f.w;   // gas pressure
  const sc = Math.max(smallc, 1e-20);

  const cx = Math.max(Math.max(cfast2(rLL, pLL, ALL, BLL, U), cfast2(rLR, pLR, ALR, BLR, U)),
    Math.max(cfast2(rRL, pRL, ARL, BRL, U), cfast2(rRR, pRR, ARR, BRR, U)));
  const cy = Math.max(Math.max(cfast2(rLL, pLL, BLL, ALL, U), cfast2(rLR, pLR, BLR, ALR, U)),
    Math.max(cfast2(rRL, pRL, BRL, ARL, U), cfast2(rRR, pRR, BRR, ARR, U)));

  const SL = Math.min(Math.min(uLL, uLR), Math.min(uRL, uRR)) - cx;
  const SR = Math.max(Math.max(uLL, uLR), Math.max(uRL, uRR)) + cx;
  const SB = Math.min(Math.min(vLL, vLR), Math.min(vRL, vRR)) - cy;
  const ST = Math.max(Math.max(vLL, vLR), Math.max(vRL, vRR)) + cy;

  const ELL = uLL * BLL - vLL * ALL;
  const ELR = uLR * BLR - vLR * ALR;
  const ERL = uRL * BRL - vRL * ARL;
  const ERR = uRR * BRR - vRR * ARR;

  const PtotLL = pLL + 0.5 * (ALL * ALL + BLL * BLL);
  const PtotLR = pLR + 0.5 * (ALR * ALR + BLR * BLR);
  const PtotRL = pRL + 0.5 * (ARL * ARL + BRL * BRL);
  const PtotRR = pRR + 0.5 * (ARR * ARR + BRR * BRR);

  const rcLLx = rLL * (uLL - SL), rcRLx = rRL * (SR - uRL);
  const rcLRx = rLR * (uLR - SL), rcRRx = rRR * (SR - uRR);
  const rcLLy = rLL * (vLL - SB), rcLRy = rLR * (ST - vLR);
  const rcRLy = rRL * (vRL - SB), rcRRy = rRR * (ST - vRR);

  const iux = 1.0 / Math.max(rcLLx + rcLRx + rcRLx + rcRRx, 1e-20);
  const ivy = 1.0 / Math.max(rcLLy + rcLRy + rcRLy + rcRRy, 1e-20);
  const ustar = (rcLLx * uLL + rcLRx * uLR + rcRLx * uRL + rcRRx * uRR
    + (PtotLL - PtotRL + PtotLR - PtotRR)) * iux;
  const vstar = (rcLLy * vLL + rcLRy * vLR + rcRLy * vRL + rcRRy * vRR
    + (PtotLL - PtotLR + PtotRL - PtotRR)) * ivy;

  const iL = 1.0 / Math.min(SL - ustar, -1e-8), iR = 1.0 / Math.max(SR - ustar, 1e-8);
  const iB = 1.0 / Math.min(SB - vstar, -1e-8), iT = 1.0 / Math.max(ST - vstar, 1e-8);
  const fLLx = (SL - uLL) * iL, fLLy = (SB - vLL) * iB;
  const fLRx = (SL - uLR) * iL, fLRy = (ST - vLR) * iT;
  const fRLx = (SR - uRL) * iR, fRLy = (SB - vRL) * iB;
  const fRRx = (SR - uRR) * iR, fRRy = (ST - vRR) * iT;

  const rstarLLx = rLL * fLLx, BstarLL = BLL * fLLx;
  const rstarLLy = rLL * fLLy, AstarLL = ALL * fLLy;
  const rstarLL = rLL * fLLx * fLLy;
  const EstarLLx = ustar * BstarLL - vLL * ALL;
  const EstarLLy = uLL * BLL - vstar * AstarLL;
  const EstarLL = ustar * BstarLL - vstar * AstarLL;

  const rstarLRx = rLR * fLRx, BstarLR = BLR * fLRx;
  const rstarLRy = rLR * fLRy, AstarLR = ALR * fLRy;
  const rstarLR = rLR * fLRx * fLRy;
  const EstarLRx = ustar * BstarLR - vLR * ALR;
  const EstarLRy = uLR * BLR - vstar * AstarLR;
  const EstarLR = ustar * BstarLR - vstar * AstarLR;

  const rstarRLx = rRL * fRLx, BstarRL = BRL * fRLx;
  const rstarRLy = rRL * fRLy, AstarRL = ARL * fRLy;
  const rstarRL = rRL * fRLx * fRLy;
  const EstarRLx = ustar * BstarRL - vRL * ARL;
  const EstarRLy = uRL * BRL - vstar * AstarRL;
  const EstarRL = ustar * BstarRL - vstar * AstarRL;

  const rstarRRx = rRR * fRRx, BstarRR = BRR * fRRx;
  const rstarRRy = rRR * fRRy, AstarRR = ARR * fRRy;
  const rstarRR = rRR * fRRx * fRRy;
  const EstarRRx = ustar * BstarRR - vRR * ARR;
  const EstarRRy = uRR * BRR - vstar * AstarRR;
  const EstarRR = ustar * BstarRR - vstar * AstarRR;

  const qLLx = irs(rstarLLx, U), qLRx = irs(rstarLRx, U);
  const qRLx = irs(rstarRLx, U), qRRx = irs(rstarRRx, U);
  const qLLy = irs(rstarLLy, U), qLRy = irs(rstarLRy, U);
  const qRLy = irs(rstarRLy, U), qRRy = irs(rstarRRy, U);
  const qLL = irs(rstarLL, U), qLR = irs(rstarLR, U);
  const qRL = irs(rstarRL, U), qRR = irs(rstarRR, U);

  const calfvenL = Math.max(sc, Math.max(Math.max(Math.abs(ALR) * qLRx, Math.abs(AstarLR) * qLR),
    Math.max(Math.abs(ALL) * qLLx, Math.abs(AstarLL) * qLL)));
  const calfvenR = Math.max(sc, Math.max(Math.max(Math.abs(ARR) * qRRx, Math.abs(AstarRR) * qRR),
    Math.max(Math.abs(ARL) * qRLx, Math.abs(AstarRL) * qRL)));
  const calfvenB = Math.max(sc, Math.max(Math.max(Math.abs(BLL) * qLLy, Math.abs(BstarLL) * qLL),
    Math.max(Math.abs(BRL) * qRLy, Math.abs(BstarRL) * qRL)));
  const calfvenT = Math.max(sc, Math.max(Math.max(Math.abs(BLR) * qLRy, Math.abs(BstarLR) * qLR),
    Math.max(Math.abs(BRR) * qRRy, Math.abs(BstarRR) * qRR)));

  const SAL = Math.min(ustar - calfvenL, 0.0), SAR = Math.max(ustar + calfvenR, 0.0);
  const SAB = Math.min(vstar - calfvenB, 0.0), SAT = Math.max(vstar + calfvenT, 0.0);
  const iSAx = 1.0 / Math.max(SAR - SAL, sc);
  const iSAy = 1.0 / Math.max(SAT - SAB, sc);

  const AstarT = (SAR * AstarRR - SAL * AstarLR) * iSAx;
  const AstarB = (SAR * AstarRL - SAL * AstarLL) * iSAx;
  const BstarR = (SAT * BstarRR - SAB * BstarRL) * iSAy;
  const BstarL = (SAT * BstarLR - SAB * BstarLL) * iSAy;

  const diag = { SL, SR, SB, ST, SAL, SAR, SAB, SAT, ustar, vstar, cx, cy };

  if (SB > 0.0) {
    if (SL > 0.0) return { E: ELL, branch: 1, diag };
    if (SR < 0.0) return { E: ERL, branch: 2, diag };
    return { E: (SAR * EstarLLx - SAL * EstarRLx + SAR * SAL * (BRL - BLL)) * iSAx, branch: 3, diag };
  }
  if (ST < 0.0) {
    if (SL > 0.0) return { E: ELR, branch: 4, diag };
    if (SR < 0.0) return { E: ERR, branch: 5, diag };
    return { E: (SAR * EstarLRx - SAL * EstarRRx + SAR * SAL * (BRR - BLR)) * iSAx, branch: 6, diag };
  }
  if (SL > 0.0) return { E: (SAT * EstarLLy - SAB * EstarLRy - SAT * SAB * (ALR - ALL)) * iSAy, branch: 7, diag };
  if (SR < 0.0) return { E: (SAT * EstarRLy - SAB * EstarRRy - SAT * SAB * (ARR - ARL)) * iSAy, branch: 8, diag };
  return {
    E: (SAL * SAB * EstarRR - SAL * SAT * EstarRL
      - SAR * SAB * EstarLR + SAR * SAT * EstarLL) * iSAx * iSAy
      - SAT * SAB * iSAy * (AstarT - AstarB)
      + SAR * SAL * iSAx * (BstarR - BstarL),
    branch: 9, diag
  };
}

#ifndef HARDENINGSTATE_HPP_INCLUDED
#define HARDENINGSTATE_HPP_INCLUDED

#include "bulkStateField.hpp"

/*
 * Strain hardening and recovery of dislocation creep (docs/REVERSIBLE_STRENGTH_PLAN.md, 4.1).
 * S in [0, 1] (0 recovered, 1 fully hardened):
 *
 *     dS/dt = (e/gamma_h) (1 - S) - S/tau_r(T),   tau_r(T) = tau_r0 exp(Q_r/R (1/T - 1/T_ref)),
 *
 * with e the dislocation strain rate (1/s), and the strain rate of dislocation creep at a given
 * stress divided by H(S)^n, H(S) = 1 + a (S - S_ref) (a factor H^-n on its inverse viscosity). At
 * S = S_ref the calibrated flow law is recovered exactly. S_ref is hard_SRef, or the steady state
 * x/(1 + x), x = e_ref tau_r(T)/gamma_h, at the reference strain rate hard_eRef (per node, from the
 * initial temperature). Steady state at rate e: S_ss = x/(1 + x), x = e tau_r/gamma_h.
 *
 * Keys: hard_type, hard_a, hard_gammaHVals/Depths, hard_tauR0 (s), hard_QR (K), hard_TRef (K),
 * hard_SRef or hard_eRef (1/s), initial hard_SVals/Depths or the file hard_S; and the common
 * hard_eTest... of BulkStateField. Needs wDislCreep = yes.
 */
class HardeningState : public BulkStateField
{
public:
  PetscScalar         _a, _tauR0, _QR, _TRef, _SRef, _eRef;
  std::vector<double> _gammaHVals, _gammaHDepths;
  Vec                 _gammaH, _SRefVec, _H;

  HardeningState(Domain& D, const Vec& T0); // T0: the temperature for S_ref from hard_eRef
  ~HardeningState();

  PetscErrorCode computeRate(const BulkInputs& in);
  PetscErrorCode pushToMaterial(PowerLaw& material);
  PetscErrorCode computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT);
  PetscErrorCode writeContextExtra(PetscViewer& ascii, PetscViewer& viewer);
  PetscErrorCode clamp();
  PetscErrorCode initialDefault(Vec& state); // S_ref
  PetscScalar tauR(const PetscScalar T) const { return _tauR0*exp(_QR*(1.0/T - 1.0/_TRef)); }
};

#endif

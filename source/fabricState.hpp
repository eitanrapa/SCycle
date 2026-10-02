#ifndef FABRICSTATE_HPP_INCLUDED
#define FABRICSTATE_HPP_INCLUDED

#include "bulkStateField.hpp"

/*
 * Fabric of the shear zone (docs/REVERSIBLE_STRENGTH_PLAN.md, 4.5), a weakening feedback. Phi in [0, 1]:
 *
 *     dPhi/dt = (e/gamma_f) (1 - Phi) [e < e_crit] - Phi/tau_c(T),   tau_c(T) = tau_c0 exp(Q_c (1/T - 1/T_ref)),
 *
 * with e the total viscous strain rate (1/s) and the bracket a gate (fabric_eCrit > 0: fabric stops
 * growing above e_crit, where recrystallization rather than planar fabric is expected). The strength
 * of every creep mechanism is multiplied by 1/(1 + beta_f Phi) (PowerLaw::setStrengthFactor: at a
 * given stress its strain rate by (1 + beta_f Phi)^n). Steady state at rate e (no gate):
 * Phi_ss = x/(1 + x), x = e tau_c/gamma_f.
 *
 * Keys: fabric_type, fabric_betaF (default 4), fabric_gammaFVals/Depths, fabric_tauC0 (s),
 * fabric_QC (K, default 0), fabric_TRef (K), fabric_eCrit (1/s, default 0: no gate), initial
 * fabric_PhiVals/Depths or the file fabric_Phi (default 0); and the common fabric_eTest...
 */
class FabricState : public BulkStateField
{
public:
  PetscScalar         _betaF, _tauC0, _QC, _TRef, _eCrit;
  std::vector<double> _gammaFVals, _gammaFDepths;
  Vec                 _gammaF, _F;

  FabricState(Domain& D);
  ~FabricState();

  bool usesTotalViscousRate() const { return true; }
  PetscErrorCode computeRate(const BulkInputs& in);
  PetscErrorCode pushToMaterial(PowerLaw& material);
  PetscErrorCode computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT);
  PetscErrorCode writeContextExtra(PetscViewer& ascii, PetscViewer& viewer);
  PetscErrorCode clamp();
  PetscScalar tauC(const PetscScalar T) const { return _tauC0*exp(_QC*(1.0/T - 1.0/_TRef)); }
};

#endif

#ifndef CEMENTSTATE_HPP_INCLUDED
#define CEMENTSTATE_HPP_INCLUDED

#include "bulkStateField.hpp"

/*
 * Cement of the shear zone (docs/REVERSIBLE_STRENGTH_PLAN.md, 4.5), a memory driven by earthquakes.
 * Psi in [0, 1]:
 *
 *     dPsi/dt = (Q/E_c) (1 - Psi) - Psi (1/tau_Psi + e/gamma_Psi),
 *
 * with Q the faults' work spread into the body over the frictional-heat kernels (kW/m^3) and e the
 * total viscous strain rate (1/s). A slip delta at traction tau deposits tau delta Gw (kJ/m^3, Gw the
 * kernel in 1/km), so E_c = tau D_cem Gw(0) makes D_cem of slip at tau cement the core of the zone.
 * The strength of every creep mechanism is multiplied by 1 + beta_c Psi.
 *
 * Keys: cement_type, cement_betaC (default 0.5), cement_Ec (kJ/m^3), cement_tauPsi (s),
 * cement_gammaPsi (default 0: no removal by deformation), cement_QTest (kW/m^3; tests: a uniform Q in
 * place of the faults' work), initial cement_PsiVals/Depths or the file cement_Psi (default 0); and
 * the common cement_eTest...
 */
class CementState : public BulkStateField
{
public:
  PetscScalar _betaC, _Ec, _tauPsi, _gammaPsi, _QTest;
  Vec         _F;

  CementState(Domain& D);
  ~CementState();

  bool usesTotalViscousRate() const { return true; }
  bool needsFaultWork() const { return _QTest < 0 && _type == "transient"; }
  PetscErrorCode computeRate(const BulkInputs& in);
  PetscErrorCode pushToMaterial(PowerLaw& material);
  PetscErrorCode computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT);
  PetscErrorCode writeContextExtra(PetscViewer& ascii, PetscViewer& viewer);
  PetscErrorCode clamp();
};

#endif

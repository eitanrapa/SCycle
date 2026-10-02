#ifndef WATERSTATE_HPP_INCLUDED
#define WATERSTATE_HPP_INCLUDED

#include "bulkStateField.hpp"

/*
 * Water content of the shear zone (docs/REVERSIBLE_STRENGTH_PLAN.md, 4.2). chi in [chi_min, 1]:
 *
 *     dchi/dt = (chi_sup - chi)/tau_hyd (1 + kappa (1 - chi)) - (e/gamma_dry) (chi - chi_min),
 *
 * with e the dislocation strain rate (1/s): supply relaxes chi toward chi_sup, deformation dries it
 * toward chi_min. It mixes the dry (disl_) and wet (disl2_) dislocation creep laws in the power law
 * (wDislWetDry = yes): log-linear, 1/eta = (1/eta_wet)^chi (1/eta_dry)^(1 - chi), or arithmetic,
 * chi/eta_wet + (1 - chi)/eta_dry (water_mix); and pressure solution (wDissPrecCreep) scales with
 * the wetted fraction max(0, (chi - chi_c)/(1 - chi_c)). Steady state at rate e:
 * chi_ss = chi_min + (chi_sup - chi_min)/(1 + x), x = tau_hyd e/gamma_dry (kappa = 0).
 *
 * Keys: water_type, water_chiSupVals/Depths (default 1), water_chiMin, water_tauHyd (s),
 * water_kappa, water_gammaDry, water_chiC, water_mix = log | arithmetic, water_supplyFrom = profile
 * (pressure and pulses are not yet implemented), initial water_chiVals/Depths or the file water_chi
 * (by default the power law's wetness: the file pl_wetDist, or 1); and the common water_eTest...
 */
class WaterState : public BulkStateField
{
public:
  std::vector<double> _chiSupVals, _chiSupDepths;
  Vec                 _chiSup;
  PetscScalar         _chiMin, _tauHyd, _kappa, _gammaDry, _chiC;
  std::string         _mix, _supplyFrom;
  Vec                 _wetDist0; // default initial state (a copy of the power law's wetDist)

  WaterState(Domain& D, const Vec& wetDist0);
  ~WaterState();

  PetscErrorCode computeRate(const BulkInputs& in);
  PetscErrorCode pushToMaterial(PowerLaw& material);
  PetscErrorCode computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT);
  PetscErrorCode writeContextExtra(PetscViewer& ascii, PetscViewer& viewer);
  PetscErrorCode clamp();
  PetscErrorCode initialDefault(Vec& state);
};

#endif

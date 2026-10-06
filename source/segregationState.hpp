#ifndef SEGREGATIONSTATE_HPP_INCLUDED
#define SEGREGATIONSTATE_HPP_INCLUDED

#include "bulkStateField.hpp"

/*
 * Phase segregation of the shear zone (Cawood & Dolan 2024, S10; docs/REVERSIBLE_STRENGTH_PLAN.md),
 * a strengthening memory of pressure solution. During dissolution-precipitation creep quartz dissolves at
 * high-stress sites and precipitates in low-stress ones, so a quartz-mica rock separates into quartz- and
 * mica-rich bands; pressure solution is fastest along quartz-mica grain boundaries, which segregation
 * removes, so it slows down (Schmidt & Platt 2022). Xi in [0, 1] (0 mixed, 1 segregated):
 *
 *     dXi/dt = (e_dp/gamma_s) (1 - Xi) - Xi (1/tau_m(T) + e_disl/gamma_m),   tau_m(T) = tau_m0 exp(Q_m (1/T - 1/T_ref)),
 *
 * with e_dp the pressure-solution strain rate (1/s), which segregates over the strain gamma_s, and two ways
 * of remixing: in time over tau_m, and by dislocation creep's strain e_disl over gamma_m (seg_gammaM = 0:
 * none). At a given stress the pressure-solution strain rate is divided by 1 + beta_s Xi
 * (PowerLaw::setDPRateFactor); the other mechanisms are unchanged. Steady state at constant rates:
 * Xi_ss = x/(1 + x + y), x = e_dp tau_m/gamma_s, y = e_disl tau_m/gamma_m.
 *
 * Needs wDissPrecCreep = yes. Keys: seg_type, seg_betaS (default 1), seg_gammaSVals/Depths,
 * seg_tauM0 (s), seg_QM (K, default 0), seg_TRef (K), seg_gammaM (default 0), seg_eMixTest (1/s, tests:
 * a prescribed e_disl; off when negative), initial seg_XiVals/Depths or the file seg_Xi (default 0); and
 * the common seg_eTest... (a prescribed e_dp).
 */
class SegregationState : public BulkStateField
{
public:
  PetscScalar         _betaS, _tauM0, _QM, _TRef, _gammaM, _eMixTest;
  std::vector<double> _gammaSVals, _gammaSDepths;
  Vec                 _gammaS, _f;

  SegregationState(Domain& D);
  ~SegregationState();

  Vec drivingRate(const BulkInputs& in) const { return in.dgVdev_dp; }
  PetscErrorCode computeRate(const BulkInputs& in);
  PetscErrorCode pushToMaterial(PowerLaw& material);
  PetscErrorCode computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT);
  PetscErrorCode writeContextExtra(PetscViewer& ascii, PetscViewer& viewer);
  PetscErrorCode clamp();
  PetscScalar tauM(const PetscScalar T) const { return _tauM0*exp(_QM*(1.0/T - 1.0/_TRef)); }
  // the remixing rate by dislocation creep (1/s) at node Jj: the prescribed test rate or e_disl/gamma_m
  PetscScalar mixRate(const PetscScalar edisl) const { return (_gammaM > 0) ? ((_eMixTest >= 0) ? _eMixTest : edisl)/_gammaM : 0.0; }
};

#endif

#ifndef PSEUDOTACHYLITESTATE_HPP_INCLUDED
#define PSEUDOTACHYLITESTATE_HPP_INCLUDED

#include "bulkStateField.hpp"

/*
 * Pseudotachylite products (Cawood & Dolan 2024, W7; docs/REVERSIBLE_STRENGTH_PLAN.md 4.9), a weakening of
 * the rock beside the faults by their seismic slip. Frictional work done at slip rates of at least v_melt
 * melts the volume Q_melt/e_melt per unit time (e_melt: the heat to bring a unit volume to melting and melt
 * it), which quenches to ultrafine grains of size d_melt; the products anneal (devitrify and coarsen into
 * the host) over tau_a. phi in [0, 1] is their volume fraction:
 *
 *     dphi/dt = (Q_melt/e_melt) (1 - phi) - phi/tau_a(T),   tau_a(T) = tau_a0 exp(Q_a (1/T - 1/T_ref)),
 *
 * with Q_melt the faults' work at |V| >= v_melt spread into the body over the frictional-heat kernels
 * (kW/m^3, so it needs wVals > 0). The grain-size-sensitive mechanisms creep as the mixture of host
 * grains d (the power law's grain size) and products, d_eff^-m = (1 - phi) d^-m + phi d_melt^-m: at a
 * given stress diffusion creep and pressure solution are multiplied by (1 - phi) + phi (d/d_melt)^m with
 * their own exponent m (PowerLaw::setMeltProducts). This leaves the host grain size (and its evolution)
 * alone and is linear in phi, where folding the products into the cell's grain size would be stiff.
 * Steady state under steady melting: phi_ss = x/(1 + x), x = Q_melt tau_a/e_melt.
 *
 * Needs wDiffCreep = yes or wDissPrecCreep = yes. Keys: pt_type, pt_dMeltVals/Depths (grain-size units),
 * pt_eMelt (kJ/m^3, default 3.6e6: 2700 kg/m^3 times 1 kJ/kg/K over 1000 K plus 330 kJ/kg), pt_vMelt
 * (m/s, default 0.1), pt_tauA0 (s), pt_QA (K, default 0), pt_TRef (K), pt_QTest (kW/m^3, tests: a uniform
 * Q_melt; off when negative), initial pt_phiPTVals/Depths or the file pt_phiPT (default 0).
 */
class PseudotachyliteState : public BulkStateField
{
public:
  PetscScalar         _eMelt, _vMelt, _tauA0, _QA, _TRef, _QTest;
  std::vector<double> _dMeltVals, _dMeltDepths;
  Vec                 _dMelt;

  PseudotachyliteState(Domain& D);
  ~PseudotachyliteState();

  PetscScalar meltSlipRate() const { return _vMelt; }
  PetscErrorCode computeRate(const BulkInputs& in);
  PetscErrorCode pushToMaterial(PowerLaw& material);
  PetscErrorCode computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT);
  PetscErrorCode writeContextExtra(PetscViewer& ascii, PetscViewer& viewer);
  PetscErrorCode clamp();
  PetscScalar tauA(const PetscScalar T) const { return _tauA0*exp(_QA*(1.0/T - 1.0/_TRef)); }
};

#endif

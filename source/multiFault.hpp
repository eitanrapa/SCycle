#ifndef MULTIFAULT_HPP_INCLUDED
#define MULTIFAULT_HPP_INCLUDED

#include <petscksp.h>
#include <petscviewerhdf5.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include "domain.hpp"
#include "fault.hpp"
#include "interiorFaultLift.hpp"
#include "sbpOps.hpp"

/*
 * Helpers shared by the quasi-dynamic mediators for models with several faults
 * (docs/TWO_FAULT_DESIGN.md): creating interior faults, step-size control over every fault, the
 * per-step fault series, and the fault positions in the context files.
 */

// Creates the interior faults named in names, each midway between the two grid rows around its
// position <name>_y (from positions), with radiation damping sqrt(mu rho)/2, and appends each
// fault and its lift to faults and lifts. Prints a note when a fault moves to the midpoint.
PetscErrorCode createInteriorFaults(Domain& D, const std::vector<std::string>& names,
  std::map<std::string,PetscScalar>& positions, const PetscScalar faultTypeScale, const bool kinkLift,
  std::vector<Fault_qd*>& faults, std::vector<InteriorFaultLift*>& lifts);

// A timeIntInds that lists slip or psi means those of every fault: <name>_slip and <name>_psi are
// added with the same scale, and slip or psi dropped if no fault is named fault. An empty list
// (every integrated variable) is left alone. Prints the final list when it changed.
PetscErrorCode extendTimeIntInds(std::vector<std::string>& inds, std::vector<double>& scale,
  const std::vector<Fault_qd*>& faults);

// mediator.txt lines (interiorFaults, interiorFaultKinkLift, <name>_y) and, in data_context.h5, the
// attributes y (km), iRow and dy (row spacing across the fault, km) of each interior fault's group.
PetscErrorCode printInteriorFaults(PetscViewer& ascii, const std::vector<std::string>& names, const int kinkLift,
  const std::vector<Fault_qd*>& faults, const std::vector<InteriorFaultLift*>& lifts);
PetscErrorCode writeInteriorFaultContext(PetscViewer& viewer, const std::vector<Fault_qd*>& faults,
  const std::vector<InteriorFaultLift*>& lifts);

// faultSeries.txt: one line per call with step, time (s), time step (s), then per fault the maximum
// slip rate (m/s), its depth (km), the potency rate (integral of the slip rate over depth, m^2/s)
// and the potency (integral of the slip, m^2). A restart appends, so steps after the last
// checkpoint can appear twice (identical lines); keep the last line of each step.
class FaultSeries
{
private:
  FaultSeries(const FaultSeries& that);
  FaultSeries& operator=(const FaultSeries& rhs);

  FILE                                   *_file;    // open on the first process only
  std::vector<Vec>                        _weights; // per fault: trapezoid weights in z (m)
  std::vector< std::vector<PetscScalar> > _depths;  // per fault: z (km) of every node, on every process

public:
  FaultSeries();
  ~FaultSeries();
  PetscErrorCode setup(const std::vector<Fault_qd*>& faults);
  PetscErrorCode write(const std::string& outputDir, const bool append, const PetscInt stepCount,
    const PetscScalar time, const PetscScalar deltaT, const std::vector<Fault_qd*>& faults);
  PetscErrorCode flush();
};

// Frictional work of the interior faults spread into the body, Q = sum_k MapV(tau_k V_k) .* Gw_k,
// with Gw_k a Gaussian of width w(z) centred on fault k (docs/REVERSIBLE_STRENGTH_PLAN.md,
// section 3.4). Each Gw_k is normalized on the grid, sum_i Wy_i Gw_k(y_i, z) = 1 at every depth with
// Wy the y-quadrature weights of the SBP norm, so the body integral of Q equals the faults' work
// integrated over depth, whatever w is relative to the cell size. Units: tau (MPa) times V (m/s)
// times Gw (1/km) gives kW/m^3, the unit of the heat equation's sources.
class FaultWorkKernel
{
private:
  FaultWorkKernel(const FaultWorkKernel& that);
  FaultWorkKernel& operator=(const FaultWorkKernel& rhs);

  std::vector<InteriorFaultLift*> _lifts; // not owned
  Vec                             _work;  // body work vector
  Vec                             _workFault; // fault-size work vector

public:
  std::vector<Vec> _Gw; // per fault: the normalized kernel (body field), NULL for a boundary fault

  FaultWorkKernel();
  ~FaultWorkKernel();
  // w: the width (km) at every body node, as HeatEquation::_w; sbp: the momentum balance's operators (for the norm)
  PetscErrorCode setup(Domain& D, SbpOps* sbp, const Vec& w, const std::vector<Fault_qd*>& faults,
    const std::vector<InteriorFaultLift*>& lifts);
  // Q = sum over interior faults of MapV(tau_k V_k) .* Gw_k, with tau_k = tauP and V_k = slipVel
  PetscErrorCode spread(const std::vector<Fault_qd*>& faults, Vec& Q);
  // each kernel as dataset Gw in its fault's group of data_context.h5
  PetscErrorCode writeContext(PetscViewer& viewer, const std::vector<Fault_qd*>& faults);
};

#endif

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
// evolving cohesion: explicit cohesion and the reseal ceiling join a non-empty timeIntInds (scaled by
// the largest ceiling and limit); implicit cohesion needs timeIntegrator = RK32_WBE or RK43_WBE
PetscErrorCode prepareCohesion(std::vector<std::string>& inds, std::vector<double>& scale,
  const std::vector<Fault_qd*>& faults, const std::string& timeIntegrator);

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
  // probes: body fields averaged near each fault at one depth, and the fault's cohesion there
  PetscScalar                             _probeDepth;   // the grid depth used (km), < 0 without probes
  std::vector<std::string>                _probeNames, _probeUnits;
  std::vector<const Vec*>                 _probeFields;
  std::vector<Vec>                        _probeWeights; // per fault: 1/count on the averaged body nodes
  std::vector<Vec>                        _probeNode;    // per fault: 1 on the fault node at that depth

public:
  FaultSeries();
  ~FaultSeries();
  PetscErrorCode setup(const std::vector<Fault_qd*>& faults);
  // probes at the grid depth nearest depth (km): body fields averaged over the nodes within width (km)
  // of each fault (always the rows next to it), and an evolving cohesion at that depth
  PetscErrorCode setProbes(Domain& D, const std::vector<Fault_qd*>& faults, const std::vector<InteriorFaultLift*>& lifts,
    const PetscScalar depth, const PetscScalar width);
  // a body field, written as <fault>_<name>@<depth>km<unit> (unit "(GPa*s)", without spaces, or "")
  PetscErrorCode addProbe(const std::string& name, const std::string& unit, const Vec* field);
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
// work = 0 where |V| < vMin (the faults' work done at seismic slip rates only)
PetscErrorCode maskBelowSlipRate(Vec& work, const Vec& V, const PetscScalar vMin);

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
  // Q = sum over interior faults of MapV(tau_k V_k) .* Gw_k, with tau_k = tauP and V_k = slipVel; with
  // vMin > 0 only the work done where |V_k| >= vMin (pseudotachylite: seismic slip)
  PetscErrorCode spread(const std::vector<Fault_qd*>& faults, Vec& Q, const PetscScalar vMin = 0.0);
  // each kernel as dataset Gw in its fault's group of data_context.h5
  PetscErrorCode writeContext(PetscViewer& viewer, const std::vector<Fault_qd*>& faults);
};

// The bottom boundary moving with the plates (momBal_bcB_qd = movingBase): the base displacement is
//   u(y, Lz, t) = (vL t / faultTypeScale) s(y) + shift(y),   s(y) = tanh((y - yc) / w),
// imposed as a Dirichlet condition. s goes from -1 to 1 across a transition centred at yc
// (momBal_bcB_center, km; default mid-domain) of width W (momBal_bcB_width, km) that holds 90% of
// the velocity change, w = W / (2 atanh 0.9); W = 0 is a step at yc. With a boundary fault
// (half-space, y >= 0) the transition is centred on the fault and s goes from 0 to 1; there a step
// is s = 1 everywhere. shift is the base row of the initial displacement (steady-state guess, file
// or zero) and is kept in the checkpoint. Where a fault meets the base the base pins it, so an
// interior fault must be locked at its deepest node, and a boundary fault too unless W = 0, when
// it must not be.
class MovingBase
{
private:
  MovingBase(const MovingBase& that);
  MovingBase& operator=(const MovingBase& rhs);

public:
  PetscScalar _center, _width; // (km); _center is unset (default) while below -1e30
  Vec         _profile, _shift; // s(y) and shift(y) (m), size Ny, laid out like Domain::_z0

  MovingBase();
  ~MovingBase();
  // the profile from the y of the bottom grid row; halfSpace: the left boundary is a fault
  PetscErrorCode setup(Domain& D, const bool halfSpace);
  PetscErrorCode checkFaults(const std::vector<Fault_qd*>& faults, const std::vector<InteriorFaultLift*>& lifts);
  PetscErrorCode setShift(Domain& D, const Vec& u); // shift = u on the bottom row
  PetscErrorCode update(Vec& bcB, const PetscScalar time, const PetscScalar vL, const PetscScalar faultTypeScale);
  PetscErrorCode rate(Vec& bcBRate, const PetscScalar vL, const PetscScalar faultTypeScale);
  PetscErrorCode writeContext(PetscViewer& ascii);   // mediator.txt lines
  PetscErrorCode writeCheckpoint(PetscViewer& viewer);
  PetscErrorCode loadCheckpoint(PetscViewer& viewer);
};

#endif

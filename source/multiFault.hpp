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

#endif

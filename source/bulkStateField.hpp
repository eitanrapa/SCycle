#ifndef BULKSTATEFIELD_HPP_INCLUDED
#define BULKSTATEFIELD_HPP_INCLUDED

#include <petscksp.h>
#include <petscviewerhdf5.h>
#include <string>
#include <vector>
#include <map>
#include "genFuncs.hpp"
#include "domain.hpp"

class PowerLaw;

/*
 * A scalar state field on the body (Ny x Nz) that evolves with the deformation of the bulk and
 * feeds back into the power-law viscosity (docs/REVERSIBLE_STRENGTH_PLAN.md, section 3.1): strain
 * hardening (hard_), water content (water_), ... Each field is an explicit integrand, varEx[<name>],
 * written to the group /<name> of data_2D.h5 and the checkpoint, with its keys <prefix><key> and the
 * context file <name>.txt. A derived class supplies the law (computeRate), its effect on the power
 * law (pushToMaterial) and a time-step bound; this class does the rest.
 *
 * Common keys: <prefix>type = off | transient | constant (held at its initial value); the initial
 * value <prefix><symbol>Vals/Depths or the file <prefix><symbol> in inputDir; and a test hook,
 * <prefix>eTest (1/s, off when negative): the law is driven by the prescribed strain rate
 * eTest (1 + eTestAmp sin(2 pi t/eTestPeriod)) instead of the dislocation strain rate, for the
 * analytic unit tests.
 */

// what the laws may use, at the current rate evaluation
struct BulkInputs
{
  PetscScalar time;
  Vec sdev;        // deviatoric stress (MPa)
  Vec dgVdev;      // viscous strain rate (1e-3/s, engineering shear)
  Vec dgVdev_disl; // its dislocation-creep part (1e-3/s)
  Vec T;           // temperature (K)
  Vec Qfault;      // the faults' work spread into the body (kW/m^3), NULL where the mediator has none
};

class BulkStateField
{
private:
  BulkStateField(const BulkStateField& that);
  BulkStateField& operator=(const BulkStateField& rhs);

public:
  Domain      *_D;
  const char  *_file;
  std::string  _delim, _inputDir, _outputDir;
  std::string  _name, _prefix, _symbol; // integrand key and group; key prefix; dataset name of the state
  std::string  _type;                   // off, transient or constant
  Vec          _state, _rate;           // the field and its rate (body)
  Vec          _e;                      // the driving strain rate of the last evaluation (1/s)
  std::vector<double> _initVals, _initDepths;
  PetscScalar  _eTest, _eTestAmp, _eTestPeriod;

  BulkStateField(Domain& D, const std::string& name, const std::string& prefix, const std::string& symbol);
  virtual ~BulkStateField();

  // the law: _rate from the inputs and the driving strain rate _e
  virtual PetscErrorCode computeRate(const BulkInputs& in) = 0;
  // its effect on the power law, from the current _state
  virtual PetscErrorCode pushToMaterial(PowerLaw& material) = 0;
  // the largest step (s) the law allows; none by default
  virtual PetscErrorCode computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT);
  // the derived class's lines of <name>.txt and fields of data_context.h5
  virtual PetscErrorCode writeContextExtra(PetscViewer& ascii, PetscViewer& viewer);
  // keep the state within its range after each update (none by default)
  virtual PetscErrorCode clamp();
  // the initial state when no <prefix><symbol>Vals are given (0 by default)
  virtual PetscErrorCode initialDefault(Vec& state);
  // the law is driven by the total viscous strain rate rather than its dislocation-creep part
  virtual bool usesTotalViscousRate() const { return false; }
  // the law needs BulkInputs::Qfault (the mediator then spreads the faults' work before its rate)
  virtual bool needsFaultWork() const { return false; }

  // life cycle, called by the mediator
  PetscErrorCode initiateIntegrand(const PetscScalar time, std::map<std::string,Vec>& varEx);
  PetscErrorCode updateFields(const PetscScalar time, const std::map<std::string,Vec>& varEx);
  PetscErrorCode d_dt(const BulkInputs& in, std::map<std::string,Vec>& dvarEx);
  PetscErrorCode addErrorControl(std::vector<std::string>& errInds, std::vector<double>& scale) const;
  PetscErrorCode writeContext(PetscViewer& viewer);
  PetscErrorCode writeStep(PetscViewer& viewer);
  PetscErrorCode writeCheckpoint(PetscViewer& viewer);
  PetscErrorCode loadCheckpoint();

protected:
  // the common keys; the derived class's parser calls this first
  bool parseCommon(const std::string& var, const std::string& rhs, const std::string& rhsFull);
  // reads the input file into (var, rhs, rhsFull) lines
  void readLines(std::vector<std::string>& vars, std::vector<std::string>& rhss, std::vector<std::string>& rhsFulls) const;
  // allocate the fields and set the initial state (keys, file or checkpoint); after the settings
  PetscErrorCode setup();
  // _e = dgVdev_disl * 1e-3 (1/s), or the test rate
  PetscErrorCode drivingStrainRate(const BulkInputs& in);
};

// Stops with a message if the input file turns on a bulk state field (hard_type, water_type, ...)
// in a mediator that does not evolve them.
PetscErrorCode refuseBulkStates(const char* file, const std::string& delim, const std::string& mediator);

#endif

#include "segregationState.hpp"
#include "powerLaw.hpp"

#define FILENAME "segregationState.cpp"

using namespace std;


SegregationState::SegregationState(Domain& D)
: BulkStateField(D,"segregation","seg_","Xi"),
  _betaS(1.0), _tauM0(-1), _QM(0), _TRef(273.15), _gammaM(0), _eMixTest(-1), _gammaS(NULL), _f(NULL)
{
  vector<string> vars, rhss, rhsFulls;
  readLines(vars,rhss,rhsFulls);
  for (size_t i = 0; i < vars.size(); i++) {
    const string& var = vars[i]; const string& rhs = rhss[i];
    if (parseCommon(var,rhs,rhsFulls[i])) { continue; }
    else if (var == "seg_betaS") { _betaS = atof(rhs.c_str()); }
    else if (var == "seg_gammaSVals") { _gammaSVals.clear(); loadVectorFromInputFile(rhsFulls[i],_gammaSVals); }
    else if (var == "seg_gammaSDepths") { _gammaSDepths.clear(); loadVectorFromInputFile(rhsFulls[i],_gammaSDepths); }
    else if (var == "seg_tauM0") { _tauM0 = atof(rhs.c_str()); }
    else if (var == "seg_QM") { _QM = atof(rhs.c_str()); }
    else if (var == "seg_TRef") { _TRef = atof(rhs.c_str()); }
    else if (var == "seg_gammaM") { _gammaM = atof(rhs.c_str()); }
    else if (var == "seg_eMixTest") { _eMixTest = atof(rhs.c_str()); }
  }
  if (_type == "off") { return; }

  if (_gammaSVals.empty() || _gammaSVals.size() != _gammaSDepths.size()) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: phase segregation needs seg_gammaSVals and seg_gammaSDepths of the same length.\n");
    assert(0);
  }
  if (!(_tauM0 > 0) || !(_TRef > 0) || !(_betaS >= 0) || !(_gammaM >= 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: phase segregation needs seg_tauM0 > 0 (s), seg_TRef > 0 (K), seg_betaS >= 0 and seg_gammaM >= 0.\n");
    assert(0);
  }
  if (_eMixTest >= 0 && !(_gammaM > 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: seg_eMixTest drives the remixing by dislocation creep: it needs seg_gammaM > 0.\n");
    assert(0);
  }
  PetscErrorCode ierr;
  VecDuplicate(D._y,&_gammaS);
  ierr = setVec(_gammaS,D._z,_gammaSVals,_gammaSDepths); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  PetscScalar gmin = 0;
  VecMin(_gammaS,NULL,&gmin);
  if (!(gmin > 0)) { PetscPrintf(PETSC_COMM_WORLD,"Error: seg_gammaSVals must be positive.\n"); assert(0); }
  PetscObjectSetName((PetscObject) _gammaS, "gammaS");
  VecDuplicate(D._y,&_f);
  PetscObjectSetName((PetscObject) _f, "f");

  ierr = setup(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}


SegregationState::~SegregationState()
{
  VecDestroy(&_gammaS);
  VecDestroy(&_f);
}


// dXi/dt = (e_dp/gamma_s)(1 - Xi) - Xi (1/tau_m(T) + e_disl/gamma_m)
PetscErrorCode SegregationState::computeRate(const BulkInputs& in)
{
  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *x, *e, *g, *T, *ed;
  PetscScalar *r;
  VecGetArrayRead(_state,&x); VecGetArrayRead(_e,&e); VecGetArrayRead(_gammaS,&g); VecGetArrayRead(in.T,&T);
  VecGetArrayRead(in.dgVdev_disl,&ed);
  VecGetArray(_rate,&r);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    r[Jj] = (fabs(e[Jj])/g[Jj])*(1.0 - x[Jj]) - x[Jj]*(1.0/tauM(T[Jj]) + mixRate(fabs(ed[Jj])*1e-3));
  }
  VecRestoreArrayRead(_state,&x); VecRestoreArrayRead(_e,&e); VecRestoreArrayRead(_gammaS,&g); VecRestoreArrayRead(in.T,&T);
  VecRestoreArrayRead(in.dgVdev_disl,&ed);
  VecRestoreArray(_rate,&r);
  return ierr;
}


// pressure solution's strain rate at a given stress times 1/(1 + beta_s Xi)
PetscErrorCode SegregationState::pushToMaterial(PowerLaw& material)
{
  PetscErrorCode ierr = 0;
  ierr = VecCopy(_state,_f); CHKERRQ(ierr);
  ierr = VecScale(_f,_betaS); CHKERRQ(ierr);
  ierr = VecShift(_f,1.0); CHKERRQ(ierr);
  ierr = VecReciprocal(_f); CHKERRQ(ierr);
  ierr = material.setDPRateFactor(_name,_f); CHKERRQ(ierr);
  return ierr;
}


// the relaxation time of the law, 1/(e_dp/gamma_s + 1/tau_m + e_disl/gamma_m), at its shortest over the nodes
PetscErrorCode SegregationState::computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT)
{
  PetscErrorCode ierr = 0;
  maxDeltaT = PETSC_MAX_REAL;
  if (_type != "transient") { return ierr; }
  ierr = drivingStrainRate(in); CHKERRQ(ierr);
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_e,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *e, *g, *T, *ed;
  VecGetArrayRead(_e,&e); VecGetArrayRead(_gammaS,&g); VecGetArrayRead(in.T,&T); VecGetArrayRead(in.dgVdev_disl,&ed);
  PetscScalar dt = PETSC_MAX_REAL;
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    dt = PetscMin(dt, 1.0/(fabs(e[Jj])/g[Jj] + 1.0/tauM(T[Jj]) + mixRate(fabs(ed[Jj])*1e-3)));
  }
  VecRestoreArrayRead(_e,&e); VecRestoreArrayRead(_gammaS,&g); VecRestoreArrayRead(in.T,&T); VecRestoreArrayRead(in.dgVdev_disl,&ed);
  ierr = MPIU_Allreduce(&dt,&maxDeltaT,1,MPIU_SCALAR,MPIU_MIN,PETSC_COMM_WORLD); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode SegregationState::clamp()
{
  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  PetscScalar *x;
  VecGetArray(_state,&x);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) { x[Jj] = PetscMin(1.0,PetscMax(0.0,x[Jj])); }
  VecRestoreArray(_state,&x);
  return ierr;
}


PetscErrorCode SegregationState::writeContextExtra(PetscViewer& ascii, PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerASCIIPrintf(ascii,"seg_betaS = %.15e\n",_betaS); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"seg_tauM0 = %.15e # (s)\n",_tauM0); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"seg_QM = %.15e # (K)\n",_QM); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"seg_TRef = %.15e # (K)\n",_TRef); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"seg_gammaM = %.15e # (0: no remixing by dislocation creep)\n",_gammaM); CHKERRQ(ierr);
  if (_eMixTest >= 0) { ierr = PetscViewerASCIIPrintf(ascii,"seg_eMixTest = %.15e # (1/s) test: prescribed e_disl\n",_eMixTest); CHKERRQ(ierr); }
  ierr = PetscViewerHDF5PushGroup(viewer,"/segregation"); CHKERRQ(ierr);
  ierr = VecView(_gammaS,viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  return ierr;
}

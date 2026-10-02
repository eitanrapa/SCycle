#include "waterState.hpp"
#include "powerLaw.hpp"

#define FILENAME "waterState.cpp"

using namespace std;


WaterState::WaterState(Domain& D, const Vec& wetDist0)
: BulkStateField(D,"water","water_","chi"),
  _chiSup(NULL), _chiMin(0.3), _tauHyd(-1), _kappa(0), _gammaDry(-1), _chiC(0),
  _mix("log"), _supplyFrom("profile"), _wetDist0(NULL)
{
  vector<string> vars, rhss, rhsFulls;
  readLines(vars,rhss,rhsFulls);
  for (size_t i = 0; i < vars.size(); i++) {
    const string& var = vars[i]; const string& rhs = rhss[i];
    if (parseCommon(var,rhs,rhsFulls[i])) { continue; }
    else if (var == "water_chiSupVals") { _chiSupVals.clear(); loadVectorFromInputFile(rhsFulls[i],_chiSupVals); }
    else if (var == "water_chiSupDepths") { _chiSupDepths.clear(); loadVectorFromInputFile(rhsFulls[i],_chiSupDepths); }
    else if (var == "water_chiMin") { _chiMin = atof(rhs.c_str()); }
    else if (var == "water_tauHyd") { _tauHyd = atof(rhs.c_str()); }
    else if (var == "water_kappa") { _kappa = atof(rhs.c_str()); }
    else if (var == "water_gammaDry") { _gammaDry = atof(rhs.c_str()); }
    else if (var == "water_chiC") { _chiC = atof(rhs.c_str()); }
    else if (var == "water_mix") { _mix = rhs; }
    else if (var == "water_supplyFrom") { _supplyFrom = rhs; }
  }
  if (_type == "off") { return; }

  if (!(_tauHyd > 0) || !(_gammaDry > 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: water content needs water_tauHyd > 0 (s) and water_gammaDry > 0.\n");
    assert(0);
  }
  if (!(_chiMin >= 0 && _chiMin <= 1) || !(_chiC >= 0 && _chiC < 1) || !(_kappa >= 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: water_chiMin must be in [0, 1], water_chiC in [0, 1), water_kappa >= 0.\n");
    assert(0);
  }
  if (_mix != "log" && _mix != "arithmetic") {
    PetscPrintf(PETSC_COMM_WORLD,"Error: water_mix must be log or arithmetic.\n");
    assert(0);
  }
  if (_supplyFrom != "profile") {
    PetscPrintf(PETSC_COMM_WORLD,"Error: water_supplyFrom = %s is not implemented yet; use profile (water_chiSupVals).\n",_supplyFrom.c_str());
    assert(0);
  }
  if (_chiSupVals.size() != _chiSupDepths.size()) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: water_chiSupVals and water_chiSupDepths must have the same length.\n");
    assert(0);
  }

  PetscErrorCode ierr;
  VecDuplicate(D._y,&_chiSup);
  PetscObjectSetName((PetscObject) _chiSup, "chiSup");
  if (_chiSupVals.empty()) { VecSet(_chiSup,1.0); }
  else { ierr = setVec(_chiSup,D._z,_chiSupVals,_chiSupDepths); CHKERRABORT(PETSC_COMM_WORLD,ierr); }
  VecDuplicate(wetDist0,&_wetDist0);
  VecCopy(wetDist0,_wetDist0);

  ierr = setup(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}


WaterState::~WaterState()
{
  VecDestroy(&_chiSup);
  VecDestroy(&_wetDist0);
}


PetscErrorCode WaterState::initialDefault(Vec& state) { return VecCopy(_wetDist0,state); }


// dchi/dt = (chi_sup - chi)/tau_hyd (1 + kappa (1 - chi)) - (e/gamma_dry)(chi - chi_min)
PetscErrorCode WaterState::computeRate(const BulkInputs& in)
{
  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *c, *e, *sup;
  PetscScalar *r;
  VecGetArrayRead(_state,&c); VecGetArrayRead(_e,&e); VecGetArrayRead(_chiSup,&sup);
  VecGetArray(_rate,&r);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    r[Jj] = (sup[Jj] - c[Jj])/_tauHyd*(1.0 + _kappa*(1.0 - c[Jj])) - (e[Jj]/_gammaDry)*(c[Jj] - _chiMin);
  }
  VecRestoreArrayRead(_state,&c); VecRestoreArrayRead(_e,&e); VecRestoreArrayRead(_chiSup,&sup);
  VecRestoreArray(_rate,&r);
  return ierr;
}


PetscErrorCode WaterState::pushToMaterial(PowerLaw& material)
{
  return material.updateWetDist(_state,_chiC,_mix);
}


// the relaxation time of the law, 1/((1 + kappa)/tau_hyd + e/gamma_dry), at its shortest
PetscErrorCode WaterState::computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT)
{
  PetscErrorCode ierr = 0;
  maxDeltaT = PETSC_MAX_REAL;
  if (_type != "transient") { return ierr; }
  ierr = drivingStrainRate(in); CHKERRQ(ierr);
  PetscScalar emax = 0;
  ierr = VecNorm(_e,NORM_INFINITY,&emax); CHKERRQ(ierr);
  maxDeltaT = 1.0/((1.0 + _kappa)/_tauHyd + emax/_gammaDry);
  return ierr;
}


PetscErrorCode WaterState::clamp()
{
  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  PetscScalar *c;
  VecGetArray(_state,&c);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) { c[Jj] = PetscMin(1.0,PetscMax(_chiMin,c[Jj])); }
  VecRestoreArray(_state,&c);
  return ierr;
}


PetscErrorCode WaterState::writeContextExtra(PetscViewer& ascii, PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerASCIIPrintf(ascii,"water_chiMin = %.15e\n",_chiMin); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"water_tauHyd = %.15e # (s)\n",_tauHyd); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"water_kappa = %.15e\n",_kappa); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"water_gammaDry = %.15e\n",_gammaDry); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"water_chiC = %.15e\n",_chiC); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"water_mix = %s\n",_mix.c_str()); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"water_supplyFrom = %s\n",_supplyFrom.c_str()); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PushGroup(viewer,"/water"); CHKERRQ(ierr);
  ierr = VecView(_chiSup,viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  return ierr;
}

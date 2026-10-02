#include "cementState.hpp"
#include "powerLaw.hpp"

#define FILENAME "cementState.cpp"

using namespace std;


CementState::CementState(Domain& D)
: BulkStateField(D,"cement","cement_","Psi"),
  _betaC(0.5), _Ec(-1), _tauPsi(-1), _gammaPsi(0), _QTest(-1), _F(NULL)
{
  vector<string> vars, rhss, rhsFulls;
  readLines(vars,rhss,rhsFulls);
  for (size_t i = 0; i < vars.size(); i++) {
    const string& var = vars[i]; const string& rhs = rhss[i];
    if (parseCommon(var,rhs,rhsFulls[i])) { continue; }
    else if (var == "cement_betaC") { _betaC = atof(rhs.c_str()); }
    else if (var == "cement_Ec") { _Ec = atof(rhs.c_str()); }
    else if (var == "cement_tauPsi") { _tauPsi = atof(rhs.c_str()); }
    else if (var == "cement_gammaPsi") { _gammaPsi = atof(rhs.c_str()); }
    else if (var == "cement_QTest") { _QTest = atof(rhs.c_str()); }
  }
  if (_type == "off") { return; }

  if (!(_Ec > 0) || !(_tauPsi > 0) || !(_betaC >= 0) || !(_gammaPsi >= 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: cement needs cement_Ec > 0 (kJ/m^3), cement_tauPsi > 0 (s), cement_betaC >= 0 and\n"
      "       cement_gammaPsi >= 0 (0: no removal by deformation).\n");
    assert(0);
  }
  VecDuplicate(D._y,&_F);
  PetscObjectSetName((PetscObject) _F, "F");

  PetscErrorCode ierr = setup(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}


CementState::~CementState()
{
  VecDestroy(&_F);
}


// dPsi/dt = (Q/E_c)(1 - Psi) - Psi (1/tau_Psi + e/gamma_Psi)
PetscErrorCode CementState::computeRate(const BulkInputs& in)
{
  PetscErrorCode ierr = 0;
  if (_QTest < 0 && in.Qfault == NULL) {
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"the cement state needs the faults' work (Qfault)");
  }
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *p, *e, *Q = NULL;
  PetscScalar *r;
  VecGetArrayRead(_state,&p); VecGetArrayRead(_e,&e);
  if (_QTest < 0) { VecGetArrayRead(in.Qfault,&Q); }
  VecGetArray(_rate,&r);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    const PetscScalar q = (Q != NULL) ? PetscMax(0.0,Q[Jj]) : _QTest;
    const PetscScalar removal = 1.0/_tauPsi + ((_gammaPsi > 0) ? fabs(e[Jj])/_gammaPsi : 0.0);
    r[Jj] = (q/_Ec)*(1.0 - p[Jj]) - p[Jj]*removal;
  }
  VecRestoreArrayRead(_state,&p); VecRestoreArrayRead(_e,&e);
  if (Q != NULL) { VecRestoreArrayRead(in.Qfault,&Q); }
  VecRestoreArray(_rate,&r);
  return ierr;
}


// strength factor 1 + beta_c Psi on every creep mechanism
PetscErrorCode CementState::pushToMaterial(PowerLaw& material)
{
  PetscErrorCode ierr = 0;
  ierr = VecCopy(_state,_F); CHKERRQ(ierr);
  ierr = VecScale(_F,_betaC); CHKERRQ(ierr);
  ierr = VecShift(_F,1.0); CHKERRQ(ierr);
  ierr = material.setStrengthFactor(_name,_F); CHKERRQ(ierr);
  return ierr;
}


// the relaxation time of the law, 1/(Q/E_c + 1/tau_Psi + e/gamma_Psi), at its shortest
PetscErrorCode CementState::computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT)
{
  PetscErrorCode ierr = 0;
  maxDeltaT = PETSC_MAX_REAL;
  if (_type != "transient") { return ierr; }
  ierr = drivingStrainRate(in); CHKERRQ(ierr);
  PetscScalar qmax = _QTest, emax = 0;
  if (_QTest < 0) {
    if (in.Qfault == NULL) { return ierr; }
    ierr = VecMax(in.Qfault,NULL,&qmax); CHKERRQ(ierr);
  }
  ierr = VecNorm(_e,NORM_INFINITY,&emax); CHKERRQ(ierr);
  maxDeltaT = 1.0/(PetscMax(0.0,qmax)/_Ec + 1.0/_tauPsi + ((_gammaPsi > 0) ? emax/_gammaPsi : 0.0));
  return ierr;
}


PetscErrorCode CementState::clamp()
{
  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  PetscScalar *p;
  VecGetArray(_state,&p);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) { p[Jj] = PetscMin(1.0,PetscMax(0.0,p[Jj])); }
  VecRestoreArray(_state,&p);
  return ierr;
}


PetscErrorCode CementState::writeContextExtra(PetscViewer& ascii, PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerASCIIPrintf(ascii,"cement_betaC = %.15e\n",_betaC); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"cement_Ec = %.15e # (kJ/m^3)\n",_Ec); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"cement_tauPsi = %.15e # (s)\n",_tauPsi); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"cement_gammaPsi = %.15e # (0: no removal by deformation)\n",_gammaPsi); CHKERRQ(ierr);
  if (_QTest >= 0) { ierr = PetscViewerASCIIPrintf(ascii,"cement_QTest = %.15e # (kW/m^3) test fault work\n",_QTest); CHKERRQ(ierr); }
  return ierr;
}

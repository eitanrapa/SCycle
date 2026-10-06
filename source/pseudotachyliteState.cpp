#include "pseudotachyliteState.hpp"
#include "powerLaw.hpp"

#define FILENAME "pseudotachyliteState.cpp"

using namespace std;


PseudotachyliteState::PseudotachyliteState(Domain& D)
: BulkStateField(D,"pseudotachylite","pt_","phiPT"),
  _eMelt(3.6e6), _vMelt(0.1), _tauA0(-1), _QA(0), _TRef(273.15), _QTest(-1), _dMelt(NULL)
{
  vector<string> vars, rhss, rhsFulls;
  readLines(vars,rhss,rhsFulls);
  for (size_t i = 0; i < vars.size(); i++) {
    const string& var = vars[i]; const string& rhs = rhss[i];
    if (parseCommon(var,rhs,rhsFulls[i])) { continue; }
    else if (var == "pt_dMeltVals") { _dMeltVals.clear(); loadVectorFromInputFile(rhsFulls[i],_dMeltVals); }
    else if (var == "pt_dMeltDepths") { _dMeltDepths.clear(); loadVectorFromInputFile(rhsFulls[i],_dMeltDepths); }
    else if (var == "pt_eMelt") { _eMelt = atof(rhs.c_str()); }
    else if (var == "pt_vMelt") { _vMelt = atof(rhs.c_str()); }
    else if (var == "pt_tauA0") { _tauA0 = atof(rhs.c_str()); }
    else if (var == "pt_QA") { _QA = atof(rhs.c_str()); }
    else if (var == "pt_TRef") { _TRef = atof(rhs.c_str()); }
    else if (var == "pt_QTest") { _QTest = atof(rhs.c_str()); }
  }
  if (_type == "off") { return; }

  if (_dMeltVals.empty() || _dMeltVals.size() != _dMeltDepths.size()) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: pseudotachylite products need pt_dMeltVals and pt_dMeltDepths of the same length.\n");
    assert(0);
  }
  if (!(_eMelt > 0) || !(_vMelt >= 0) || !(_tauA0 > 0) || !(_TRef > 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: pseudotachylite products need pt_eMelt > 0 (kJ/m^3), pt_vMelt >= 0 (m/s), pt_tauA0 > 0 (s) and pt_TRef > 0 (K).\n");
    assert(0);
  }
  PetscErrorCode ierr;
  VecDuplicate(D._y,&_dMelt);
  ierr = setVec(_dMelt,D._z,_dMeltVals,_dMeltDepths); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  PetscScalar dmin = 0;
  VecMin(_dMelt,NULL,&dmin);
  if (!(dmin > 0)) { PetscPrintf(PETSC_COMM_WORLD,"Error: pt_dMeltVals must be positive.\n"); assert(0); }
  PetscObjectSetName((PetscObject) _dMelt, "dMelt");

  ierr = setup(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}


PseudotachyliteState::~PseudotachyliteState()
{
  VecDestroy(&_dMelt);
}


// dphi/dt = (Q_melt/e_melt)(1 - phi) - phi/tau_a(T)
PetscErrorCode PseudotachyliteState::computeRate(const BulkInputs& in)
{
  PetscErrorCode ierr = 0;
  if (_QTest < 0 && in.Qmelt == NULL) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"pseudotachylite products need the faults' seismic work (Qmelt)"); }
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *p, *T, *Q = NULL;
  PetscScalar *r;
  VecGetArrayRead(_state,&p); VecGetArrayRead(in.T,&T);
  if (_QTest < 0) { VecGetArrayRead(in.Qmelt,&Q); }
  VecGetArray(_rate,&r);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    const PetscScalar q = PetscMax(0.0, (Q != NULL) ? Q[Jj] : _QTest);
    r[Jj] = (q/_eMelt)*(1.0 - p[Jj]) - p[Jj]/tauA(T[Jj]);
  }
  VecRestoreArrayRead(_state,&p); VecRestoreArrayRead(in.T,&T);
  if (Q != NULL) { VecRestoreArrayRead(in.Qmelt,&Q); }
  VecRestoreArray(_rate,&r);
  return ierr;
}


// grain-size-sensitive creep as the mixture of host grains and products
PetscErrorCode PseudotachyliteState::pushToMaterial(PowerLaw& material)
{
  return material.setMeltProducts(_state,_dMelt);
}


// the relaxation time of the law, 1/(Q_melt/e_melt + 1/tau_a), at its shortest over the nodes
PetscErrorCode PseudotachyliteState::computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT)
{
  PetscErrorCode ierr = 0;
  maxDeltaT = PETSC_MAX_REAL;
  if (_type != "transient") { return ierr; }
  if (_QTest < 0 && in.Qmelt == NULL) { return ierr; }
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *T, *Q = NULL;
  VecGetArrayRead(in.T,&T);
  if (_QTest < 0) { VecGetArrayRead(in.Qmelt,&Q); }
  PetscScalar dt = PETSC_MAX_REAL;
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    const PetscScalar q = PetscMax(0.0, (Q != NULL) ? Q[Jj] : _QTest);
    dt = PetscMin(dt, 1.0/(q/_eMelt + 1.0/tauA(T[Jj])));
  }
  VecRestoreArrayRead(in.T,&T);
  if (Q != NULL) { VecRestoreArrayRead(in.Qmelt,&Q); }
  ierr = MPIU_Allreduce(&dt,&maxDeltaT,1,MPIU_SCALAR,MPIU_MIN,PETSC_COMM_WORLD); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode PseudotachyliteState::clamp()
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


PetscErrorCode PseudotachyliteState::writeContextExtra(PetscViewer& ascii, PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerASCIIPrintf(ascii,"pt_eMelt = %.15e # (kJ/m^3)\n",_eMelt); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"pt_vMelt = %.15e # (m/s)\n",_vMelt); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"pt_tauA0 = %.15e # (s)\n",_tauA0); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"pt_QA = %.15e # (K)\n",_QA); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"pt_TRef = %.15e # (K)\n",_TRef); CHKERRQ(ierr);
  if (_QTest >= 0) { ierr = PetscViewerASCIIPrintf(ascii,"pt_QTest = %.15e # (kW/m^3) test: uniform melting power\n",_QTest); CHKERRQ(ierr); }
  ierr = PetscViewerHDF5PushGroup(viewer,"/pseudotachylite"); CHKERRQ(ierr);
  ierr = VecView(_dMelt,viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  return ierr;
}

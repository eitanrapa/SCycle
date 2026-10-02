#include "hardeningState.hpp"
#include "powerLaw.hpp"

#define FILENAME "hardeningState.cpp"

using namespace std;


HardeningState::HardeningState(Domain& D, const Vec& T0)
: BulkStateField(D,"hard","hard_","S"),
  _a(0.5), _tauR0(-1), _QR(0), _TRef(273.15), _SRef(-1), _eRef(-1),
  _gammaH(NULL), _SRefVec(NULL), _H(NULL)
{
  vector<string> vars, rhss, rhsFulls;
  readLines(vars,rhss,rhsFulls);
  for (size_t i = 0; i < vars.size(); i++) {
    const string& var = vars[i]; const string& rhs = rhss[i];
    if (parseCommon(var,rhs,rhsFulls[i])) { continue; }
    else if (var == "hard_a") { _a = atof(rhs.c_str()); }
    else if (var == "hard_gammaHVals") { _gammaHVals.clear(); loadVectorFromInputFile(rhsFulls[i],_gammaHVals); }
    else if (var == "hard_gammaHDepths") { _gammaHDepths.clear(); loadVectorFromInputFile(rhsFulls[i],_gammaHDepths); }
    else if (var == "hard_tauR0") { _tauR0 = atof(rhs.c_str()); }
    else if (var == "hard_QR") { _QR = atof(rhs.c_str()); }
    else if (var == "hard_TRef") { _TRef = atof(rhs.c_str()); }
    else if (var == "hard_SRef") { _SRef = atof(rhs.c_str()); }
    else if (var == "hard_eRef") { _eRef = atof(rhs.c_str()); }
  }
  if (_type == "off") { return; }

  if (_gammaHVals.empty() || _gammaHVals.size() != _gammaHDepths.size()) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: hardening needs hard_gammaHVals and hard_gammaHDepths of the same length.\n");
    assert(0);
  }
  if (!(_tauR0 > 0) || !(_TRef > 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: hardening needs hard_tauR0 > 0 (s) and hard_TRef > 0 (K).\n");
    assert(0);
  }
  if ((_SRef < 0) == (_eRef < 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: give exactly one of hard_SRef (in [0, 1]) and hard_eRef (1/s).\n");
    assert(0);
  }

  PetscErrorCode ierr;
  VecDuplicate(D._y,&_gammaH);
  ierr = setVec(_gammaH,D._z,_gammaHVals,_gammaHDepths); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  PetscScalar gmin = 0;
  VecMin(_gammaH,NULL,&gmin);
  if (!(gmin > 0)) { PetscPrintf(PETSC_COMM_WORLD,"Error: hard_gammaHVals must be positive.\n"); assert(0); }
  PetscObjectSetName((PetscObject) _gammaH, "gammaH");

  // S_ref: given, or the steady state at e_ref with the initial temperature
  VecDuplicate(D._y,&_SRefVec);
  PetscObjectSetName((PetscObject) _SRefVec, "SRef");
  if (_SRef >= 0) { VecSet(_SRefVec,_SRef); }
  else {
    PetscInt Istart, Iend;
    VecGetOwnershipRange(_SRefVec,&Istart,&Iend);
    PetscScalar *s;
    const PetscScalar *T, *g;
    VecGetArray(_SRefVec,&s); VecGetArrayRead(T0,&T); VecGetArrayRead(_gammaH,&g);
    for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
      const PetscScalar x = _eRef*tauR(T[Jj])/g[Jj];
      s[Jj] = x/(1.0 + x);
    }
    VecRestoreArray(_SRefVec,&s); VecRestoreArrayRead(T0,&T); VecRestoreArrayRead(_gammaH,&g);
  }
  VecDuplicate(D._y,&_H);
  PetscObjectSetName((PetscObject) _H, "H");

  // initial state: S_ref unless given (setup calls initialDefault)
  ierr = setup(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}


HardeningState::~HardeningState()
{
  VecDestroy(&_gammaH);
  VecDestroy(&_SRefVec);
  VecDestroy(&_H);
}


// dS/dt = (e/gamma_h) (1 - S) - S/tau_r(T)
PetscErrorCode HardeningState::computeRate(const BulkInputs& in)
{
  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *S, *e, *T, *g;
  PetscScalar *r;
  VecGetArrayRead(_state,&S); VecGetArrayRead(_e,&e); VecGetArrayRead(in.T,&T); VecGetArrayRead(_gammaH,&g);
  VecGetArray(_rate,&r);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    r[Jj] = (e[Jj]/g[Jj])*(1.0 - S[Jj]) - S[Jj]/tauR(T[Jj]);
  }
  VecRestoreArrayRead(_state,&S); VecRestoreArrayRead(_e,&e); VecRestoreArrayRead(in.T,&T); VecRestoreArrayRead(_gammaH,&g);
  VecRestoreArray(_rate,&r);
  return ierr;
}


// H = 1 + a (S - S_ref), to the power law
PetscErrorCode HardeningState::pushToMaterial(PowerLaw& material)
{
  PetscErrorCode ierr = 0;
  ierr = VecWAXPY(_H,-1.0,_SRefVec,_state); CHKERRQ(ierr);
  ierr = VecScale(_H,_a); CHKERRQ(ierr);
  ierr = VecShift(_H,1.0); CHKERRQ(ierr);
  ierr = material.updateHardening(_H); CHKERRQ(ierr);
  return ierr;
}


// the relaxation time of the law, 1/(e/gamma_h + 1/tau_r), at its shortest
PetscErrorCode HardeningState::computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT)
{
  PetscErrorCode ierr = 0;
  maxDeltaT = PETSC_MAX_REAL;
  if (_type != "transient") { return ierr; }
  ierr = drivingStrainRate(in); CHKERRQ(ierr);
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *e, *T, *g;
  VecGetArrayRead(_e,&e); VecGetArrayRead(in.T,&T); VecGetArrayRead(_gammaH,&g);
  PetscScalar local = PETSC_MAX_REAL;
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    local = PetscMin(local,1.0/(PetscAbsScalar(e[Jj])/g[Jj] + 1.0/tauR(T[Jj])));
  }
  VecRestoreArrayRead(_e,&e); VecRestoreArrayRead(in.T,&T); VecRestoreArrayRead(_gammaH,&g);
  ierr = MPI_Allreduce(&local,&maxDeltaT,1,MPIU_SCALAR,MPI_MIN,PETSC_COMM_WORLD); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode HardeningState::initialDefault(Vec& state) { return VecCopy(_SRefVec,state); }


PetscErrorCode HardeningState::clamp()
{
  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  PetscScalar *S;
  VecGetArray(_state,&S);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) { S[Jj] = PetscMin(1.0,PetscMax(0.0,S[Jj])); }
  VecRestoreArray(_state,&S);
  return ierr;
}


PetscErrorCode HardeningState::writeContextExtra(PetscViewer& ascii, PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerASCIIPrintf(ascii,"hard_a = %.15e\n",_a); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"hard_tauR0 = %.15e # (s)\n",_tauR0); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"hard_QR = %.15e # (K)\n",_QR); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"hard_TRef = %.15e # (K)\n",_TRef); CHKERRQ(ierr);
  if (_SRef >= 0) { ierr = PetscViewerASCIIPrintf(ascii,"hard_SRef = %.15e\n",_SRef); CHKERRQ(ierr); }
  else { ierr = PetscViewerASCIIPrintf(ascii,"hard_eRef = %.15e # (1/s)\n",_eRef); CHKERRQ(ierr); }
  ierr = PetscViewerHDF5PushGroup(viewer,"/hard"); CHKERRQ(ierr);
  ierr = VecView(_gammaH,viewer); CHKERRQ(ierr);
  ierr = VecView(_SRefVec,viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  return ierr;
}

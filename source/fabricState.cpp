#include "fabricState.hpp"
#include "powerLaw.hpp"

#define FILENAME "fabricState.cpp"

using namespace std;


FabricState::FabricState(Domain& D)
: BulkStateField(D,"fabric","fabric_","Phi"),
  _betaF(4.0), _tauC0(-1), _QC(0), _TRef(273.15), _eCrit(0), _anisotropic(0), _gammaF(NULL), _F(NULL)
{
  vector<string> vars, rhss, rhsFulls;
  readLines(vars,rhss,rhsFulls);
  for (size_t i = 0; i < vars.size(); i++) {
    const string& var = vars[i]; const string& rhs = rhss[i];
    if (parseCommon(var,rhs,rhsFulls[i])) { continue; }
    else if (var == "fabric_betaF") { _betaF = atof(rhs.c_str()); }
    else if (var == "fabric_gammaFVals") { _gammaFVals.clear(); loadVectorFromInputFile(rhsFulls[i],_gammaFVals); }
    else if (var == "fabric_gammaFDepths") { _gammaFDepths.clear(); loadVectorFromInputFile(rhsFulls[i],_gammaFDepths); }
    else if (var == "fabric_tauC0") { _tauC0 = atof(rhs.c_str()); }
    else if (var == "fabric_QC") { _QC = atof(rhs.c_str()); }
    else if (var == "fabric_TRef") { _TRef = atof(rhs.c_str()); }
    else if (var == "fabric_eCrit") { _eCrit = atof(rhs.c_str()); }
    else if (var == "fabric_anisotropic") { _anisotropic = atoi(rhs.c_str()); }
  }
  if (_type == "off") { return; }

  if (_gammaFVals.empty() || _gammaFVals.size() != _gammaFDepths.size()) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: fabric needs fabric_gammaFVals and fabric_gammaFDepths of the same length.\n");
    assert(0);
  }
  if (!(_tauC0 > 0) || !(_TRef > 0) || !(_betaF >= 0) || !(_eCrit >= 0)) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: fabric needs fabric_tauC0 > 0 (s), fabric_TRef > 0 (K), fabric_betaF >= 0 and fabric_eCrit >= 0.\n");
    assert(0);
  }
  if (_anisotropic != 0 && _anisotropic != 1) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: fabric_anisotropic must be 0 (every direction) or 1 (fault-parallel shear only).\n");
    assert(0);
  }
  PetscErrorCode ierr;
  VecDuplicate(D._y,&_gammaF);
  ierr = setVec(_gammaF,D._z,_gammaFVals,_gammaFDepths); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  PetscScalar gmin = 0;
  VecMin(_gammaF,NULL,&gmin);
  if (!(gmin > 0)) { PetscPrintf(PETSC_COMM_WORLD,"Error: fabric_gammaFVals must be positive.\n"); assert(0); }
  PetscObjectSetName((PetscObject) _gammaF, "gammaF");
  VecDuplicate(D._y,&_F);
  PetscObjectSetName((PetscObject) _F, "F");

  ierr = setup(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}


FabricState::~FabricState()
{
  VecDestroy(&_gammaF);
  VecDestroy(&_F);
}


// dPhi/dt = (e/gamma_f)(1 - Phi) [e < e_crit] - Phi/tau_c(T)
PetscErrorCode FabricState::computeRate(const BulkInputs& in)
{
  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_state,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *p, *e, *g, *T;
  PetscScalar *r;
  VecGetArrayRead(_state,&p); VecGetArrayRead(_e,&e); VecGetArrayRead(_gammaF,&g); VecGetArrayRead(in.T,&T);
  VecGetArray(_rate,&r);
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
    const PetscScalar ea = fabs(e[Jj]);
    const PetscScalar growth = (_eCrit > 0 && !(ea < _eCrit)) ? 0.0 : (ea/g[Jj])*(1.0 - p[Jj]);
    r[Jj] = growth - p[Jj]/tauC(T[Jj]);
  }
  VecRestoreArrayRead(_state,&p); VecRestoreArrayRead(_e,&e); VecRestoreArrayRead(_gammaF,&g); VecRestoreArrayRead(in.T,&T);
  VecRestoreArray(_rate,&r);
  return ierr;
}


// strength factor 1/(1 + beta_f Phi) on every creep mechanism, in every direction or (fabric_anisotropic)
// on fault-parallel shear only
PetscErrorCode FabricState::pushToMaterial(PowerLaw& material)
{
  PetscErrorCode ierr = 0;
  ierr = VecCopy(_state,_F); CHKERRQ(ierr);
  ierr = VecScale(_F,_betaF); CHKERRQ(ierr);
  ierr = VecShift(_F,1.0); CHKERRQ(ierr);
  ierr = VecReciprocal(_F); CHKERRQ(ierr);
  if (_anisotropic == 1) { ierr = material.setXYStrengthFactor(_name,_F); CHKERRQ(ierr); }
  else { ierr = material.setStrengthFactor(_name,_F); CHKERRQ(ierr); }
  return ierr;
}


// the relaxation time of the law, 1/(e/gamma_f + 1/tau_c), at its shortest over the nodes
PetscErrorCode FabricState::computeMaxTimeStep(const BulkInputs& in, PetscScalar& maxDeltaT)
{
  PetscErrorCode ierr = 0;
  maxDeltaT = PETSC_MAX_REAL;
  if (_type != "transient") { return ierr; }
  ierr = drivingStrainRate(in); CHKERRQ(ierr);
  PetscInt Istart, Iend;
  ierr = VecGetOwnershipRange(_e,&Istart,&Iend); CHKERRQ(ierr);
  const PetscScalar *e, *g, *T;
  VecGetArrayRead(_e,&e); VecGetArrayRead(_gammaF,&g); VecGetArrayRead(in.T,&T);
  PetscScalar dt = PETSC_MAX_REAL;
  for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) { dt = PetscMin(dt, 1.0/(fabs(e[Jj])/g[Jj] + 1.0/tauC(T[Jj]))); }
  VecRestoreArrayRead(_e,&e); VecRestoreArrayRead(_gammaF,&g); VecRestoreArrayRead(in.T,&T);
  ierr = MPIU_Allreduce(&dt,&maxDeltaT,1,MPIU_SCALAR,MPIU_MIN,PETSC_COMM_WORLD); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode FabricState::clamp()
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


PetscErrorCode FabricState::writeContextExtra(PetscViewer& ascii, PetscViewer& viewer)
{
  PetscErrorCode ierr = 0;
  ierr = PetscViewerASCIIPrintf(ascii,"fabric_betaF = %.15e\n",_betaF); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"fabric_tauC0 = %.15e # (s)\n",_tauC0); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"fabric_QC = %.15e # (K)\n",_QC); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"fabric_TRef = %.15e # (K)\n",_TRef); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"fabric_eCrit = %.15e # (1/s, 0: no gate)\n",_eCrit); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"fabric_anisotropic = %i # (1: fault-parallel shear only)\n",_anisotropic); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PushGroup(viewer,"/fabric"); CHKERRQ(ierr);
  ierr = VecView(_gammaF,viewer); CHKERRQ(ierr);
  ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  return ierr;
}

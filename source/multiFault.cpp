#include "multiFault.hpp"

#define FILENAME "multiFault.cpp"

using namespace std;


PetscErrorCode createInteriorFaults(Domain& D, const vector<string>& names, map<string,PetscScalar>& positions,
  const PetscScalar faultTypeScale, const bool kinkLift, vector<Fault_qd*>& faults, vector<InteriorFaultLift*>& lifts)
{
  PetscErrorCode ierr = 0;
  for (size_t i = 0; i < names.size(); i++) {
    const string& name = names[i];
    PetscInt iRow = 0;
    PetscScalar yMid = 0;
    ierr = InteriorFaultLift::locate(D,positions[name],iRow,yMid); CHKERRQ(ierr);
    if (PetscAbsScalar(yMid - positions[name]) > 1e-9*PetscMax(1.0,PetscAbsScalar(yMid))) {
      PetscPrintf(PETSC_COMM_WORLD,"Note: interior fault %s is placed at y = %.6g km, midway between grid rows %i and %i (%s_y = %g).\n",
        name.c_str(),yMid,iRow,iRow+1,name.c_str(),positions[name]);
    }
    VecScatter *rows = NULL;
    ierr = D.makeRowScatter(iRow,rows); CHKERRQ(ierr);
    Fault_qd *f = new Fault_qd(D,*rows,faultTypeScale,name);
    f->setEtaScale(0.5); // between two identical half-spaces: eta = sqrt(mu*rho)/2
    faults.push_back(f);
    lifts.push_back(new InteriorFaultLift(D,iRow,kinkLift));
  }
  return ierr;
}


PetscErrorCode extendTimeIntInds(vector<string>& inds, vector<double>& scale, const vector<Fault_qd*>& faults)
{
  PetscErrorCode ierr = 0;
  if (inds.empty()) { return ierr; }
  vector<double> scale0 = scale;
  while (scale0.size() < inds.size()) { scale0.push_back(1.0); } // integrators default missing scales to 1
  vector<string> newInds;
  vector<double> newScale;
  for (size_t k = 0; k < inds.size(); k++) {
    vector<string> keys(1,inds[k]);
    if (inds[k] == "slip" || inds[k] == "psi") {
      keys.clear();
      for (size_t i = 0; i < faults.size(); i++) { keys.push_back(inds[k] == "slip" ? faults[i]->_slipKey : faults[i]->_psiKey); }
    }
    for (size_t j = 0; j < keys.size(); j++) {
      if (std::find(newInds.begin(),newInds.end(),keys[j]) == newInds.end()) { newInds.push_back(keys[j]); newScale.push_back(scale0[k]); }
    }
  }
  if (newInds != inds) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: timeIntInds = %s, for the slip and state of every fault.\n",vector2str(newInds).c_str()); CHKERRQ(ierr);
    inds = newInds;
    scale = newScale;
  }
  return ierr;
}


PetscErrorCode prepareCohesion(vector<string>& inds, vector<double>& scale, const vector<Fault_qd*>& faults,
  const string& timeIntegrator)
{
  PetscErrorCode ierr = 0;
  const vector<string> inds0 = inds;
  for (size_t i = 0; i < faults.size(); i++) {
    const Fault_qd *f = faults[i];
    if (f->_cohesionEvolution == "implicit" && timeIntegrator != "RK32_WBE" && timeIntegrator != "RK43_WBE") {
      PetscPrintf(PETSC_COMM_WORLD,"Error: %scohesionEvolution = implicit needs timeIntegrator = RK32_WBE or RK43_WBE.\n",f->_prefix.c_str());
      assert(0);
    }
    if (inds.empty() || !f->cohesionEvolves()) { continue; }
    while (scale.size() < inds.size()) { scale.push_back(1.0); } // integrators default missing scales to 1
    // scales from the input, not the state, so that a restart controls the steps as the run it continues
    PetscScalar m = f->_cohesionMaxVals.empty() ? 0.0 : *std::max_element(f->_cohesionMaxVals.begin(),f->_cohesionMaxVals.end());
    if (f->_cohesionReseal) { m = PetscMax(m,f->_cohesionLim); }
    if (!(m > 0)) { m = 1.0; }
    if (f->_cohesionEvolution == "explicit" && std::find(inds.begin(),inds.end(),f->_cohesionKey) == inds.end()) {
      inds.push_back(f->_cohesionKey); scale.push_back(m);
    }
    if (f->_cohesionReseal && std::find(inds.begin(),inds.end(),f->_cohesionMaxKey) == inds.end()) {
      inds.push_back(f->_cohesionMaxKey); scale.push_back(m);
    }
  }
  if (inds != inds0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: timeIntInds = %s, with the evolving cohesion.\n",vector2str(inds).c_str()); CHKERRQ(ierr);
  }
  return ierr;
}


PetscErrorCode printInteriorFaults(PetscViewer& ascii, const vector<string>& names, const int kinkLift,
  const vector<Fault_qd*>& faults, const vector<InteriorFaultLift*>& lifts)
{
  PetscErrorCode ierr = 0;
  if (names.empty()) { return ierr; }
  ierr = PetscViewerASCIIPrintf(ascii,"interiorFaults = %s\n",vector2str(names).c_str());CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(ascii,"interiorFaultKinkLift = %i # 1: B+ (second-order fault traction)\n",kinkLift);CHKERRQ(ierr);
  for (size_t i = 0; i < faults.size(); i++) {
    if (lifts[i] == NULL) { continue; }
    ierr = PetscViewerASCIIPrintf(ascii,"%s_y = %.15e # (km) midway between grid rows %i and %i\n",
      faults[i]->_name.c_str(),lifts[i]->_yFault,lifts[i]->_iRow,lifts[i]->_iRow+1);CHKERRQ(ierr);
  }
  return ierr;
}


PetscErrorCode writeInteriorFaultContext(PetscViewer& viewer, const vector<Fault_qd*>& faults, const vector<InteriorFaultLift*>& lifts)
{
  PetscErrorCode ierr = 0;
  for (size_t i = 0; i < faults.size(); i++) {
    if (lifts[i] == NULL) { continue; }
    const string g = faults[i]->group();
    const PetscScalar y = lifts[i]->_yFault, dy = 2.0*lifts[i]->_dyPlus;
    const PetscInt iRow = lifts[i]->_iRow;
    ierr = PetscViewerHDF5WriteAttribute(viewer, g.c_str(), "y", PETSC_SCALAR, &y); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, g.c_str(), "iRow", PETSC_INT, &iRow); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, g.c_str(), "dy", PETSC_SCALAR, &dy); CHKERRQ(ierr);
  }
  return ierr;
}


FaultSeries::FaultSeries() : _file(NULL), _probeDepth(-1) {}

FaultSeries::~FaultSeries()
{
  for (size_t i = 0; i < _probeWeights.size(); i++) { VecDestroy(&_probeWeights[i]); }
  for (size_t i = 0; i < _probeNode.size(); i++) { VecDestroy(&_probeNode[i]); }
  for (size_t i = 0; i < _weights.size(); i++) { VecDestroy(&_weights[i]); }
  if (_file != NULL) { PetscFClose(PETSC_COMM_WORLD,_file); }
}


// the depth of every fault node, and trapezoid weights for the depth integrals
PetscErrorCode FaultSeries::setup(const vector<Fault_qd*>& faults)
{
  PetscErrorCode ierr = 0;
  for (size_t i = 0; i < faults.size(); i++) {
    VecScatter toAll;
    Vec zAll, w;
    ierr = VecScatterCreateToAll(faults[i]->_z,&toAll,&zAll); CHKERRQ(ierr);
    ierr = VecScatterBegin(toAll,faults[i]->_z,zAll,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(toAll,faults[i]->_z,zAll,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    PetscInt n = 0;
    ierr = VecGetSize(zAll,&n); CHKERRQ(ierr);
    const PetscScalar *z;
    ierr = VecGetArrayRead(zAll,&z); CHKERRQ(ierr);
    _depths.push_back(vector<PetscScalar>(z,z+n));
    ierr = VecRestoreArrayRead(zAll,&z); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&toAll); CHKERRQ(ierr);
    ierr = VecDestroy(&zAll); CHKERRQ(ierr);

    const vector<PetscScalar>& zz = _depths.back();
    ierr = VecDuplicate(faults[i]->_slip,&w); CHKERRQ(ierr);
    PetscInt Istart, Iend;
    ierr = VecGetOwnershipRange(w,&Istart,&Iend); CHKERRQ(ierr);
    for (PetscInt Jj = Istart; Jj < Iend; Jj++) {
      PetscScalar dz = 0;
      if (Jj > 0) { dz += 0.5*(zz[Jj] - zz[Jj-1]); }
      if (Jj < n - 1) { dz += 0.5*(zz[Jj+1] - zz[Jj]); }
      ierr = VecSetValue(w,Jj,1e3*dz,INSERT_VALUES); CHKERRQ(ierr); // km -> m
    }
    ierr = VecAssemblyBegin(w); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(w); CHKERRQ(ierr);
    _weights.push_back(w);
  }
  return ierr;
}


PetscErrorCode FaultSeries::write(const string& outputDir, const bool append, const PetscInt stepCount,
  const PetscScalar time, const PetscScalar deltaT, const vector<Fault_qd*>& faults)
{
  PetscErrorCode ierr = 0;

  if (_file == NULL) { // (PetscFOpen stops with an error if the file cannot be opened)
    string name = outputDir + "faultSeries.txt";
    ierr = PetscFOpen(PETSC_COMM_WORLD,name.c_str(),append ? "a" : "w",&_file); CHKERRQ(ierr);
    if (!append) {
      ierr = PetscFPrintf(PETSC_COMM_WORLD,_file,"# step time(s) dt(s)"); CHKERRQ(ierr);
      for (size_t i = 0; i < faults.size(); i++) {
        const char *n = faults[i]->_name.c_str();
        ierr = PetscFPrintf(PETSC_COMM_WORLD,_file," %s_maxV(m/s) %s_zMaxV(km) %s_potencyRate(m^2/s) %s_potency(m^2)",n,n,n,n); CHKERRQ(ierr);
        if (_probeDepth >= 0) {
          for (size_t k = 0; k < _probeNames.size(); k++) {
            ierr = PetscFPrintf(PETSC_COMM_WORLD,_file," %s_%s@%.4gkm%s",n,_probeNames[k].c_str(),_probeDepth,_probeUnits[k].c_str()); CHKERRQ(ierr);
          }
          if (faults[i]->cohesionEvolves()) { ierr = PetscFPrintf(PETSC_COMM_WORLD,_file," %s_cohesion@%.4gkm(MPa)",n,_probeDepth); CHKERRQ(ierr); }
        }
      }
      ierr = PetscFPrintf(PETSC_COMM_WORLD,_file,"\n"); CHKERRQ(ierr);
    }
  }

  ierr = PetscFPrintf(PETSC_COMM_WORLD,_file,"%i %.15e %.6e",stepCount,time,deltaT); CHKERRQ(ierr);
  for (size_t i = 0; i < faults.size(); i++) {
    PetscInt loc = 0;
    PetscScalar maxV = 0, potRate = 0, pot = 0;
    ierr = VecMax(faults[i]->_slipVel,&loc,&maxV); CHKERRQ(ierr);
    ierr = VecDot(faults[i]->_slipVel,_weights[i],&potRate); CHKERRQ(ierr);
    ierr = VecDot(faults[i]->_slip,_weights[i],&pot); CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD,_file," %.9e %.6g %.9e %.9e",maxV,_depths[i][loc],potRate,pot); CHKERRQ(ierr);
    if (_probeDepth >= 0) {
      for (size_t k = 0; k < _probeFields.size(); k++) {
        PetscScalar v = 0;
        ierr = VecDot(*_probeFields[k],_probeWeights[i],&v); CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,_file," %.9e",v); CHKERRQ(ierr);
      }
      if (faults[i]->cohesionEvolves()) {
        PetscScalar c = 0;
        ierr = VecDot(faults[i]->_cohesion,_probeNode[i],&c); CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,_file," %.9e",c); CHKERRQ(ierr);
      }
    }
  }
  ierr = PetscFPrintf(PETSC_COMM_WORLD,_file,"\n"); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode FaultSeries::setProbes(Domain& D, const vector<Fault_qd*>& faults, const vector<InteriorFaultLift*>& lifts,
  const PetscScalar depth, const PetscScalar width)
{
  PetscErrorCode ierr = 0;
  if (_depths.empty()) { SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,"FaultSeries::setProbes needs setup first"); }
  // the grid depth nearest the requested one (the z of a fault's nodes, the same in every row)
  const vector<PetscScalar>& z = _depths[0];
  PetscInt iz0 = 0;
  for (size_t iz = 0; iz < z.size(); iz++) { if (fabs(z[iz] - depth) < fabs(z[iz0] - depth)) { iz0 = (PetscInt) iz; } }
  _probeDepth = z[iz0];
  const PetscInt Nz = D._Nz;
  for (size_t i = 0; i < faults.size(); i++) {
    // body nodes at depth z[iz0] within width of the fault, and always the rows next to it
    const bool interior = (i < lifts.size() && lifts[i] != NULL);
    const PetscScalar yf = interior ? lifts[i]->_yFault : 0.0;
    Vec w;
    ierr = VecDuplicate(D._y,&w); CHKERRQ(ierr);
    ierr = VecSet(w,0.0); CHKERRQ(ierr);
    PetscInt Istart, Iend;
    ierr = VecGetOwnershipRange(w,&Istart,&Iend); CHKERRQ(ierr);
    const PetscScalar *y;
    PetscScalar *wa;
    ierr = VecGetArrayRead(D._y,&y); CHKERRQ(ierr);
    ierr = VecGetArray(w,&wa); CHKERRQ(ierr);
    for (PetscInt Ii = Istart; Ii < Iend; Ii++) {
      const PetscInt iy = Ii/Nz, iz = Ii - iy*Nz;
      if (iz != iz0) { continue; }
      const bool next = interior ? (iy == lifts[i]->_iRow || iy == lifts[i]->_iRow + 1) : (iy == 0);
      if (next || fabs(y[Ii - Istart] - yf) <= width) { wa[Ii - Istart] = 1.0; }
    }
    ierr = VecRestoreArrayRead(D._y,&y); CHKERRQ(ierr);
    ierr = VecRestoreArray(w,&wa); CHKERRQ(ierr);
    PetscScalar count = 0;
    ierr = VecSum(w,&count); CHKERRQ(ierr);
    ierr = VecScale(w,1.0/count); CHKERRQ(ierr);
    _probeWeights.push_back(w);
    // the fault node at that depth
    Vec n;
    ierr = VecDuplicate(faults[i]->_slip,&n); CHKERRQ(ierr);
    ierr = VecSet(n,0.0); CHKERRQ(ierr);
    ierr = VecSetValue(n,iz0,1.0,INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(n); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(n); CHKERRQ(ierr);
    _probeNode.push_back(n);
  }
  return ierr;
}


PetscErrorCode FaultSeries::addProbe(const string& name, const string& unit, const Vec* field)
{
  _probeNames.push_back(name);
  _probeUnits.push_back(unit);
  _probeFields.push_back(field);
  return 0;
}


PetscErrorCode FaultSeries::flush()
{
  if (_file != NULL) { fflush(_file); }
  return 0;
}


FaultWorkKernel::FaultWorkKernel() : _work(NULL), _workFault(NULL) {}

FaultWorkKernel::~FaultWorkKernel()
{
  for (size_t i = 0; i < _Gw.size(); i++) { VecDestroy(&_Gw[i]); }
  VecDestroy(&_work);
  VecDestroy(&_workFault);
}


PetscErrorCode FaultWorkKernel::setup(Domain& D, SbpOps* sbp, const Vec& w, const vector<Fault_qd*>& faults,
  const vector<InteriorFaultLift*>& lifts)
{
  PetscErrorCode ierr = 0;
  _lifts = lifts;

  // Quadrature of the body: (H J) x = H (J .* x), with J = (dy/dq)(dz/dr) on a curvilinear grid.
  // H J is separable, (Hy Jy)(Hz Jz), and the y-part integrates dy/dq = Dq y exactly (SBP), so the
  // depth weights are Wz = (sum over rows of H J 1)/Ly; normalizing each kernel so that the sum over
  // rows of H J Gw equals Wz gives sum_i (Hy Jy)_i Gw_i = 1 at every depth.
  Vec Jd = NULL, hj = NULL, Wz = NULL, Nk = NULL;
  ierr = VecDuplicate(D._y,&_work); CHKERRQ(ierr);
  ierr = VecDuplicate(D._y,&hj); CHKERRQ(ierr);
  if (D._gridSpacingType == "variableGridSpacing") {
    Mat J,Jinv,qy,rz,yq,zr;
    ierr = sbp->getCoordTrans(J,Jinv,qy,rz,yq,zr); CHKERRQ(ierr);
    ierr = VecDuplicate(D._y,&Jd); CHKERRQ(ierr);
    ierr = MatGetDiagonal(J,Jd); CHKERRQ(ierr);
  }

  InteriorFaultLift *any = NULL; // the row sums use a lift's fault-to-rows scatter (the same for every fault)
  for (size_t k = 0; k < lifts.size(); k++) { if (lifts[k] != NULL) { any = lifts[k]; break; } }

  for (size_t k = 0; k < faults.size(); k++) {
    if (lifts[k] == NULL) { _Gw.push_back(NULL); continue; }
    if (_workFault == NULL) {
      ierr = VecDuplicate(faults[k]->_slip,&_workFault); CHKERRQ(ierr);
      ierr = VecDuplicate(faults[k]->_slip,&Wz); CHKERRQ(ierr);
      ierr = VecDuplicate(faults[k]->_slip,&Nk); CHKERRQ(ierr);
      // Wz = (sum over rows of H J 1) / Ly
      ierr = VecSet(_work,1.0); CHKERRQ(ierr);
      if (Jd != NULL) { ierr = VecPointwiseMult(_work,_work,Jd); CHKERRQ(ierr); }
      ierr = sbp->H(_work,hj); CHKERRQ(ierr);
      ierr = VecSet(Wz,0.0); CHKERRQ(ierr);
      ierr = VecScatterBegin(any->_fault2body,hj,Wz,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecScatterEnd(any->_fault2body,hj,Wz,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecScale(Wz,1.0/D._Ly); CHKERRQ(ierr);
    }
    // Gaussian in the distance from the fault, scaled by its value at the nearest rows (half a row
    // spacing away) so that a width far below the spacing does not underflow; the normalization
    // below removes any constant factor
    Vec g;
    ierr = VecDuplicate(D._y,&g); CHKERRQ(ierr);
    const PetscScalar yf = lifts[k]->_yFault, d0 = lifts[k]->_dyPlus;
    PetscInt Istart, Iend;
    ierr = VecGetOwnershipRange(g,&Istart,&Iend); CHKERRQ(ierr);
    const PetscScalar *y, *wA;
    PetscScalar *gA;
    ierr = VecGetArrayRead(D._y,&y); CHKERRQ(ierr);
    ierr = VecGetArrayRead(w,&wA); CHKERRQ(ierr);
    ierr = VecGetArray(g,&gA); CHKERRQ(ierr);
    for (PetscInt Jj = 0; Jj < Iend - Istart; Jj++) {
      if (!(wA[Jj] > 0)) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"the frictional-heat width w must be > 0 at every depth"); }
      const PetscScalar d = y[Jj] - yf;
      gA[Jj] = exp(-(d*d - d0*d0)/(2.*wA[Jj]*wA[Jj]));
    }
    ierr = VecRestoreArrayRead(D._y,&y); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(w,&wA); CHKERRQ(ierr);
    ierr = VecRestoreArray(g,&gA); CHKERRQ(ierr);

    // N = (sum over rows of H J g) / Wz; g /= N
    if (Jd != NULL) { ierr = VecPointwiseMult(_work,g,Jd); CHKERRQ(ierr); }
    else { ierr = VecCopy(g,_work); CHKERRQ(ierr); }
    ierr = sbp->H(_work,hj); CHKERRQ(ierr);
    ierr = VecSet(Nk,0.0); CHKERRQ(ierr);
    ierr = VecScatterBegin(lifts[k]->_fault2body,hj,Nk,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(lifts[k]->_fault2body,hj,Nk,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecPointwiseDivide(Nk,Nk,Wz); CHKERRQ(ierr);
    ierr = VecScatterBegin(lifts[k]->_fault2body,Nk,_work,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(lifts[k]->_fault2body,Nk,_work,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecPointwiseDivide(g,g,_work); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) g, "Gw"); CHKERRQ(ierr);
    _Gw.push_back(g);
  }
  ierr = VecDestroy(&Jd); CHKERRQ(ierr);
  ierr = VecDestroy(&hj); CHKERRQ(ierr);
  ierr = VecDestroy(&Wz); CHKERRQ(ierr);
  ierr = VecDestroy(&Nk); CHKERRQ(ierr);
  return ierr;
}


PetscErrorCode FaultWorkKernel::spread(const vector<Fault_qd*>& faults, Vec& Q)
{
  PetscErrorCode ierr = 0;
  ierr = VecSet(Q,0.0); CHKERRQ(ierr);
  for (size_t k = 0; k < faults.size(); k++) {
    if (_Gw[k] == NULL) { continue; }
    ierr = VecPointwiseMult(_workFault,faults[k]->_tauP,faults[k]->_slipVel); CHKERRQ(ierr);
    ierr = VecScatterBegin(_lifts[k]->_fault2body,_workFault,_work,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(_lifts[k]->_fault2body,_workFault,_work,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecPointwiseMult(_work,_work,_Gw[k]); CHKERRQ(ierr);
    ierr = VecAXPY(Q,1.0,_work); CHKERRQ(ierr);
  }
  return ierr;
}


PetscErrorCode FaultWorkKernel::writeContext(PetscViewer& viewer, const vector<Fault_qd*>& faults)
{
  PetscErrorCode ierr = 0;
  for (size_t k = 0; k < faults.size(); k++) {
    if (_Gw[k] == NULL) { continue; }
    ierr = PetscViewerHDF5PushGroup(viewer, faults[k]->group().c_str()); CHKERRQ(ierr);
    ierr = VecView(_Gw[k], viewer); CHKERRQ(ierr);
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
  }
  return ierr;
}

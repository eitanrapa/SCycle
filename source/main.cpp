
#include <petscts.h>
#include <petscviewerhdf5.h>
#include <string>
#include <petscdmda.h>

#include "genFuncs.hpp"
#include "spmat.hpp"
#include "domain.hpp"
#include "sbpOps.hpp"
#include "fault.hpp"
#include "linearElastic.hpp"
#include "powerLaw.hpp"
#include "pressureEq.hpp"
#include "heatEquation.hpp"
#include "linearElastic.hpp"
#include "powerLaw.hpp"
#include "problemContext.hpp"
#include "strikeSlip_linearElastic_qd.hpp"
#include "strikeSlip_linearElastic_fd.hpp"
#include "strikeSlip_linearElastic_qd_fd.hpp"
#include "strikeSlip_powerLaw_qd.hpp"
#include "strikeSlip_powerLaw_qd_fd.hpp"

using namespace std;


int runMMSTests(const char * inputFile)
{
  PetscErrorCode ierr = 0;

  PetscPrintf(PETSC_COMM_WORLD,"%-3s %-2s %-10s %-10s %-22s %-10s %-22s %-10s %-22s\n", "ord","Ny","dy","errL2u","log2(errL2u)","errL2gxy","log2(errL2gxy)", "errL2gxz","log2(errL2gxz)");

for (PetscInt Ny = 21; Ny < 82; Ny = (Ny - 1) * 2 + 1)
  {
    Domain d(inputFile,Ny,Ny);
    // Domain d(inputFile,Ny,1);
    //~ d.write();

    StrikeSlip_LinearElastic_qd m(d);
    ierr = m.writeContext(); CHKERRQ(ierr);
    ierr = m.integrate(); CHKERRQ(ierr);

    ierr = m.view(); CHKERRQ(ierr);
    ierr = m.measureMMSError();CHKERRQ(ierr);
  }

  return ierr;
}


// calculate Green's function mapping fault slip to surface displacement
// written to file "G"
// also write bcL and surfDisp into file


// Green's function mapping a unit displacement of each fault (left-boundary) node to surface
// displacement. Written as a dense Ny x Nz matrix to <outputDir>G_fault (PETSc binary).
int computeGreensFunction_fault(Domain& d)
{
  PetscErrorCode ierr = 0;
  PetscPrintf(PETSC_COMM_WORLD,"Running computeGreensFunction_fault\n");

  // linear elastic body: displacement prescribed on the fault and the right boundary, free surface on top and bottom
  LinearElastic le(d,"Dirichlet","Neumann","Dirichlet","Neumann");
  Mat A;
  le._sbp->getA(A);
  ierr = le.setupKSP(le._ksp,le._pc,A,le._linSolverSS); CHKERRQ(ierr);

  ierr = VecSet(le._bcT,0.0); CHKERRQ(ierr);
  ierr = VecSet(le._bcB,0.0); CHKERRQ(ierr);
  ierr = VecSet(le._bcR,0.0); CHKERRQ(ierr);

  Mat G;
  ierr = MatCreateDense(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,d._Ny,d._Nz,NULL,&G); CHKERRQ(ierr);
  ierr = MatSetUp(G); CHKERRQ(ierr);

  // each rank inserts the surface points it owns
  PetscInt sStart,sEnd,bStart,bEnd;
  ierr = VecGetOwnershipRange(le._surfDisp,&sStart,&sEnd); CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(le._bcL,&bStart,&bEnd); CHKERRQ(ierr);
  std::vector<PetscInt> rows(sEnd-sStart);
  for (PetscInt i = 0; i < sEnd-sStart; i++) { rows[i] = sStart + i; }

  // every rank takes part in every solve (the old loop ran over each rank's own bcL entries,
  // so collective calls were mismatched and parallel runs hung)
  for (PetscInt Ii = 0; Ii < d._Nz; Ii++) {
    PetscPrintf(PETSC_COMM_WORLD,"Ii = %i\n",Ii);
    ierr = VecSet(le._bcL,0.0); CHKERRQ(ierr);
    if (Ii >= bStart && Ii < bEnd) { ierr = VecSetValue(le._bcL,Ii,1.0,INSERT_VALUES); CHKERRQ(ierr); }
    ierr = VecAssemblyBegin(le._bcL); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(le._bcL); CHKERRQ(ierr);

    ierr = le._sbp->setRhs(le._rhs,le._bcL,le._bcR,le._bcT,le._bcB); CHKERRQ(ierr);
    ierr = KSPSolve(le._ksp,le._rhs,le._u); CHKERRQ(ierr);
    ierr = le.setSurfDisp(); CHKERRQ(ierr);

    const PetscScalar *si;
    ierr = VecGetArrayRead(le._surfDisp,&si); CHKERRQ(ierr);
    if (sEnd > sStart) { ierr = MatSetValues(G,sEnd-sStart,rows.data(),1,&Ii,si,INSERT_VALUES); CHKERRQ(ierr); }
    ierr = VecRestoreArrayRead(le._surfDisp,&si); CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(G,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(G,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = writeMat(G, d._outputDir + "G_fault"); CHKERRQ(ierr);
  ierr = MatDestroy(&G); CHKERRQ(ierr);
  return ierr;
}

// Green's function mapping a unit viscous strain in each body node (gVxy then gVxz, 2*Ny*Nz
// columns) to surface displacement; can map viscous strain rate to surface velocity.
// Written to <outputDir>G_offFault.h5 one column per time step; a rerun continues where the
// file ends.
int computeGreensFunction_offFault(Domain& d)
{
  PetscErrorCode ierr = 0;
  PetscPrintf(PETSC_COMM_WORLD,"Running computeGreensFunction_offFault\n");

  PetscViewer viewer;
  string outFileName = d._outputDir + "G_offFault.h5";
  PetscInt startIi = 0;
  if (doesFileExist(outFileName)) {
    PetscPrintf(PETSC_COMM_WORLD,"%s exists: continuing it\n",outFileName.c_str());
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, outFileName.c_str(), FILE_MODE_APPEND, &viewer);CHKERRQ(ierr);
    ierr = PetscViewerHDF5SetBaseDimension2(viewer, PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerHDF5PushTimestepping(viewer);                     CHKERRQ(ierr);
    PetscBool hasIi = PETSC_FALSE;
    ierr = PetscViewerHDF5HasAttribute(viewer, "surfDisp", "Ii", &hasIi); CHKERRQ(ierr);
    if (hasIi) { // otherwise the file is empty (a previous run stopped before its first column)
      ierr = PetscViewerHDF5ReadAttribute(viewer, "surfDisp", "Ii", PETSC_INT, NULL, &startIi); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"previous Ii = %i\n",startIi);
      startIi++;
    }
    ierr = PetscViewerHDF5SetTimestep(viewer, startIi); CHKERRQ(ierr);
  }
  else {
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, outFileName.c_str(), FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
    ierr = PetscViewerHDF5SetBaseDimension2(viewer, PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerHDF5PushTimestepping(viewer);                  CHKERRQ(ierr);
  }

  PowerLaw pl(d,"Dirichlet","Neumann","Dirichlet","Neumann");
  HeatEquation he(d);
  pl.updateTemperature(he._T);

  Mat A;
  pl._sbp->getA(A);
  ierr = pl.setupKSP(pl._ksp,pl._pc,A,pl._linSolverTrans); CHKERRQ(ierr);

  ierr = VecSet(pl._bcR,0.0); CHKERRQ(ierr);
  ierr = VecSet(pl._bcT,0.0); CHKERRQ(ierr);
  ierr = VecSet(pl._bcL,0.0); CHKERRQ(ierr);
  ierr = VecSet(pl._bcB,0.0); CHKERRQ(ierr);

  Vec viscSource;
  ierr = VecDuplicate(pl._gVxy,&viscSource); CHKERRQ(ierr);
  ierr = VecSet(viscSource,0.0); CHKERRQ(ierr);

  // one-entry Vec holding Ii, so a stopped job and its output stay consistent
  Vec test;
  ierr = VecCreate(PETSC_COMM_WORLD,&test); CHKERRQ(ierr);
  ierr = VecSetSizes(test,PETSC_DECIDE,1); CHKERRQ(ierr);
  ierr = VecSetFromOptions(test); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) test, "test"); CHKERRQ(ierr);

  const PetscInt N = d._Ny*d._Nz;
  PetscInt gStart,gEnd;
  ierr = VecGetOwnershipRange(pl._gVxy,&gStart,&gEnd); CHKERRQ(ierr);
  // all columns (the loop had been cut to startIi+10 while debugging)
  for (PetscInt Ii = startIi; Ii < 2*N; Ii++) {
    PetscPrintf(PETSC_COMM_WORLD,"Ii = %i...",Ii);
    ierr = VecSet(test,Ii); CHKERRQ(ierr);
    ierr = VecSet(pl._gVxy,0.0); CHKERRQ(ierr);
    ierr = VecSet(pl._gVxz,0.0); CHKERRQ(ierr);

    // a unit value in one entry of gVxy (Ii < N) or gVxz
    Vec g = (Ii < N) ? pl._gVxy : pl._gVxz;
    const PetscInt Jj = (Ii < N) ? Ii : Ii - N;
    if (Jj >= gStart && Jj < gEnd) { ierr = VecSetValue(g,Jj,1.0,INSERT_VALUES); CHKERRQ(ierr); }
    ierr = VecAssemblyBegin(g); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(g); CHKERRQ(ierr);

    // rhs source terms d/dy(mu*gVxy) + d/dz(mu*gVxz), then solve (also updates surfDisp)
    ierr = pl.computeViscStrainSourceTerms(viscSource); CHKERRQ(ierr);
    pl.setRHS();
    ierr = VecAXPY(pl._rhs,1.0,viscSource); CHKERRQ(ierr);
    ierr = pl.computeU(); CHKERRQ(ierr);

    if (Ii > startIi) {
      ierr = PetscViewerHDF5IncrementTimestep(viewer);                  CHKERRQ(ierr);
    }
    ierr = VecView(pl._surfDisp, viewer);                               CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, "surfDisp", "Ii", PETSC_INT, &Ii); CHKERRQ(ierr);
    ierr = VecView(test, viewer);                                       CHKERRQ(ierr);
    ierr = flushHDF5Viewer(viewer);                                     CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"finished.\n");
  }

  ierr = VecDestroy(&viscSource); CHKERRQ(ierr);
  ierr = VecDestroy(&test); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  return ierr;
}


int runEqCycle(Domain& d)
{
  PetscErrorCode ierr = 0;
  // report number of processors
  PetscMPIInt numCores;
  MPI_Comm_size(PETSC_COMM_WORLD,&numCores);
  PetscPrintf(PETSC_COMM_WORLD,"Total number of processors: %i\n\n",numCores);

  if (d._bulkDeformationType == "linearElastic") {
    ProblemContext *m = NULL;
    if (d._momentumBalanceType == "quasidynamic") { m = new StrikeSlip_LinearElastic_qd(d); }
    else if (d._momentumBalanceType == "dynamic") { m = new StrikeSlip_LinearElastic_fd(d); }
    else if (d._momentumBalanceType == "quasidynamic_and_dynamic") { m = new StrikeSlip_LinearElastic_qd_fd(d); }
    else {
      PetscPrintf(PETSC_COMM_WORLD,"momentumBalanceType = %s is not supported with bulkDeformationType = linearElastic\n",d._momentumBalanceType.c_str());
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"unsupported momentumBalanceType");
    }

    if (d._restartFromChkpt == 0) { ierr = m->writeContext(); CHKERRQ(ierr); }
    PetscPrintf(PETSC_COMM_WORLD,"\n\n\n");
    ierr = m->integrate(); CHKERRQ(ierr);
    ierr = m->view(); CHKERRQ(ierr);
    delete m;
  }

  else if (d._bulkDeformationType == "powerLaw") {
    ProblemContext *m = NULL;
    if (d._momentumBalanceType == "quasidynamic") { m = new StrikeSlip_PowerLaw_qd(d); }
    else if (d._momentumBalanceType == "quasidynamic_and_dynamic") { m = new StrikeSlip_PowerLaw_qd_fd(d); }
    else {
      PetscPrintf(PETSC_COMM_WORLD,"momentumBalanceType = %s is not supported with bulkDeformationType = powerLaw\n",d._momentumBalanceType.c_str());
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"unsupported momentumBalanceType");
    }

    if (d._restartFromChkpt == 0) { ierr = m->writeContext(); CHKERRQ(ierr); }
    PetscPrintf(PETSC_COMM_WORLD,"\n\n\n");

    if (d._systemEvolutionType == "steadyStateIts") { ierr = m->integrateSS(); CHKERRQ(ierr); }
    if (d._systemEvolutionType == "transient") {
      ierr = m->initiateIntegrand(); CHKERRQ(ierr);
      ierr = m->integrate(); CHKERRQ(ierr);
    }

    ierr = m->view(); CHKERRQ(ierr);
    delete m;
  }
  else {
    PetscPrintf(PETSC_COMM_WORLD,"bulkDeformationType = %s is not supported\n",d._bulkDeformationType.c_str());
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"unsupported bulkDeformationType");
  }

  return ierr;
}


int main(int argc,char **args)
{
  PetscErrorCode ierr = 0;
  ierr = PetscInitialize(&argc,&args,NULL,NULL);
  if (ierr) { return 1; }

  const char * inputFile;

  if (argc > 1) { inputFile = args[1]; }
  else { inputFile = "init.in"; }

  {
    Domain d(inputFile);
    if (d._isMMS) { ierr = runMMSTests(inputFile); } // builds its own domain for each resolution
    else if (d._computeGreensFunction_fault) { ierr = computeGreensFunction_fault(d); }
    else if (d._computeGreensFunction_offFault) { ierr = computeGreensFunction_offFault(d); }
    else { ierr = runEqCycle(d); }
  }

  PetscFinalize();
  return (ierr != 0); // nonzero exit status if the run reported a PETSc error
}

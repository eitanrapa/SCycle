#include "interiorFaultLift.hpp"

#define FILENAME "interiorFaultLift.cpp"

using namespace std;


InteriorFaultLift::InteriorFaultLift(Domain& D, const PetscInt iRow, const bool kinkLift)
: _D(&D), _iRow(iRow), _yFault(0.), _kinkLift(kinkLift),
  _step(NULL), _band(NULL), _Ut(NULL), _U(NULL), _c(NULL), _cFault(NULL), _Aq(NULL), _halfDy2(NULL), _dy(NULL), _work(NULL), _work2(NULL),
  _tauPlus(NULL), _cPlus(NULL), _muPlus(NULL), _dyPlus(0.),
  _fault2body(NULL), _rowMinus(NULL), _rowPlus(NULL)
{
  PetscErrorCode ierr = 0;

  // y of every grid row on every rank (y along the top boundary is D._z0): the fault position
  {
    VecScatter toAll; Vec yAll;
    VecScatterCreateToAll(D._z0,&toAll,&yAll);
    VecScatterBegin(toAll,D._z0,yAll,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(toAll,D._z0,yAll,INSERT_VALUES,SCATTER_FORWARD);
    const PetscScalar *y;
    VecGetArrayRead(yAll,&y);
    _yFault = 0.5*(y[iRow] + y[iRow+1]);
    _dyPlus = y[iRow+1] - _yFault;
    VecRestoreArrayRead(yAll,&y);
    VecScatterDestroy(&toAll);
    VecDestroy(&yAll);
  }

  // body fields: the step (1 on rows past the fault), the band of rows whose D_y stencil crosses
  // the fault (interior stencils reach 2 rows for 4th order, 1 for 2nd), and with B+ the distance
  // from the fault
  VecDuplicate(D._y,&_step);
  VecDuplicate(D._y,&_band);
  VecDuplicate(D._y,&_Ut); VecSet(_Ut,0.);
  VecDuplicate(D._y,&_U);  VecSet(_U,0.);
  VecDuplicate(D._y,&_work);
  VecDuplicate(D._y,&_work2);
  PetscInt Istart, Iend;
  VecGetOwnershipRange(_step,&Istart,&Iend);
  PetscScalar *s;
  VecGetArray(_step,&s);
  for (PetscInt Ii = Istart; Ii < Iend; Ii++) { s[Ii-Istart] = (Ii/D._Nz > iRow) ? 1.0 : 0.0; }
  VecRestoreArray(_step,&s);
  VecGetArray(_band,&s);
  for (PetscInt Ii = Istart; Ii < Iend; Ii++) {
    const PetscInt iy = Ii/D._Nz;
    s[Ii-Istart] = (iy >= iRow - 2 && iy <= iRow + 3) ? 1.0 : 0.0;
  }
  VecRestoreArray(_band,&s);
  if (_kinkLift) {
    VecDuplicate(D._y,&_c); VecSet(_c,0.);
    VecDuplicate(D._y,&_dy);
    VecCopy(D._y,_dy);
    VecShift(_dy,-_yFault);                     // y - yf
    VecDuplicate(D._y,&_halfDy2);
    VecPointwiseMult(_halfDy2,_dy,_dy);
    VecScale(_halfDy2,0.5);                     // (y - yf)^2 / 2
  }

  // scatter copying a fault vector into every row: body entry Ii takes fault entry Ii % Nz
  {
    const PetscInt n = Iend - Istart;
    PetscInt *from;
    PetscMalloc1(n,&from);
    for (PetscInt i = 0; i < n; i++) { from[i] = (Istart + i) % D._Nz; }
    IS isFrom, isTo;
    ISCreateGeneral(PETSC_COMM_WORLD,n,from,PETSC_OWN_POINTER,&isFrom);
    ISCreateStride(PETSC_COMM_WORLD,n,Istart,1,&isTo);
    VecScatterCreate(D._y0,isFrom,_Ut,isTo,&_fault2body);
    ISDestroy(&isFrom);
    ISDestroy(&isTo);
  }

  // the rows on either side of the fault
  ierr = D.makeRowScatter(iRow,_rowMinus); assert(ierr == 0);
  ierr = D.makeRowScatter(iRow+1,_rowPlus); assert(ierr == 0);
  VecDuplicate(D._y0,&_tauPlus);
  VecDuplicate(D._y0,&_cPlus); VecSet(_cPlus,0.);
  VecDuplicate(D._y0,&_muPlus);
}


InteriorFaultLift::~InteriorFaultLift()
{
  VecDestroy(&_step);
  VecDestroy(&_band);
  VecDestroy(&_Ut);
  VecDestroy(&_U);
  VecDestroy(&_c);
  VecDestroy(&_cFault);
  VecDestroy(&_Aq);
  VecDestroy(&_halfDy2);
  VecDestroy(&_dy);
  VecDestroy(&_work);
  VecDestroy(&_work2);
  VecDestroy(&_tauPlus);
  VecDestroy(&_cPlus);
  VecDestroy(&_muPlus);
  VecScatterDestroy(&_fault2body);
}


// row iRow with y between rows iRow and iRow + 1, and the midpoint of those rows
PetscErrorCode InteriorFaultLift::locate(Domain& D, const PetscScalar yWanted, PetscInt& iRow, PetscScalar& yMid)
{
  PetscErrorCode ierr = 0;

  // y of every grid row, on every rank (y along the top boundary is D._z0)
  VecScatter toAll; Vec yAll;
  ierr = VecScatterCreateToAll(D._z0,&toAll,&yAll); CHKERRQ(ierr);
  ierr = VecScatterBegin(toAll,D._z0,yAll,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(toAll,D._z0,yAll,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  const PetscScalar *y;
  ierr = VecGetArrayRead(yAll,&y); CHKERRQ(ierr);

  iRow = -1;
  for (PetscInt i = 0; i < D._Ny - 1; i++) {
    if (y[i] <= yWanted && yWanted < y[i+1]) { iRow = i; break; }
  }
  if (iRow >= 0) { yMid = 0.5*(y[iRow] + y[iRow+1]); }
  const PetscScalar yFirst = y[0], yLast = y[D._Ny-1];

  ierr = VecRestoreArrayRead(yAll,&y); CHKERRQ(ierr);
  ierr = VecScatterDestroy(&toAll); CHKERRQ(ierr);
  ierr = VecDestroy(&yAll); CHKERRQ(ierr);

  if (iRow < 0) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: interior fault position y = %g km is outside the grid [%g, %g] km\n",yWanted,yFirst,yLast);
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"interior fault outside the grid");
  }
  if (iRow < 6 || iRow + 1 > D._Ny - 7) {
    PetscPrintf(PETSC_COMM_WORLD,"Error: an interior fault must lie at least 6 grid rows from the left and right boundaries;\n"
      "       y = %g km falls between rows %i and %i of %i\n",yWanted,iRow,iRow+1,D._Ny);
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"interior fault too close to a boundary");
  }
  return ierr;
}


// Ut = slip (+ Kt with B+) on every row, U = Ut on rows past the fault
PetscErrorCode InteriorFaultLift::setSlip(const Vec& slip, const Mat& A, const Vec* source)
{
  PetscErrorCode ierr = 0;
  ierr = VecScatterBegin(_fault2body,slip,_Ut,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(_fault2body,slip,_Ut,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  if (_kinkLift) {
    if (_Aq == NULL) { // A q on row _iRow: the y-part of A on the quadratic, mu times the row weights
      ierr = VecDuplicate(_tauPlus,&_Aq); CHKERRQ(ierr);
      ierr = VecDuplicate(_tauPlus,&_cFault); CHKERRQ(ierr);
      ierr = MatMult(A,_halfDy2,_work); CHKERRQ(ierr);
      ierr = VecScatterBegin(*_rowMinus,_work,_Aq,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(*_rowMinus,_work,_Aq,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    }
    // c = -(A delta - [S])/(A q) on row _iRow, the jump of w_yy across the fault, copied to every row
    ierr = MatMult(A,_Ut,_work); CHKERRQ(ierr);
    ierr = VecScatterBegin(*_rowMinus,_work,_cFault,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(*_rowMinus,_work,_cFault,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    if (source != NULL) { // [S] = S(_iRow + 1) - S(_iRow)
      ierr = VecScatterBegin(*_rowPlus,*source,_tauPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(*_rowPlus,*source,_tauPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterBegin(*_rowMinus,*source,_cPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(*_rowMinus,*source,_cPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecAXPY(_tauPlus,-1.0,_cPlus); CHKERRQ(ierr);
      ierr = VecAXPY(_cFault,-1.0,_tauPlus); CHKERRQ(ierr);
    }
    ierr = VecPointwiseDivide(_cFault,_cFault,_Aq); CHKERRQ(ierr);
    ierr = VecScale(_cFault,-1.0); CHKERRQ(ierr);
    ierr = VecScatterBegin(_fault2body,_cFault,_c,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(_fault2body,_cFault,_c,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    // Kt = c (y - yf)^2 / 2, added to the slip
    ierr = VecPointwiseMult(_work,_c,_halfDy2); CHKERRQ(ierr);
    ierr = VecAXPY(_Ut,1.0,_work); CHKERRQ(ierr);
  }
  ierr = VecPointwiseMult(_U,_step,_Ut); CHKERRQ(ierr);
  return ierr;
}


// rhs += A U - step .* (A Ut): the right-hand side for the displacement u with this fault's jump
PetscErrorCode InteriorFaultLift::addToRhs(const Mat& A, Vec& rhs)
{
  PetscErrorCode ierr = 0;
  ierr = MatMult(A,_U,_work); CHKERRQ(ierr);
  ierr = VecAXPY(rhs,1.0,_work); CHKERRQ(ierr);
  ierr = MatMult(A,_Ut,_work); CHKERRQ(ierr);
  ierr = VecPointwiseMult(_work,_step,_work); CHKERRQ(ierr);
  ierr = VecAXPY(rhs,-1.0,_work); CHKERRQ(ierr);
  return ierr;
}


// On the band rows, sxy = mu D_y u becomes mu D_y (u - U) + mu c (y - yf) step: the y-derivative
// of the continuous part, plus the exact dK/dy of the B+ lift. Off the band D_y U = 0 (U is constant
// in y there), and mu D_y u is kept: it is more accurate than D_y applied to the large quadratic K.
PetscErrorCode InteriorFaultLift::correctStress(Vec& sxy, SbpOps* sbp, const Vec& mu)
{
  PetscErrorCode ierr = 0;
  ierr = sbp->muxDy(_U,_work); CHKERRQ(ierr);           // mu D_y U
  if (_kinkLift) {
    ierr = VecPointwiseMult(_work2,_c,_dy); CHKERRQ(ierr);
    ierr = VecPointwiseMult(_work2,_work2,_step); CHKERRQ(ierr);
    ierr = VecPointwiseMult(_work2,_work2,mu); CHKERRQ(ierr);
    ierr = VecAXPY(_work,-1.0,_work2); CHKERRQ(ierr);   // - mu c (y - yf) past the fault
  }
  ierr = VecPointwiseMult(_work,_work,_band); CHKERRQ(ierr);
  ierr = VecAXPY(sxy,-1.0,_work); CHKERRQ(ierr);
  return ierr;
}


// The strain version of correctStress: gxy = D_y u becomes D_y (u - U) + c (y - yf) step on the band.
PetscErrorCode InteriorFaultLift::correctStrain(Vec& gxy, SbpOps* sbp)
{
  PetscErrorCode ierr = 0;
  ierr = sbp->Dy(_U,_work); CHKERRQ(ierr);              // D_y U
  if (_kinkLift) {
    ierr = VecPointwiseMult(_work2,_c,_dy); CHKERRQ(ierr);
    ierr = VecPointwiseMult(_work2,_work2,_step); CHKERRQ(ierr);
    ierr = VecAXPY(_work,-1.0,_work2); CHKERRQ(ierr);   // - c (y - yf) past the fault
  }
  ierr = VecPointwiseMult(_work,_work,_band); CHKERRQ(ierr);
  ierr = VecAXPY(gxy,-1.0,_work); CHKERRQ(ierr);
  return ierr;
}


// Fault traction: the average of sxy (after correctStress) on the two rows, less this lift's own
// dK/dy there (mu c (y - yf) on row iRow + 1 only, so half of it in the average); what remains is
// the average of mu D_y w, which includes the stress of every other fault.
PetscErrorCode InteriorFaultLift::traction(const Vec& sxy, const Vec& mu, Vec& tau)
{
  PetscErrorCode ierr = 0;
  ierr = VecScatterBegin(*_rowMinus,sxy,tau,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(*_rowMinus,sxy,tau,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterBegin(*_rowPlus,sxy,_tauPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(*_rowPlus,sxy,_tauPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecAXPY(tau,1.0,_tauPlus); CHKERRQ(ierr);
  ierr = VecScale(tau,0.5); CHKERRQ(ierr);
  if (_kinkLift) {
    ierr = VecScatterBegin(*_rowPlus,_c,_cPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(*_rowPlus,_c,_cPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterBegin(*_rowPlus,mu,_muPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(*_rowPlus,mu,_muPlus,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecPointwiseMult(_cPlus,_cPlus,_muPlus); CHKERRQ(ierr);
    ierr = VecAXPY(tau,-0.5*_dyPlus,_cPlus); CHKERRQ(ierr);
  }
  return ierr;
}

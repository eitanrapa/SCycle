#ifndef INTERIORFAULTLIFT_HPP_INCLUDED
#define INTERIORFAULTLIFT_HPP_INCLUDED

#include <petscksp.h>
#include "domain.hpp"
#include "sbpOps.hpp"

/*
 * Interior fault in the elastic momentum balance by a jump-corrected right-hand side
 * ("Route B" of docs/TWO_FAULT_DESIGN.md).
 *
 * The fault lies midway between grid rows iRow and iRow + 1, and its slip delta(z) is the jump of
 * the displacement u across it. With the step field U = delta(z) on rows > iRow (0 elsewhere) and
 * its smooth extension Ut = delta(z) on every row, the displacement solves
 *
 *     A u = rhs_BC + J,    J = A U - step .* (A Ut),
 *
 * with the unchanged operator A. J is nonzero only on rows whose y-stencil crosses the fault: the
 * y-derivatives annihilate Ut, and the z-derivatives and the boundary terms act within a row or see
 * the same field in U and Ut. u is the physical displacement, jump included. w = u - U is
 * continuous across the fault, so the y-strain on rows whose D_y stencil crosses the fault is taken
 * from w (correctStress); elsewhere, and for z-strains, u is used as it is. The fault traction, the
 * average of mu D_y w on rows iRow and iRow + 1, is first order in the grid spacing when the slip
 * varies with depth, because w_yy then jumps by c = -(mu delta')'/mu across the fault.
 *
 * "B+" (kinkLift, the default): the lift also carries K = c(z) (y - yf)^2 / 2 on rows past the
 * fault (Kt on every row), which removes that jump: w = u - U - K is C2 and the traction becomes
 * second order. Near the fault the stress is mu D_y w + mu c (y - yf) on rows past it, and the fault
 * traction is taken from mu D_y w alone (dK/dy vanishes at the fault). Far from it, mu D_y u is
 * used: D_y of the quadratic K is not exact on stretched grids or at the boundary closures, and K
 * grows with the distance from the fault. c is computed each time the slip changes from the
 * operator itself, so that A Ut vanishes on the fault rows: c = -(A delta)/(A q) on row iRow, with
 * delta and q = (y - yf)^2 / 2 extended to every row. A acts on delta through its z-part (with the
 * top and bottom boundary terms) and on q through its y-part, so the ratio is the discrete
 * -(mu delta')'/mu, consistent with A up to the free surface.
 *
 * Requires mu continuous across the fault and the fault at least 6 rows from the y-boundaries
 * (the boundary closures of the y-operators span 6 rows).
 */
class InteriorFaultLift
{
private:
  // disable default copy constructor and assignment operator
  InteriorFaultLift(const InteriorFaultLift& that);
  InteriorFaultLift& operator=(const InteriorFaultLift& rhs);

public:
  Domain      *_D;
  PetscInt     _iRow;        // the fault lies between grid rows _iRow and _iRow + 1
  PetscScalar  _yFault;      // its position (km), the midpoint of the two rows
  bool         _kinkLift;    // B+: also lift the jump in w_yy (second-order traction)
  Vec          _step;        // body field: 1 on rows > _iRow, 0 elsewhere
  Vec          _band;        // body field: 1 on rows _iRow - 2 .. _iRow + 3 (D_y stencils that cross the fault), 0 elsewhere
  Vec          _Ut, _U;      // lift on every row (slip + Kt) and on rows > _iRow only (body fields)
  Vec          _c;           // B+: c(z) = -(mu delta')'/mu on every row (body field)
  Vec          _cFault, _Aq; // B+: c on the fault nodes, and (A q) on row _iRow (computed once)
  Vec          _halfDy2, _dy; // B+: (y - yf)^2 / 2 and y - yf (body fields)
  Vec          _work, _work2; // body work vectors
  Vec          _tauPlus, _cPlus, _muPlus; // fault-size work vectors (row _iRow + 1)
  PetscScalar  _dyPlus;      // y of row _iRow + 1 minus yf
  VecScatter   _fault2body;  // copies a fault vector (size Nz) into every row of a body field
  VecScatter  *_rowMinus, *_rowPlus; // body field -> rows _iRow and _iRow + 1 (owned by the Domain)

  InteriorFaultLift(Domain& D, const PetscInt iRow, const bool kinkLift = true);
  ~InteriorFaultLift();

  // set U and Ut from the slip (size Nz); with B+ also c, from the momentum-balance operator A
  PetscErrorCode setSlip(const Vec& slip, const Mat& A);
  PetscErrorCode addToRhs(const Mat& A, Vec& rhs);   // rhs += A U - step .* (A Ut)
  // sxy = mu D_y u on entry; on the band rows it becomes mu D_y (u - U) (+ mu c (y - yf) past the fault)
  PetscErrorCode correctStress(Vec& sxy, SbpOps* sbp, const Vec& mu);
  // fault traction from the physical sxy: average of the two rows, less this lift's own dK/dy term
  PetscErrorCode traction(const Vec& sxy, const Vec& mu, Vec& tau);

  // the row iRow such that y lies between rows iRow and iRow + 1, and the midpoint of those rows;
  // stops with a message if the fault would be closer than 6 rows to the left or right boundary
  static PetscErrorCode locate(Domain& D, const PetscScalar y, PetscInt& iRow, PetscScalar& yMid);
};

#endif

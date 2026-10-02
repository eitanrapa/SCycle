# SCycle audit, October 2026

## 1. Summary

This document records an audit of SCycle, the C++/PETSc code for 2D antiplane earthquake-cycle simulations, and the fixes made on branch `audit/fixes-2026-10`. The audit started from an unmodified clone of upstream master (github.com/kali-allison/SCycle, last upstream commit `74a132f` of 19 December 2024).

The branch adds 67 commits to master: 58 fixes, 3 tool commits, and one commit each for the build, the repository contents, the example inputs, a code comment, the SEAS BP1 benchmark with an input-key checker, and this documentation. Together they fix 94 distinct findings (32 Critical, 25 High, 25 Medium, 12 Low), listed in section 3. Three are fixed only in part (A-35, A-45, A-57), and the MATLAB tools of A-93 were not run. Section 5 lists 16 open items and documented limitations. Each commit message records its finding: what was wrong, the failure it caused, the fix and the test evidence (`git show <hash>`). A few spot tests from early in the audit appear only in this document; their entries say so.

How the fixes were verified:

- `tools/regress.sh` runs `examples/ex1.in` (1D) and `examples/ex2.in` (2D) and compares every HDF5 output file with a stored baseline using `h5diff`, bit for bit unless a tolerance is set. Unless its message says otherwise, each code commit was checked this way. The persistent baseline is in `data/regress-baseline/` (ignored by git). It was generated from commit `35ba52f`, is bit-identical to the output of the build at `06bada2`, and was still bit-identical at `aa53679`.
- `tools/mms.in` runs the manufactured-solution (MMS) convergence test of the linear-elastic quasi-dynamic momentum balance on Ny = Nz = 21, 41 and 81.
- `tools/checkkeys.py` lists keys in an input file that no part of the code reads.
- Targeted spot runs for configurations that ex1 and ex2 do not exercise. Each commit message describes them; the evidence lines below quote their numbers.

Build: PETSc 3.26.0 debug build (`PETSC_ARCH=arch-darwin-c-debug`) on macOS, real double-precision scalars, 32-bit indices, with MUMPS, HYPRE and HDF5. The regression baseline uses one MPI rank; entries say when a test used two.

## 2. Severity scale

| Severity | Meaning |
|---|---|
| Critical | Crash, hang, or corrupted or lost output in a supported configuration. |
| High | Silently wrong results. |
| Medium | Wrong in rare configurations, wasted resources, or misleading output. |
| Low | Hygiene, documentation, tooling. |

How the scale is applied:

- An assertion failure or a PETSc error that stops the run counts as a crash. A requested output field or restart file that is missing, truncated or unreadable counts as lost output.
- "Latent" means the defect was found by reading the code and the failure was not reproduced. The entry keeps the severity of the failure the defect would cause.
- "Partly fixed" entries have a matching open item in section 5.

## 3. Fixed findings

| Group | Entries | Critical | High | Medium | Low |
|---|---|---|---|---|---|
| 3.1 Crashes, undefined behaviour and memory | A-01 to A-22 | 17 | 0 | 5 | 0 |
| 3.2 Physics and numerics (wrong results) | A-23 to A-47 | 1 | 17 | 4 | 3 |
| 3.3 Time integration, checkpoints and I/O | A-48 to A-68 | 11 | 3 | 6 | 1 |
| 3.4 Input handling and error reporting | A-69 to A-79 | 3 | 3 | 4 | 1 |
| 3.5 Build, repository and tooling | A-80 to A-87 | 0 | 0 | 1 | 7 |
| 3.6 Examples, benchmarks and post-processing | A-88 to A-94 | 0 | 2 | 5 | 0 |
| Total | 94 | 32 | 25 | 25 | 12 |

"No run reported" means the commit message describes the fix without a dedicated test; the ex1/ex2 regression comparison still applied to it.

### 3.1 Crashes, undefined behaviour and memory

**A-01. Unsupported problem types ran through an uninitialized pointer.** Critical. `e9cdc49`. `source/main.cpp`.
- Defect: `runEqCycle()` assigned its `ProblemContext *m` only inside independent `if` blocks, so combinations that pass Domain's checks but have no problem class (power law with `momentumBalanceType = dynamic`, or `momentumBalanceType = steadyStateIts`, which the manual lists) called methods through an uninitialized pointer.
- Evidence: no run reported; the cases are now an `else if` chain that ends in a PETSc error naming the combination (checked in the code).

**A-02. BracketedNewton evaluated its first residual at an unset member.** Critical (latent). `b7aa7b7`. `source/rootFinder.cpp`.
- Defect: `findRoot()` checked the initial guess at the member `_x`, which was never initialized (the iterate was a local variable that shadowed it), so every rate-and-state velocity solve (one per fault node per call) read stack garbage; a NaN or Inf there trips the asserts.
- Evidence: no run numbers; results are unchanged in practice, because a guess within tolerance reached the same answer through the main loop.

**A-03. Flash heating gave NaN at locked fault nodes.** Critical. `5c2606b`. `source/fault.cpp`.
- Defect: `flashHeating_psi()` evaluates log(|V|/v0), which is -inf at V = 0, so `stateLaw = flashHeating` with any region held locked (`lockedVals > 0.5`) produced a NaN state rate and an assertion failure.
- Evidence: spot test during the audit, not recorded in the commit message: flash heating with a locked depth range aborts on the NaN assert with the old code and runs 50 steps with the fix. The restored guard is the one in `slipLaw_psi`, and the rate's limit as V -> 0 is 0, so returning 0 is exact.

**A-04. Fault_fd::setPhi re-acquired arrays instead of restoring them.** Critical. `0802469`. `source/fault.cpp`.
- Defect: the cleanup block called `VecGetArrayRead` a second time on six vectors instead of `VecRestoreArrayRead`, leaving them locked after every fully dynamic step; PETSc debug builds report this as an error at the next write access.
- Evidence: no dedicated test; later fully dynamic runs on the same debug build complete (for example 14958 dynamic steps in `8954864`).

**A-05. The pressure state handed to the integrator was freed.** Critical. `5c87a81`. `source/pressureEq.cpp`.
- Defect: `PressureEq::initiateIntegrand()` inserted new Vecs for pressure and permeability into the integrand maps and then destroyed them, so every run with `hydraulicCoupling != no` had the integrator duplicate and update freed memory.
- Evidence: spot test during the audit, not recorded in the commit message: ex2 with coupled pore pressure stops with PETSc "Corrupt argument: Object already free" on the old code and runs with the fix. `04298a8` later reproduced its own failure on the original code plus this fix only.

**A-06. Grain-size flags and output strides unset in the power-law qd_fd class.** Critical. `4931d1c`. `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: `StrikeSlip_PowerLaw_qd_fd` never initialized `_grainSizeEvCoupling`, `_grainSizeEvCouplingSS` or its six per-phase output strides, so an input without `grainSizeEvCoupling` failed the assert in `checkInput` and an omitted stride key left the output cadence to stack garbage.
- Evidence: no run reported; defaults now follow the sibling classes ("no" and stride 10).

**A-07. Fully dynamic runs wrote to time Vecs that were never created.** Critical. `a54f046`. `source/strikeSlip_linearElastic_fd.cpp`.
- Defect: `_time1DVec`, `_dtime1DVec`, `_time2DVec` and `_dtime2DVec` stayed NULL, so every `momentumBalanceType = dynamic` run failed at its first output (step 0).
- Evidence: spot test during the audit, not recorded in the commit message: ex2 with `momentumBalanceType = dynamic` stops with a PETSc null-pointer error at the first output on the old code and runs with the fix. The mode is also among the runs checked with `-objects_dump` in `8916df7`.

**A-08. Members read before being set in six classes.** Critical. `bcf2fed`. `source/linearElastic.cpp`, `source/fault.cpp`, `source/sbpOps_m_varGrid.cpp`, `source/grainSizeEvolution.cpp`, `source/strikeSlip_linearElastic_qd.cpp`, `source/strikeSlip_powerLaw_qd.cpp`.
- Defect: LinearElastic passed unset tolerances (`_akspTol`, `_rkspTol`) to the CG solver; Fault_fd had no default `_deltaT` and evaluated an unset prestress perturbation even with the default `timeMode`, which is NaN when the width holds 0; GrainSizeEvolution asserted on an unset `_c`; `SbpOps_m_varGrid` tested unset `_y`/`_z` against NULL; and two mediators read `_deltaT` before setting it.
- Evidence: defaults added (tolerances 1e-10, perturbation parameters 0 and widths 1, `_c` = 0 so a missing `grainSizeEv_c` fails with a message); ex1 and ex2 bit-identical.

**A-09. Radiogenic decay length unset; ex4 produced NaN temperatures.** Critical. `0e844ac`. `source/heatEquation.cpp`.
- Defect: `_Lrad` was set only from `he_Lrad`, and without `he_A0Vals` the source A0 exp(-z/_Lrad) became 0 * exp(-0/0) = NaN when that memory held 0, so examples/ex4.in produced NaN temperatures and a NaN stress in its second steady-state iteration; `_wMax` was also unset.
- Evidence: `_Lrad` now defaults to 10 km (irrelevant when A0 = 0); later ex4 runs complete their viscosity iterations (`9745bb2`, `c366a1b`).

**A-10. A zero shear-zone width divided by zero.** Critical. `25e29c9`. `source/heatEquation.cpp`.
- Defect: `constructMapV()` built the Gaussian heat kernel whenever `wVals` was given, so `wVals = [0 0]` (examples/ex5.in) divided by zero and aborted on the NaN assert before the first step.
- Evidence: the kernel is built only when max(w) > 0, a mix of zero and positive widths stops with a message, and ex5 runs (`6b2315e`).

**A-11. Matrix preallocation arrays on the stack.** Critical (latent). `009d188`. `source/spmat.cpp`, `source/genFuncs.cpp`.
- Defect: `kronConvert()` declared variable-length arrays `d_nnz[m]` and `o_nnz[m]` on the stack, and 2 x 4 bytes x Ny*Nz rows on one rank exceeds the default 8 MB stack near Ny*Nz = 1e6 (BP1 at full resolution, 401 x 2401, has 962,801 rows).
- Evidence: threshold computed in the message, no overflow reproduced; with `2005eed` the code compiles without warnings, so `WERROR=1` works.

**A-12. The heat equation was used without being created.** Critical. `aece9f9`. `source/strikeSlip_linearElastic_qd.cpp`, `source/strikeSlip_linearElastic_qd_fd.cpp`.
- Defect: the quasi-dynamic class created the heat equation only for `thermalCoupling != no` but used it for `evolveTemperature = 1` and `computeSSHeatEq = 1` (null dereference), and the qd_fd class did the opposite.
- Evidence: both now create it when any of the three asks for it; ex1 and ex2 bit-identical.

**A-13. Power-law qd_fd built the ice-stream forcing before the material existed.** Critical. `e77dcda`. `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: the forcing term was built in `allocateFields`, before `_material` existed, which dereferenced a null pointer for `forcingType = iceStream`.
- Evidence: no run reported; it is now built in the constructor after the material.

**A-14. The ice-stream context was written to a destroyed viewer.** Critical. `2febba7`. `source/strikeSlip_powerLaw_qd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: both power-law classes wrote the ice-stream forcing term to an ASCII viewer destroyed a few lines earlier, so `writeContext` failed for `forcingType = iceStream` and the run never integrated.
- Evidence: no ice-stream run reported; ex1 and ex2 bit-identical.

**A-15. Pressure fields missing from the integrand at phase changes and in steady-state iterations.** Critical. `04298a8`. `source/strikeSlip_linearElastic_qd_fd.cpp`, `source/strikeSlip_powerLaw_qd.cpp`.
- Defect: `prepare_qd2fd` copied into `_varFD["pressure"]`, which `initiateIntegrand_fd` never created (that code was commented out), so the first quasi-dynamic to dynamic switch passed a NULL Vec and pressure and permeability were not carried through dynamic phases; `StrikeSlip_PowerLaw_qd::solveSStau` also left the pressure fields out of its integrand.
- Evidence: coupled pore pressure on the ex2 geometry with kL_p = 1 mm through the first event: implicit quasidynamic_and_dynamic 17544 dynamic steps, exit 0, and explicit quasidynamic_and_dynamic exit 0.

**A-16. Implicit pore pressure aborted on constant grids.** Critical. `ca11678`. `source/pressureEq.cpp`.
- Defect: `PressureEq::be` asked the operator for its coordinate transform on every grid and `SbpOps_m_constGrid::getCoordTrans` asserts, so implicit pore pressure aborted on the first step of any `constantGridSpacing` run.
- Evidence: ex2 geometry on a constant grid with implicit coupled pressure aborts at step 0 with the old code and runs 50 steps with the fix.

**A-17. updateBCMats destroyed the wrong boundary matrices.** Medium (latent). `ca11678`. `source/sbpOps_m_constGrid.cpp`.
- Defect: the top Dirichlet branch destroyed the right-boundary Dirichlet matrices (`_AR_D`, `_rhsR_D`) instead of the unused top Neumann ones, leaving `_AR`/`_rhsR` dangling; reaching it needs a 2D constant-grid operator whose coefficient is updated with Dirichlet right and top boundaries.
- Evidence: a build with only this line reverted gives identical output for the pore-pressure operator, also under `-malloc_debug`; the variable-grid version was already correct.

**A-18. LinearElastic::setupKSP created one solver and configured another.** Medium. `ca11678`. `source/linearElastic.cpp`.
- Defect: it created `_ksp` but configured its `ksp` argument, so a repeated call leaked a solver and its factorization, and a different `ksp` argument was never created.
- Evidence: it now destroys and recreates the context it is given; ex1 and ex2 unaffected.

**A-19. The slip-rate ceiling was applied while the array was checked out.** Medium. `8954864`. `source/fault.cpp`.
- Defect: with `limitSlipVel`, `Fault_qd::computeVel` called `imposeSlipVelCeiling`, which takes the slip-rate array again, before restoring it.
- Evidence: no `limitSlipVel` run reported; the call now follows the restore, and ex1 and ex2 are bit-identical.

**A-20. Power-law qd_fd copied a temperature entry that does not exist.** Medium. `28b40f9`. `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: at each switch to the dynamic phase it copied `"Temp"` whenever thermal coupling was on, but that entry is integrated, and present in the maps, only with `evolveTemperature = 1`.
- Evidence: no run reported; the copy now depends on `evolveTemperature`.

**A-21. The fault Green's-function driver failed in parallel.** Critical. `89a6b0b`. `source/main.cpp`.
- Defect: `computeGreensFunction_fault` looped over each rank's own boundary entries while calling collective solves inside ("Object is in wrong state" on 2 ranks), never assembled `bcL`, restored the wrong array and inserted all Ny rows from the local part of `surfDisp`.
- Evidence: at Ny = Nz = 21 the serial G is identical to the old code and 2 ranks agree to 1.5e-16; G(surface at fault, top node) = 0.998.

**A-22. PETSc objects leaked.** Medium. `2febba7`, `9866e9b`, `f8fa0e4`, `4942bd8`, `aece9f9`, `e77dcda`, `89a6b0b`, `8916df7`. Several files.
- Defect: thirteen places created a binary viewer and then overwrote its handle with `PetscViewerHDF5Open`; several Vecs were never destroyed (`_sdev`, PowerLaw's `_wetDist`, the heat equation's `maxdT`, qd_fd's `SS_index`, temperature Vecs in three constructors, `viscSource` in the off-fault Green's-function driver); PressureEq allocated five Vecs twice.
- Evidence: after `8916df7`, PETSc's `-objects_dump` lists no SCycle objects at exit for ex2, ex5, the 2D quasidynamic_and_dynamic case, implicit pore pressure and the fully dynamic mode; the original code left 44 unfreed blocks after 20 steps of ex2.

### 3.2 Physics and numerics (wrong results)

Section 4 summarizes how much each of these changes results.

**A-23. Wrong Jacobian in the quasi-dynamic velocity solve.** Medium. `da28ff5`. `source/fault.cpp`.
- Defect: for F(V) = A asinh(B V) - (tauQS - eta V) the code used A V / sqrt(1 + B^2 V^2) + eta instead of A B / sqrt(1 + B^2 V^2) + eta, many orders of magnitude too small at interseismic slip rates, so BracketedNewton rejected most Newton steps and fell back to bisection.
- Evidence: measured during the audit, not recorded in the commit message: the slip-velocity solve time in ex2 fell from 4.88 s to 0.67 s. Together with A-24, final slip changes by about 1e-6 relative and event onsets agree to 5 or 6 digits; the converged velocity is unchanged up to rootTol.

**A-24. The root finder never accepted a bracket collapsed to round-off.** Critical. `74fdbb9`. `source/rootFinder.cpp`.
- Defect: BracketedNewton stopped only when |f| < rootTol, an absolute tolerance on the stress residual (default 1e-12 MPa); when the residual's round-off floor (about eps x |stress|) is larger, it ran to maxNumIts (1e4) and aborted on an assert.
- Evidence: seen on examples/ex4.in with the old and current code (tauQS about -3.8e4 MPa, bracket at adjacent doubles, |f| = 2.3e-12); the iterate is now accepted once the bracket is within 4 ulps, which changes V by at most a few ulps where the old code converged.

**A-25. Bisect::setBounds lost the bracket when the bounds were reversed.** Low. `9be7f3d`. `source/rootFinder.cpp`.
- Defect: with left > right, both bounds were set to right; `Bisect` appears only in commented-out code, so no run was affected.
- Evidence: code fix only; results unchanged.

**A-26. Absorbing boundaries were disabled in dynamic phases.** High. `4777e9a`. `source/strikeSlip_linearElastic_qd_fd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: a 2018 edit (`896a5e3`) dropped the body of an `if` in `computePenaltyVectors()`, so the `if` guarded the loop increment and the damping vector `_ay` stayed essentially zero; in every fully dynamic phase of quasidynamic_and_dynamic runs (elastic and power law), outgoing waves reflected off the right and bottom boundaries.
- Evidence: no run reported; the code now matches `StrikeSlip_LinearElastic_fd::computePenaltyVectors`.

**A-27. The transient heat solve dropped the heat source on constant grids.** High. `e665cde`. `source/heatEquation.cpp`.
- Defect: `be_transient` filled the temporary holding J*Q only for `variableGridSpacing`, so with `constantGridSpacing` frictional heating in a finite-width shear zone and viscous shear heating were dropped (or garbage was added).
- Evidence: spot test during the audit, not recorded in the commit message: ex2 with the heat equation on a constant grid keeps T at 300 K with the old code and reaches a maximum of 308.4 K with the fix, which now matches the explicit `d_dt` path (J is the identity).

**A-28. The transient heat solve reused the first step's factorization.** High. `b8f197e`. `source/heatEquation.cpp`.
- Defect: with `linSolver_heateq = MUMPSLU` or `MUMPSCHOLESKY`, `KSPSetReusePreconditioner` kept the factorization of I - dt*D2 from the first step, and with `KSPPREONLY` nothing corrects it, so every later temperature update used the first step's dt.
- Evidence: no run reported; reuse is now disabled, as in the CG and AMG branches.

**A-29. The fault heat flux had the wrong sign in two heat solvers.** High. `24f6732`. `source/heatEquation.cpp`.
- Defect: with zero shear-zone width, `be_steadyState` (`heatEquationType = steadyState`) and the explicit `d_dt` (fully dynamic phases) applied the frictional heat flux with the wrong sign and cooled the fault; a July 2022 sign change had fixed `be_transient` but broken these two.
- Evidence: 1D test (uniform 300 K, k = 1.89 W/m/K): steady-state guess 565.5 K at the fault (565 K predicted), `be_transient` unchanged (328.9 K during events), `steadyState` now heats (+5e11 K at coseismic rates instead of -5e11 K). The explicit `d_dt` path was checked after the commit: a 2D quasidynamic_and_dynamic run with the fault heat flux (w = 0) through the first event raises the maximum temperature at every one of the 1755 saved dynamic steps, from 0.33 K to 2.22 K, consistent with the frictional heat spread over the first grid cell.

**A-30. Two-point temperature profiles lost the geotherm.** High. `06bada2`. `source/heatEquation.cpp`.
- Defect: with two-entry `TVals`/`TDepths`, the base of the lithosphere was set to the surface, so the whole domain took the top temperature.
- Evidence: ex2 with TVals = [283 900] over 0 to 30 km gives 283 K everywhere with the old code and 283, 591.5, 900 K at top, middle and bottom with the fix.

**A-31. Viscous dissipation used wrong factors.** High. `574a54b`. `source/heatEquation.cpp`, `source/grainSizeEvolution.cpp`.
- Defect: since 2019 the heat source used sqrt(2) x sdev x dgdev (41% more heat than the work done) and the paleowattmeter grain-size reduction used 0.5 x sdev x dgdev_disl (half the work), while with the code's definitions the dissipation is sdev x dgdev, so energy was not conserved.
- Evidence: no run numbers; marked PHYSICS CHANGE: power-law runs with viscous shear heating or grain-size evolution change, linear-elastic runs do not, and reverting this commit alone restores the old factors.

**A-32. The shear modulus profile was interpolated across the fault-normal direction.** High. `9608cec`. `source/linearElastic.cpp`.
- Defect: `LinearElastic` built mu from `muVals`/`muDepths` as a function of y, while PowerLaw, Fault (radiation damping) and the manual use depth, so a depth-varying mu gave a bulk stiffness that varied with distance from the fault.
- Evidence: no run with non-uniform mu reported; uniform mu (all shipped examples) is unaffected.

**A-33. Fourth-order operators depended on the number of MPI ranks.** High. `b65b5e0`. `source/sbpOps_m_constGrid.cpp`, `source/sbpOps_m_varGrid.cpp`, `source/spmat.cpp`.
- Defect: the midpoint coefficient B3 of the 4th-order R term was averaged using locally owned values only, so the last local row (two rows for R_z) on every rank but the last kept half the coefficient.
- Evidence: ex2 on 2 vs 1 ranks differed by up to 3e-6 (4e-8 relative) at the partition boundary; now the difference is 3.5e-10 (4e-12 relative) spread over the domain; ex1 bit-identical.

**A-34. R_y averaged the wrong neighbours.** High. `b65b5e0`. `source/sbpOps_m_constGrid.cpp`, `source/sbpOps_m_varGrid.cpp`.
- Defect: R_y averaged z-neighbours (Ii, Ii+1) instead of y-neighbours (Ii, Ii+Nz).
- Evidence: ex2 serial changes by up to 4e-8 relative; MMS at order 4 (Ny = Nz = 21, 41, 81, variable and constant grids) keeps the L2 error in u within 0.05% and the rate at 3.5; order 2 is unaffected.

**A-35. Compatible 4th-order boundary rows had coefficient errors.** High. Partly fixed. `e892279`. `source/spmat.cpp`, `source/domain.cpp`.
- Defect: in `sbp_Spmat` (order 4, compatible) the right-boundary derivative row lacked its 1/h factor, and the modified D1 boundary row set column N-4 twice and never N-3.
- Evidence: after the fix the MMS displacement errors are still about 0.2 at N = 21, 41 and 81 (fullyCompatible, order 2: 6.9e-5, 1.4e-5, 3.1e-6), so the option is now refused (O-02).

**A-36. getMus returned mu for the metric-weighted moduli.** Low (latent). `dd9de2f`. `source/sbpOps_m_varGrid.cpp`.
- Defect: `SbpOps_m_varGrid::getMus` returned `_mu` for all three outputs instead of mu x q_y and mu x r_z; the only caller (`PowerLaw::initializeMomBalMats`) uses mu alone.
- Evidence: ex1 and ex2 bit-identical.

**A-37. The documented cohesion had no effect.** High. `6b3dad7`, `aa53679`. `source/fault.cpp`.
- Defect: `cohesionVals`/`cohesionDepths` (MPa) were read, written to output and passed to the velocity solve, but never entered the fault strength, in the quasi-dynamic solve or in the fully dynamic one.
- Evidence: zero cohesion leaves ex1 and ex2 bit-identical; with 5 MPa in ex1 the fault starts at steady state (tau = 38.45 = 33.45 + 5 MPa) and the first event occurs at the same time. In the dynamic phase (`aa53679`, 2D quasidynamic_and_dynamic through the first event), zero cohesion is bit-identical to the previous build; with 5 MPa the old code reached a peak V of 28.93 m/s and 20.96 m of slip with the dynamic strength down to 4.00 MPa, below the cohesion itself, while the fix gives 25.13 m/s and 19.72 m, close to the 25.09 m/s and 19.70 m without cohesion, as a uniform cohesion that only shifts the stress level should.

**A-38. Creeping nodes slipped freely in dynamic phases.** High. `8954864`. `source/fault.cpp`.
- Defect: `Fault_fd::d_dt` updated every node faster than 1e-14 m/s from the wave-equation balance, so nodes marked as creeping (`lockedVals < -0.5`) were not held at vL during fully dynamic phases.
- Evidence: 2D quasidynamic_and_dynamic run with creep below 20 km through the first event (14958 dynamic steps): old code up to V = 1.67e-9 m/s at those depths, new code V = vL = 1e-9 exactly; ex1 and ex2 bit-identical.

**A-39. The initial effective normal stress ignored the pore pressure.** High. `61bbac0`. `source/strikeSlip_linearElastic_qd.cpp`.
- Defect: the coupling statement was an `else if` behind `if (_hydraulicCoupling != "no")` and could never run, so with `hydraulicCoupling = coupled` the fault and the steady-state guess started from the uncoupled normal stress.
- Evidence: no run reported; it is now an independent `if`, as in the qd_fd and power-law classes.

**A-40. Power-law runs did not evolve slip-dependent permeability.** High. `325945f`. `source/strikeSlip_powerLaw_qd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: `d_dt` computed the pressure rates before the fault slip rate, so fault valving always saw V = 0, and with implicit pressure the permeability was never evolved.
- Evidence: no run reported; the call order now matches the linear-elastic classes.

**A-41. The second dislocation mechanism was only partly wired in.** Medium. `9745bb2`. `source/powerLaw.cpp`.
- Defect: with `wDislCreep2 = yes`, the dislocation strain rate that drives grain-size reduction and the heating partition, and the initial effective-viscosity guess, ignored the second mechanism; `wLinearMaxwell = yes` left it switched on.
- Evidence: no run with two mechanisms reported.

**A-42. The initial viscosity guess was off by a units factor.** Low. `9745bb2`. `source/dislocationCreep.cpp`, `source/diffusionCreep.cpp`.
- Defect: `guessInvEffVisc` omitted the factor 1e3 that `computeInvEffVisc` uses for strain rates in 1e-3/s, so the guessed stress was 1e3^(1/n) too large (10 times the viscosity for n = 3); it is only an initial guess.
- Evidence: ex4 at 101 x 121 still converges in 44 steady-state viscosity iterations, to the same tolerance.

**A-43. Dissolution-precipitation creep returned a strain rate as an inverse viscosity.** High. `b42c53d`. `source/dissolutionPrecipitationCreep.cpp`.
- Defect: `computeInvEffVisc` omitted the division by stress in the formula documented in the class header, so the mechanism was wrong by a factor of the stress in MPa; the gas constant was also stored in single precision.
- Evidence: it now uses expm1(k s)/s, which tends to k as s -> 0; ex4 at 101 x 121 with the mechanism on runs without NaN or asserts.

**A-44. Dissolution-precipitation creep did nothing by default.** Medium. `b42c53d`. `source/powerLaw.cpp`.
- Defect: the wet fraction defaulted to 0 unless a `pl_wetDist` input file existed, so `wDissPrecCreep = yes` silently had no effect.
- Evidence: it now defaults to 1 (wet everywhere) when no file is given; same ex4 check as A-43.

**A-45. The steady-state viscosity iteration diverged without warning.** High. Partly fixed. `c366a1b`. `source/strikeSlip_powerLaw_qd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: the fixed-point iteration for the steady-state effective viscosity used a fixed damping (`fss_EffVisc`, default 0.2) and stopped silently at its iteration limit; ex4 at 301 x 351 reached a relative change of 0.34 and then grew.
- Evidence: the damping now halves (down to 1% of `fss_EffVisc`) when the change grows, and a warning is printed at the limit; ex4 at full resolution converges in 43 iterations, a 101 x 121 variant with dissolution-precipitation creep converges in 72 after one halving (not in 50 before), and the converged solution is unchanged; outer iteration 1 of ex4 still hits the iteration limit (O-05).

**A-46. Top-boundary loading used depths as horizontal coordinates.** High. `e77dcda`. `source/strikeSlip_powerLaw_qd.cpp`.
- Defect: `updateBCT_atan_u` read `Domain::_y0` (length Nz, the depths along the fault) as the y-coordinates of the top boundary (length Ny), giving wrong boundary data and out-of-range reads when Ny > Nz; `updateBCT_atan_v` restored its array on the wrong Vec.
- Evidence: no run reported; y is now scattered onto the top boundary, as `updateBCT_atan_v` does.

**A-47. The implicit pressure Picard loop ignored its tolerance.** Medium. `4942bd8`. `source/pressureEq.cpp`.
- Defect: the backward-Euler Picard loop computed its relative change but never stopped on `minBeDifference`, so it always ran `maxBeIteration` iterations (default 1), and it could divide by a zero norm.
- Evidence: no separate run reported.

### 3.3 Time integration, checkpoints and I/O

**A-48. The grain-size limit could override the Maxwell time-step limit.** Medium. `00b3b33`. `source/strikeSlip_powerLaw_qd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: with transient grain-size evolution, `timeMonitor()` replaced the cap of 0.9 Maxwell times by min(maxDeltaT, 0.9 grain-size relaxation times), so steps could exceed the Maxwell time.
- Evidence: no run reported; the minimum of both limits is used.

**A-49. The next step size could be unset or stale.** Critical. `39f8663`. `source/odeSolver.cpp`, `source/odeSolverImex.cpp`.
- Defect: when the error estimate was exactly zero the RK integrators did not update `_newDeltaT`; the IMEX integrators never initialized it, so a 1D constant-grid thermal run stopped at t = 1e-3 s and looped forever with dt = 0 (depending on memory layout), and the explicit ones stalled at the initial step size.
- Evidence: `_newDeltaT` is initialized and a zero error grows the step by the factor-5 cap (bounded by maxDeltaT); ex1 and ex2 never see a zero error and are bit-identical.

**A-50. Rejected steps were accepted with the wrong step size.** High. `fbe6ee0`. `source/odeSolver.cpp`, `source/odeSolverImex.cpp`.
- Defect: in all four adaptive RK integrators, when a shrunk step hit minDeltaT (or the 100th attempt) the loop accepted the state computed with the larger step while time advanced by the smaller one; ex2 hits the minimum step hundreds of times per event.
- Evidence: ex2 keeps its event onset (54.55 yr), its peak slip rate goes from 6.351 to 6.346 m/s and its final slip changes by 6e-7 relative; ex1 is unchanged.

**A-51. Checkpoints stored a never-set error estimate.** Low. `fbe6ee0`. `source/odeSolver.cpp`, `source/odeSolverImex.cpp`.
- Defect: locals in `integrate()` and `computeError()` shadowed the `_totErr` member, so checkpoints wrote an unset value.
- Evidence: no separate test; the member now holds the last step's error estimate.

**A-52. The relative error norms were broken.** Critical. `d823c7b`. `source/odeSolver.cpp`, `source/odeSolverImex.cpp`, `source/genFuncs.cpp`.
- Defect: with `normType = max_relative`, RK32 compared against a Vec it never fills (a PETSc error) and RK32 and RK43 divided by the signed solution (negative or infinite errors); `L2_relative` divided by a zero norm for an identically zero variable (NaN or Inf time steps).
- Evidence: a new helper computes max |a - b| / |b|, and zero norms fall back to the absolute error; ex1 runs with L2_relative and max_relative for RK32 and RK43, and ex1 and ex2 are bit-identical.

**A-53. The PID controller produced NaN in integrators created mid-run.** Critical. `d823c7b`. `source/odeSolver.cpp`, `source/odeSolverImex.cpp`.
- Defect: PID, the default `timeControlType`, used the error history whenever stepCount >= 4, but an integrator created mid-run (each phase of a quasidynamic_and_dynamic run) has no history, which gave 0 x inf = NaN in the step-size ratio.
- Evidence: proportional control is used until the history exists; no qd_fd run in this commit, but ex2 later runs as quasidynamic_and_dynamic with the default PID (`b78b1cb`).

**A-54. FEuler ignored stop requests and could not checkpoint in parallel.** Medium. `d823c7b`. `source/odeSolver.cpp`.
- Defect: FEuler ignored `stopIntegration` from `timeMonitor`, and its checkpoint Vec had local size 1 on every rank, which fails on more than one rank.
- Evidence: ex1 runs with FEuler.

**A-55. IMEX coupling lagged one step.** High. `185a2ac`. `source/odeSolverImex.cpp` and the four quasi-dynamic `source/strikeSlip_*.cpp` classes.
- Defect: in RK32_WBE and RK43_WBE every stage of step n+1 used the old implicit state (temperature T^n, effective normal stress from p^n), so thermal and hydraulic feedback lagged one step and a restart diverged from the uninterrupted run.
- Evidence: IMEX restarts now continue bit for bit (implicit pore pressure on ex2; ex5 with RK43_WBE and RK32_WBE, 300 + 300 vs 600 steps); the ex5 peak fault temperature goes from 289.120 to 289.141 K; ex1 and ex2 bit-identical.

**A-56. The quasi-dynamic to dynamic hand-off lost the particle velocity.** High. `b78b1cb`. `source/strikeSlip_linearElastic_qd_fd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: with the explicit integrators, the single quasi-dynamic step in `prepare_qd2fd` took no step (they were built with maxNumSteps = 1, a limit `setInitialStepCount` does not raise for them), and the wave solver started with the n-1 level equal to level n, so each dynamic phase began without the quasi-dynamic particle velocity.
- Evidence: ex2 run as quasidynamic_and_dynamic: first-event peak slip rate 27.11 -> 26.57 m/s at the same step.

**A-57. Slip-dependent permeability went unstable during events.** Critical. Partly fixed. `04298a8`, `35ba52f`. `source/pressureEq.cpp` and the four quasi-dynamic `source/strikeSlip_*.cpp` classes.
- Defect: dk/dt = -|V|/L (k - kmax) - (k - kmin)/T was always integrated explicitly and never entered step-size control; with kL_p = 1 mm on the ex2 geometry, k oscillated to 1e297 in magnitude during the first event, the next implicit pressure solve failed inside HYPRE, and the run continued until it segfaulted at exit (exit 59).
- Evidence: permeability is now relaxed exactly over each step (implicit mode and dynamic phases) and explicit hydraulic fields join the step-size control; with kL_p = 1 mm through the first event at 54.55 yr the implicit quasidynamic run takes 1709 steps (exit 0), k stays within [1e-19, 1e-17], a restart from step 400 matches the uninterrupted run at step 450 bit for bit, and the explicit quasidynamic run still diverges (O-01). Follow-up `35ba52f`: `04298a8` had replaced the default of an empty `timeIntInds` (all explicit variables) with the hydraulic fields alone, dropping psi and slip from step-size control; an empty list is now left alone and the integrator again reports [permeability pressure psi slip] with scales [1 1 1 1] (no example was affected).

**A-58. RK32 restarts read the wrong HDF5 group.** Critical. `ab5973c`. `source/odeSolver.cpp`.
- Defect: `RK32::writeCheckpoint` stored its state under /odeSolver/odeSolver_chkpt_data, but `RK32::loadCheckpoint` looked in /time1D, so every restart with `timeIntegrator = RK32` failed.
- Evidence: spot test during the audit, not recorded in the commit message: ex1 with RK32 could not restart from its checkpoint with the old code; with the fix a restarted run is bit-identical to an uninterrupted one.

**A-59. Checkpoint and steady-state files stayed open after loading.** Medium. `ab5973c`, `9866e9b`. `source/odeSolver.cpp`, `source/odeSolverImex.cpp`, `source/odeSolver_WaveEq.cpp`, `source/powerLaw.cpp`.
- Defect: no integrator's `loadCheckpoint` destroyed its viewer on checkpoint.h5, which stayed open while the problem class reopened the file to write the next checkpoint; `OdeSolver_WaveEq` pushed an HDF5 group without popping it; `PowerLaw::loadCheckpointSS` left data_steadyState.h5 open for the whole run.
- Evidence: no run reported; restarts are exercised in `3949db3` and `04298a8`.

**A-60. qd_fd regime flags were checkpointed through the wrong type.** Critical. `da39e0e`. `source/strikeSlip_linearElastic_qd_fd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: the bool members `_inDynamic` and `_allowed` were written and read as PETSC_INT attributes through their addresses, so writing stored padding bytes in the file and reading wrote 4 bytes into a 1-byte member, overwriting its neighbour.
- Evidence: no run reported; PetscInt temporaries are used and the file format is unchanged.

**A-61. Restart flags without a checkpoint file aborted the run.** Critical. `14a1331`. `source/domain.cpp`.
- Defect: with `restartFromChkpt = 1` and `restartFromChkptSS = 1` but no checkpoint.h5, the grid was never built and every component tried to load steady-state files that did not exist ("Unable to open file"); correcting the misspelled restart key in ex1 and ex2 would have triggered this.
- Evidence: spot test during the audit: with both restart flags set and no checkpoint file the old code aborted on "Unable to open file"; the fixed code falls through to the steady-state restart and then to a fresh start, and a restart prints the file it continues from.

**A-62. Pore pressure was never written or checkpointed.** Critical. `4942bd8`. `source/pressureEq.cpp`.
- Defect: the output held only material properties (not p, its rate or k), checkpoints wrote the material properties twice and none of the evolving state, and the constructor asserted on any restart; `_kmin2_p` was also named "meanP_p".
- Evidence: ex2 with `hydraulicCoupling = coupled` for 50 steps writes the pressure to data_1D.h5, and with explicit pressure a restart at step 25 reproduces the uninterrupted run bit for bit (implicit mode after `185a2ac`).

**A-63. Steady-state restarts overwrote the first dislocation mechanism.** Medium. `9866e9b`. `source/powerLaw.cpp`.
- Defect: `loadCheckpointSS()` loaded the second mechanism's 1/eta (group /momBal/dislocationCreep2) into `_disl`, so steady-state restarts with two mechanisms lost the first.
- Evidence: no run reported.

**A-64. sxz was never written; sdev was misnamed.** Critical. `f8fa0e4`. `source/linearElastic.cpp`.
- Defect: with `momBal_computeSxz = 1`, `writeStep2D`, `writeCheckpoint` and both checkpoint loaders handled `_sxy` a second time instead of `_sxz`, so sxz was never written and a restart loaded sxy twice; with `momBal_computeSdev = 1` the code named `_sxz` instead of `_sdev`, an error when sxz is not allocated.
- Evidence: no run reported.

**A-65. A job killed during a checkpoint write could lose the only restart file.** Critical (latent). `3949db3`. `source/genFuncs.cpp` and the four quasi-dynamic `source/strikeSlip_*.cpp` classes.
- Defect: checkpoint.h5 was overwritten in place through a viewer open for the whole run, and data_1D.h5 and data_2D.h5 were never flushed (`PetscViewerFlush` is a no-op for HDF5 viewers).
- Evidence: checkpoints now go to checkpoint.h5.tmp, which is closed and renamed after the output files are flushed; ex1 run to step 1000 and restarted to 1500 equals an uninterrupted 1500-step run bit for bit, and ex2 killed with SIGKILL at step 631 restarts from step 620; the old build also survived a kill at step 654, so the failure was not reproduced.

**A-66. The off-fault Green's function stopped after 10 columns.** Critical. `89a6b0b`. `source/main.cpp`.
- Defect: the loop in `computeGreensFunction_offFault` had been cut to startIi+10 while debugging, so a run produced 10 of its 2 x Ny x Nz columns.
- Evidence: Ny = Nz = 11 now gives all 242 columns and a rerun on the finished file does no more work; the unit source is now set by its owning rank, the file is flushed after each column, and a file without columns restarts at column 0.

**A-67. The 2D output stride was taken from the 1D stride.** Medium. `28b40f9`, `2febba7`. `source/strikeSlip_linearElastic_qd_fd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`, `source/strikeSlip_linearElastic_fd.cpp`, `source/strikeSlip_powerLaw_qd.cpp`.
- Defect: both qd_fd classes set the 2D stride of the fully dynamic phase from the 1D stride when switching phases, and in three classes the "stride2D" line of the context file printed the 1D stride.
- Evidence: no run reported; ex1 and ex2 bit-identical.

**A-68. The power-law qd_fd class ignored initTime.** Medium. `28b40f9`. `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: only `_initTime` was set from `initTime`, but phases start at `_currTime`, so a fresh run always started at t = 0.
- Evidence: no run reported; both are now set, as in the linear-elastic qd_fd class.

### 3.4 Input handling and error reporting

**A-69. Failed runs exited with status 0.** Medium. `e9cdc49`. `source/main.cpp`.
- Defect: `main()` ignored the return codes of every driver, so failed runs looked successful to scripts.
- Evidence: the code is propagated and the exit status is 1 on error (checked in the code).

**A-70. The list parser dropped or invented entries.** High. `7232307`. `source/genFuncs.cpp`.
- Defect: a one-element numeric list such as `muVals = [30]` came back empty, and repeated or leading or trailing spaces inside the brackets became extra entries (0 for numbers).
- Evidence: one whitespace tokenizer now serves all three overloads; well-formed lists with single spaces parse exactly as before.

**A-71. Depth profiles were not validated.** Critical. `aece9f9`. `source/genFuncs.cpp`, `source/fault.cpp`.
- Defect: `setVec`, which builds every `<name>Vals`/`<name>Depths` profile, read element 0 of an empty list, read past the end of lists of different lengths, looped about 2^64 times on an empty depth list, divided by zero on a repeated last depth and silently built garbage from non-monotonic depths (as from ex4's duplicated sN lists); `VwType = function_of_Tw` called it on an empty `VwVals`.
- Evidence: it now stops with a message listing the values and depths, and a repeated depth gives a step in the profile; ex1 and ex2 bit-identical.

**A-72. A missing input file was not reported.** Medium. `aece9f9`. `source/domain.cpp`.
- Defect: Domain kept every setting at its default and the run failed later on an unrelated assert.
- Evidence: the error now names the file; no run numbers.

**A-73. The NaN and Inf checks in dissolution-precipitation creep could not fire.** Medium. `2005eed`. `source/dissolutionPrecipitationCreep.cpp`.
- Defect: `assert(~PetscIsNanReal(expVal))` applies bitwise NOT to a bool, which is always nonzero, so an overflowing exponential silently gave an infinite (or, with zero wet fraction, NaN) inverse viscosity.
- Evidence: logical NOT is used, and the clang `-Wbool-operation` warning is gone.

**A-74. A rigidFault left boundary was always rejected in fully dynamic runs.** Critical. `28b40f9`. `source/strikeSlip_linearElastic_fd.cpp`.
- Defect: `checkInput` tested `_bcRType` where it meant `_bcLType`, so `momBal_bcL_fd = rigidFault` always failed the assert; `outGoingCharacteristics`, which the manual lists for bcL, is now accepted too.
- Evidence: no run reported.

**A-75. The regime-switch fallback compared limit_qd with the wrong trigger.** Medium. `b78b1cb`. `source/strikeSlip_linearElastic_qd_fd.cpp`, `source/strikeSlip_powerLaw_qd_fd.cpp`.
- Defect: `checkInput` reset `limit_qd` when it exceeded `trigger_fd2qd`, but `limit_qd` re-arms the quasi-dynamic to dynamic switch and must be compared with `trigger_qd2fd`; the power-law class had neither fallback, and the step log printed "maxReq", which has been 0 since the switch criterion became the slip rate (`74a132f`).
- Evidence: no run reported.

**A-76. Root-finder iteration counts were never recorded.** Low. `9be7f3d`. `source/rootFinder.cpp`, `source/fault.cpp`.
- Defect: every `findRoot` counted in a local variable, so `getNumIts()` returned 0 and the fault's root-iteration statistics stayed empty; builds with `-DVERBOSE=4` or higher failed in RegulaFalsi.
- Evidence: the fault's counter is now 64-bit and printed by `Fault::view` (ex1: 48299 iterations in 2201 steps); results unchanged.

**A-77. The initial pore pressure pVals was ignored.** High. `4942bd8`. `source/pressureEq.cpp`.
- Defect: the required key `pVals` (initial pore pressure, per the manual) was never used.
- Evidence: no separate run reported; it is used when `guessSteadyStateICs = 0`.

**A-78. The mediators did not parse hydraulicTimeIntType; implicit pressure could freeze.** High. `04298a8`. The four quasi-dynamic `source/strikeSlip_*.cpp` classes.
- Defect: the mediators never set their copy of `hydraulicTimeIntType`, and implicit pressure paired with an explicit-only integrator (RK32 or RK43, which never update the implicit variables) left the pressure frozen without a message.
- Evidence: the key is parsed and that pairing now stops with a message (checked in the code).

**A-79. The mediators ignored errors from the pressure and grain-size rates.** Critical. `04298a8`. The four quasi-dynamic `source/strikeSlip_*.cpp` classes.
- Defect: the return codes of the pressure and grain-size `d_dt` calls were not checked, so the failed pressure solve of A-57 let the run continue on bad state until it segfaulted at exit.
- Evidence: the codes are now checked; see A-57 for the runs.

### 3.5 Build, repository and tooling

**A-80. The Makefile named a nonexistent linker and hid link failures.** Low. `44a81b1`. `source/Makefile`.
- Defect: it linked with `openmpicc` (which PETSc's variables overrode anyway), ignored link errors through a leading '-', and deleted main.o after every link.
- Evidence: it now links with PETSc's `CXXLINKER` and `PETSC_LIB`; every build in this audit used it with PETSc 3.26.0.

**A-81. Stale header dependencies and mandatory -Werror.** Low. `44a81b1`. `source/Makefile`.
- Defect: the hand-written dependency list named files that do not exist (`mainLinearElastic.cpp`, `sbpOps_sc.hpp`) and missed headers, -Werror was always on although new compilers warn on this code, and an FDP target pointed at sources in `source/sfsbp/`.
- Evidence: dependencies come from `-MMD -MP` and -Werror is opt-in (`WERROR=1`); after `009d188` and `2005eed` the code compiles without warnings.

**A-82. Binaries and OS files were tracked.** Low. `8c44ddc`. `.gitignore`.
- Defect: two Linux ELF test executables, generated test data and .DS_Store files were in git, and `.gitignore` ignored itself.
- Evidence: untracked (the files stay on disk) and covered by `.gitignore`.

**A-83. There was no automated regression check.** Low. `26139e7`, `9c0cc61`. `tools/regress.sh`.
- Defect: nothing compared output across code changes; the first version of the script also passed runs whose datasets had a different shape, because h5diff exits 0 when objects are not comparable.
- Evidence: used after every later code commit; non-comparable datasets now count as differences.

**A-84. The manufactured-solution test restarted from the previous resolution.** Medium. `14a1331`, `e0b9041`. `source/domain.cpp`, `tools/mms.in`.
- Defect: the Domain constructor used by MMS runs defaulted `restartFromChkptSS` to 1 and append mode, so each resolution tried to restart from the previous resolution's files; the repository had no input file for the test.
- Evidence: `tools/mms.in` records the expected rates: order 4 about 3.5 for u and 2.5 for sigma_xy (errors 4.3e-9 and 1.2e-4 at Ny = 81), order 2 about 2.2 for u.

**A-85. Misspelled or retired input keys were ignored without a warning.** Low. `df64c96`. `tools/checkkeys.py`.
- Defect: SCycle silently ignores keys it does not read (for example `retartFromChkptSS` in ex1 and ex2); the code still does, but `tools/checkkeys.py` now lists them.
- Evidence: all examples, `SEAS_benchmarks/BP1/BP1.in` and `tools/mms.in` pass; for its limit see O-16.

**A-86. Unreachable debug drivers with hard-coded paths.** Low. `89a6b0b`. `source/main.cpp`.
- Defect: `runTests`, `initiateFields`, `runFirstStep`, `runSecondStep`, `testHDF5` and `computeGreensFunction_test` were referenced only in comments and contained hard-coded /Users/kallison paths.
- Evidence: removed; main.cpp went from 624 to 280 lines.

**A-87. Misleading comment on Domain::_y0 and _z0.** Low. `c844c59`. `source/domain.hpp`.
- Defect: the comment "q(y), r(z)" described neither: `_y0` is the size-Nz fault template and holds z, `_z0` the size-Ny surface template and holds y (the variable misused in A-46).
- Evidence: comment change only.

### 3.6 Examples, benchmarks and post-processing

**A-88. ex1 and ex2 started from zero shear stress.** High. `6b2315e`. `examples/ex1.in`, `examples/ex2.in`.
- Defect: since commit `089727a` (2023) the steady-state momentum balance behind `guessSteadyStateICs` runs only with `computeSSMomBal = 1`, which the examples lacked, and `restartFromChkptSS` was misspelled; ex2 had no earthquake in its 207-year window (peak slip rate 3e-11 m/s).
- Evidence: with the flag, ex2 starts at steady state and has an event at 54.5 yr.

**A-89. ex3 ran the wrong problem class.** High. `6b2315e`. `examples/ex3.in`.
- Defect: its header describes quasi-dynamic interseismic periods with fully dynamic ruptures, but `momentumBalanceType` was commented out, so it ran the quasi-dynamic class and ignored all switching and per-phase stride keys.
- Evidence: enabled, with `computeSSMomBal = 1` and a note that the trigger and limit values are slip velocities in m/s; ex3 runs (10 steps checked).

**A-90. ex4 and ex5 failed input validation.** Medium. `6b2315e`. `examples/ex4.in`, `examples/ex5.in`.
- Defect: ex4 used retired power-law key names, lacked `systemEvolutionType = steadyStateIts` and `computeSSHeatEq = 1`, and repeated its sN lists, which the parser concatenated into a non-monotonic profile; ex5 gave the flash-heating parameters as scalars instead of depth profiles and called 1732 K 900 C (it is about 1459 C).
- Evidence: ex5 runs; ex4 ran but its viscosity iteration did not converge at that commit (A-45, O-05).

**A-91. The examples' momentum-balance solver keys were ignored.** Medium. `6b2315e`. `examples/ex2.in`, `examples/ex4.in`, `examples/ex5.in`.
- Defect: these examples set `linSolver`, which only the pore-pressure solver reads; the momentum balance reads `linSolverSS` and `linSolverTrans`.
- Evidence: renamed; `tools/checkkeys.py` passes on all examples (`df64c96`).

**A-92. The Python loader did not read current output.** Medium. `f3b5829`. `examples/loadFuncs.py`, `examples/visualize_ex1.ipynb`, `examples/visualize_ex2.ipynb`.
- Defect: its comment stripping was overwritten by the next line, a mutable default dict made separate calls share results, the default key lists named datasets that do not exist (dt1D, tauP, slipVelocity) and skipped others, and the notebooks pointed at /Users/kallison/scycle/data.
- Evidence: with h5py 3.16 on ex1 and ex2 output all four files load with the documented shapes, `Ly` parses as 100.0, and separate calls return separate dicts.

**A-93. The MATLAB loaders called functions that were never committed.** Medium. `f3b5829`. `matlab/visualizePetsc/`, `examples/visualize_ex1.m`, `examples/visualize_ex2.m`.
- Defect: the HDF5 loaders called `h5read_reshape`, `recursiveGroupFinder`, `catField` and `findEvents_dyn`, which were not in the repository, `loadData1D_hdf5` appended time on each load, `loadStruct` mis-parsed lines, and the example scripts read PETSc binary files the code no longer writes.
- Evidence: not run, because the MATLAB license on the test machine has expired; reviewed against MATLAB's reversed HDF5 dimension order and every call site (O-09).

**A-94. SEAS benchmark BP1 could not be run.** Medium. `df64c96`. `SEAS_benchmarks/BP1/BP1.in`, `SEAS_benchmarks/BP1/createICs.py`.
- Defect: BP1.in pointed at a cluster path and used the retired keys `problemType` and `linSolver`, and its initial conditions came from createICs.m, which calls a function that was never committed (`constructCoord_constL1`) and writes to a laptop path.
- Evidence: the new createICs.py (numpy only) at `--dz 0.1` (251 x 751, 30 steps): all six files load, tau = 26.5461 MPa everywhere, V = 1e-9 m/s at t = 0, psi in the velocity-weakening zone rises at 1.9e-6 /s as the benchmark intends, peak memory 1.4 GB; the full 25 m grid is 401 x 2401.

## 4. Changes that alter results

Unless its message says otherwise, each code commit was followed by a bit-for-bit comparison of ex1 and ex2 output with `tools/regress.sh`. The commits below change results in at least one configuration. "Not reported" means the message gives no size. Commits that only turn a crash, a hang or a NaN into a result are not listed.

| Commit | Change | Runs affected | Reported size |
|---|---|---|---|
| `da28ff5` | Newton Jacobian (A-23) | every rate-and-state velocity solve | converged V unchanged up to rootTol |
| `74fdbb9` | Collapsed bracket accepted (A-24) | solves whose residual floor exceeds rootTol | at most a few ulps of V where the old code converged |
| `4777e9a` | Absorbing boundaries (A-26) | fully dynamic phases of quasidynamic_and_dynamic | not reported |
| `00b3b33` | Time-step cap (A-48) | power law with transient grain-size evolution | not reported |
| `e665cde` | Heat source on constant grids (A-27) | `constantGridSpacing` with transient heat and shear heating | not reported |
| `b8f197e` | No factorization reuse (A-28) | `linSolver_heateq` = MUMPSLU or MUMPSCHOLESKY | not reported |
| `7232307` | List parsing (A-70) | single-element or irregularly spaced lists | well-formed lists unchanged |
| `9608cec` | mu interpolated in depth (A-32) | linear elastic with non-uniform `muVals` | not reported; uniform mu unaffected |
| `325945f` | Pressure rates after the slip rate (A-40) | power law with slip-dependent permeability | not reported |
| `61bbac0` | Coupled initial sNEff (A-39) | linear-elastic quasi-dynamic, `hydraulicCoupling = coupled` | not reported |
| `b65b5e0` | 4th-order R coefficient (A-33, A-34) | order 4 | ex2 serial up to 4e-8 relative; 1 vs 2 ranks 3.5e-10; ex1 bit-identical; MMS error within 0.05% |
| `6b2315e` | Example inputs (A-88 to A-91) | ex1 to ex5 | ex2 now has an event at 54.5 yr (none in 207 yr before) |
| `24f6732` | Fault heat-flux sign (A-29) | `heatEquationType = steadyState`; explicit heat `d_dt` | -5e11 K -> +5e11 K at coseismic rates; dynamic phases now heat the fault (max dT 0.33 -> 2.22 K over one event) |
| `39f8663` | Step after a zero error estimate (A-49) | runs with an exactly zero error estimate | ex1 and ex2 bit-identical |
| `4942bd8` | Picard stop; pVals used (A-47, A-77) | implicit pressure with `maxBeIteration` > 1; `guessSteadyStateICs = 0` | not reported |
| `185a2ac` | IMEX coupling (A-55) | RK32_WBE, RK43_WBE with thermal or hydraulic coupling | ex5 peak fault T 289.120 -> 289.141 K; peak V and end time unchanged to 6 and 10 digits |
| `9745bb2` | Second mechanism; guess units (A-41, A-42) | `wDislCreep2 = yes`; initial viscosity guess | not reported; ex4 at 101 x 121 converges in 44 iterations to the same tolerance |
| `b42c53d` | Inverse viscosity; wet by default (A-43, A-44) | `wDissPrecCreep = yes` | not reported |
| `6b3dad7` | Cohesion, quasi-dynamic (A-37) | nonzero `cohesionVals` | 5 MPa in ex1: initial tau 38.45 MPa, first event at the same time |
| `fbe6ee0` | Step acceptance (A-50) | adaptive RK runs that reach minDeltaT | ex2 peak V 6.351 -> 6.346 m/s, final slip 6e-7 relative; ex1 unchanged |
| `c366a1b` | Adaptive damping (A-45) | steady-state viscosity iteration | converged solution unchanged; runs that stopped unconverged now converge or warn |
| `d823c7b` | Error norms, PID, FEuler (A-52 to A-54) | `max_relative`, `L2_relative` with zero fields, PID in new integrators | ex1 and ex2 bit-identical |
| `28b40f9` | initTime, dynamic-phase 2D stride (A-67, A-68) | power-law qd_fd; qd_fd 2D output | not reported |
| `e77dcda` | Top-boundary coordinates (A-46) | power law with the atan top-boundary loading | not reported |
| `b78b1cb` | Hand-off; limit_qd fallback (A-56, A-75) | quasidynamic_and_dynamic | ex2: first-event peak V 27.11 -> 26.57 m/s |
| `574a54b` | Viscous dissipation, marked PHYSICS CHANGE (A-31) | power law with viscous shear heating or grain-size evolution | not reported; revert this commit alone to reproduce old runs |
| `06bada2` | Geotherm (A-30) | two-point `TVals` | ex2 with TVals = [283 900]: 283 K everywhere -> 283 to 900 K |
| `8954864` | Creeping nodes in dynamic phases (A-38) | dynamic phases with `lockedVals < -0.5` | creeping V 1.67e-9 -> 1e-9 m/s exactly |
| `9866e9b` | disl2 loaded into its own field (A-63) | steady-state restarts with two mechanisms | not reported |
| `04298a8` | Permeability relaxation, error control (A-57) | pore pressure with slip-dependent permeability | kL_p = 1 mm runs now finish except explicit quasi-dynamic (O-01); k within [1e-19, 1e-17] in implicit runs |
| `35ba52f` | Empty timeIntInds left alone (A-57) | explicit pore pressure without `timeIntInds` | restores the behaviour before `04298a8` |
| `aa53679` | Cohesion, fully dynamic (A-37) | quasidynamic_and_dynamic with nonzero `cohesionVals` | 5 MPa: peak V 28.93 -> 25.13 m/s, max slip 20.96 -> 19.72 m |

## 5. Open items and documented limitations

These are not fixed, or are handled only by a startup message or a refusal.

**O-01. Explicit pore pressure in pure quasi-dynamic runs** (`04298a8`). Slip-dependent permeability integrated explicitly is unstable when |V| x dt > about 2.8 x kL_p, and quasi-dynamic steps never go below minDeltaT (one grid wave-transit time, about 14 ms for ex2's grid). A startup note points users to `hydraulicTimeIntType = implicit`, which relaxes permeability exactly; the run with kL_p = 1 mm still diverges in explicit mode (k -> -1e158).

**O-02. `sbpCompatibilityType = compatible` is refused at startup** (`e892279`). Its manufactured-solution errors stay near 0.2 for orders 2 and 4 on constant and variable grids even after the boundary-row typos were fixed (A-35). The root cause, in how the compatible boundary derivative is assembled into the D2/SAT terms, was not found. `fullyCompatible`, the default used by all examples, is unaffected.

**O-03. Top and bottom remote loading stay fixed** (`232bbe2`, `df64c96`). In the quasi-dynamic classes (`source/strikeSlip_linearElastic_qd.cpp`, `source/strikeSlip_powerLaw_qd.cpp`), `momBal_bcT_qd`/`momBal_bcB_qd = remoteLoading` hold that boundary at its initial (steady-state) displacement, because the time-dependent update with vL*t is commented out upstream in `d_dt`. A startup note says so; the physics is left as the authors had it.

**O-04. `timeControlType = PI` is refused** (`232bbe2`). No integrator implements it: `computeStepSize` printed "timeControlType not understood" and asserted at the first adaptive step. P and PID are supported.

**O-05. `systemEvolutionType = steadyStateIts` is impractically slow for ex4.** Every outer iteration runs a full transient simulation to maxTime to find the steady-state fault stress. ex4 at full resolution (301 x 351, maxTime = 1e12 s, time step capped by a Maxwell time of about 300 s) does not finish in practical time on one core: outer iteration 0's viscosity loop converged in 43 iterations, iteration 1's hit its 50-iteration cap, and the transient then took more than 1500 steps of about 266 s while writing 813 MB of 2D output in 28 minutes.

**O-06. The aging law freezes the state variable in two regimes** (`source/fault.cpp`, `agingLaw_psi`). The state rate is set to 0 where b <= 1e-3 and where exp((f0 - psi)/b) overflows, that is where psi is more than about 709 x b below f0. This is a modeling choice; it is documented here and was not changed.

**O-07. With hydraulic coupling, `sNVals` is the total normal stress** (`06bada2`). With hydraulic coupling and gravity, p is the total pore pressure and sNEff = sN - p, so `sNVals` must be the total normal stress. The examples describe `sNVals` as effective stress, which holds only without hydraulic coupling. A startup note says so.

**O-08. `maxStepCount` on restart.** On restart, `maxStepCount` counts the steps taken after the restart, not the total: the mediators add the checkpointed step count to it.

**O-09. MATLAB tools not run** (`f3b5829`). The MATLAB loaders and example scripts (A-93) were rewritten but not run: the MATLAB license on the test machine has expired.

**O-10. `tests/` is not maintained.** It holds old copies of source files that have diverged from `source/` (last changed upstream in March 2019) and cluster scripts such as `tests/memoryLeak/runMazama.sh`. It does not build against the current code. The maintained checks are `tools/regress.sh`, `tools/mms.in` and `tools/checkkeys.py`.

**O-11. Unused and duplicated code.** About half of `source/genFuncs.cpp` is unused (a name search finds 21 of the 43 functions declared in `source/genFuncs.hpp` with no caller outside `genFuncs.cpp`). The five problem classes (`source/strikeSlip_*.cpp`) are largely copy-pasted, so a fix in one usually has to be mirrored in the others; many entries in section 3 touch two or four of them.

**O-12. Per-step Vec allocation** (`8916df7`). `VecDuplicate`/`VecDestroy` pairs in `propagateWaves` and `checkSwitchRegime` run every step. They cost time but do not leak.

**O-13. Constant pore-pressure fields in the time series.** The pore-pressure 1D output writes the constant material fields `beta_p`, `eta_p`, `n_p` and `rho_f` at every output step.

**O-14. `quasidynamic_and_dynamic` with Nz = 1 is refused** (`b78b1cb`). On a 1D grid (a spring slider) the fully dynamic phase runs away (the slip rate grows linearly past 1e4 m/s), in the original code as well; 2D grids work. The run now stops at startup with a message, and `momentumBalanceType = dynamic` with Nz = 1 prints a warning.

**O-15. The off-fault Green's function cannot resume from a crashed run.** `computeGreensFunction_offFault` cannot continue a `G_offFault.h5` that a crashed run left unclosed, because HDF5 cannot open it. Delete the file and rerun.

**O-16. Limit of `tools/checkkeys.py`** (`df64c96`). It cannot detect a key that is read by a different component than the user intended, for example `linSolver`, which only the pore-pressure solver reads.



## 6. Stale statements in scycle-manual.pdf (2020)

The manual (dated May 23, 2020) was re-read for this audit from its text (extracted with pdftotext) and compared with the current code.

- **PETSc version.** The manual says SCycle was tested with PETSc 3.3.2 and 3.12.2. The code needs PETSc 3.14 or newer for the HDF5 timestepping API (`PetscViewerHDF5PushTimestepping`). It was built here with PETSc 3.26.0, real scalars and 32-bit indices, with HDF5 and MUMPS; HYPRE is needed for the AMG and CG solver options.
- **Operator selection.** The manual's key is `sbpType` (mc, mfc, mfc_coordTrans). The code does not read `sbpType`. It reads `sbpCompatibilityType` (`fullyCompatible`, the default; `compatible` is refused, O-02), `gridSpacingType` (`variableGridSpacing`, the default, or `constantGridSpacing`) and `operatorType` (`matrix-based`).
- **Time-step tolerance.** The manual's key is `atol` (default 1e-9). The code reads `timeStepTol` (default 1e-8) and ignores `atol`.
- **`bCoordTrans`.** The manual gives a default of 5 for the coordinate-transformed operators. In `source/domain.cpp` the default is -1, which means no stretching: y = Ly x q, a uniform grid in y, unless a `y` file is given in `inputDir`. A positive value gives y = Ly sinh(b q)/sinh(b) and forces `gridSpacingType = variableGridSpacing`.
- **Steady-state iterations.** The manual lists `steadyStateIts` as a value of `momentumBalanceType`. Domain accepts it, but `main` stops with "unsupported momentumBalanceType" (A-01). Steady-state iterations are selected with `systemEvolutionType = steadyStateIts` and are implemented for power-law runs only.

Other differences found while comparing the manual with the code:

- Momentum-balance solver: the manual's `linSolver` is read only by the pore-pressure solver; the momentum balance reads `linSolverSS` and `linSolverTrans`. The heat-equation key `linSolver_heateq` is as documented.
- Power-law parameters: the manual's `AVals`, `BVals`, `nVals` are now `disl_AVals`, `disl_QRVals`, `disl_nVals` with matching `Depths` keys (`disl2_...` for the second mechanism).
- `rootTol`: the manual calls it a relative tolerance with default 1e-9. It is an absolute tolerance on the stress residual in MPa, default 1e-12.
- `normType`: the code accepts `L2_absolute` (default), `L2_relative` and `max_relative`.
- `timeControlType`: the manual lists PI, which is refused (O-04).
- `minDeltaT`: the manual's default is 1e-3 s. If the key is not set, the code uses one grid wave-transit time, min(dy/cs, dz/cs).
- `lockedVals`: the manual's thresholds are garbled. In the code, values > 0.5 hold V = 0, values < -0.5 creep at vL, and friction decides in between.
- Regime switching in quasidynamic_and_dynamic runs: the manual describes the ratio R = max(eta V / tau_qs). Since `74a132f` the code compares the maximum slip velocity with `trigger_qd2fd`, `trigger_fd2qd`, `limit_qd` and `limit_fd`, which are slip velocities in m/s.
- Fixed-point defaults: the manual gives `fss_T` 0.1, `fss_EffVisc` 0.1 and `gss_t` 1e-8; the code uses 0.15, 0.2 and 1e-10.
- Radiogenic heating: the manual's `he_LVals`/`he_LDepths` do not exist; the code reads a scalar `he_Lrad` (km, default 10 since `0e844ac`).
- Pore pressure: the manual's `maxBeDifference` is `minBeDifference` in the code (default 0.01). `sNVals` is the effective normal stress only without hydraulic coupling (O-07).
- Examples: the manual describes ex3 as an ice-stream simulation; the current `examples/ex3.in` is a quasi-dynamic and fully dynamic earthquake-cycle example. `examples/ex5.in` (flash heating) is not described.

## 7. How to re-verify

Environment used for this audit:

```bash
export PETSC_DIR=$HOME/opt/petsc-3.26.0
export PETSC_ARCH=arch-darwin-c-debug
```

Commands, from the repository root:

```bash
make -C source                                   # builds source/main; add WERROR=1 to fail on warnings
tools/regress.sh compare data/regress-baseline   # ex1 and ex2 against the stored baseline
tools/regress.sh baseline <dir>                  # store a new baseline before a change
tools/regress.sh compare <dir>                   # rerun and compare after it; exit status 1 on any difference
./source/main tools/mms.in                       # MMS convergence on Ny = Nz = 21, 41, 81
tools/checkkeys.py examples/*.in                 # keys no component reads (also BP1.in, tools/mms.in)
```

Notes:

- `tools/regress.sh` needs `h5diff` on the PATH. `REGRESS_NP` sets the number of MPI ranks, `REGRESS_CASES` the examples (default `ex1 ex2`), `REGRESS_DELTA` an absolute tolerance (default: exact), `REGRESS_WORK` the scratch directory and `SCYCLE_BIN` the executable. Each case starts fresh, and checkpoint.h5 is not compared.
- `tools/mms.in` prints the L2 errors of u and sigma_xy and their rates. Expected for order 4: about 3.5 for u and 2.5 for sigma_xy (errors 4.3e-9 and 1.2e-4 at Ny = 81); for order 2: about 2.2 for u. Change `order` and `gridSpacingType` in the file to test the other cases.
- For a leak check on a debug PETSc build, add PETSc's `-objects_dump` option after the input file, as in `8916df7`.
- Configurations that ex1 and ex2 do not cover were checked with the spot runs described in each commit message; repeat those for changes in the matching code.

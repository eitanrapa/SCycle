# Two-fault extension: design note

Status: design agreed, not yet implemented. Stage 0 is done on branch `audit/fixes-2026-10`.
Line numbers refer to upstream commit `74a132f` and will drift; function names are the stable
reference.

## 1. Goal

Determine with forward models whether two neighbouring sub-parallel strike-slip faults that share
the far-field load trade slip, and under what conditions, and produce the matching surface-velocity
predictions (SAF–SJF, CF–AF). Three stages of increasing physics:

1. **Elastic.** Two rate-and-state faults in an elastic half-space. Verify stress transfer and slip
   partitioning, and recover the single-fault limit when one fault is suppressed.
2. **Power-law viscoelastic.** Off-fault power-law creep only. Test whether trading emerges from
   viscoelastic relaxation alone.
3. **Reversible strength change.** Grain-size evolution, fault valving (pore pressure) and shear heating,
   each in turn.

Deliverable: a library of two-fault models spanning the trading regimes, each with a long-term
partitioning history and present-day surface-velocity profiles.

## 2. Geometry: full domain, both faults interior

Today the single fault is the left boundary `y = 0` with the `symmFault` condition `u(0+) = slip/2`,
which is an antisymmetric half-space: the solution satisfies `u(-y) = -u(y)`.

Adding a second fault at `y = +y_f` to that model silently adds its mirror image at `-y_f`, with the
same slip. The model becomes a three-fault system whose long-term rates satisfy `V1 + 2 V2 = vL`,
not `V1 + V2 = vL`. The `rigidFault` option avoids the image but doubles the stiffness of fault 1 and
flattens the surface-velocity profile on one side. Neither represents SAF–SJF.

**Decision:** model the full domain `y ∈ [-L1, L2]` with **both** faults interior and treated by the
same code. Load both far boundaries with remote displacement:

    bcL = -vL t / 2 + bcLShift        bcR = +vL t / 2 + bcRShift

Consequences:
- Unknowns double relative to the half-space. The momentum-balance matrix is still factored once
  (MUMPS Cholesky with a reused preconditioner), so cost per step grows roughly linearly.
- The built-in `sinh` grid stretching (`Domain::setFields`) refines only near `y = 0`. Supply the
  y-coordinates from a file instead (the existing `inputDir + "y"` PETSc-binary hook in
  `Domain::setFields`), refined near both faults. Metrics are computed numerically from `y(q)`
  (`SbpOps_m_varGrid::constructJacobian`), so any smooth monotone `y(q)` works.
- Each fault must sit at the midpoint between two grid lines `i_f` and `i_f + 1`, at least 6 nodes
  from the left and right boundaries (the Neumann closure of `D1` spans 6 columns).
- The left-boundary assert in `StrikeSlip_LinearElastic_qd::checkInput` currently rejects
  `momBal_bcL_qd = remoteLoading`; `d_dt` must update `bcL` like `bcR`.

## 3. Mechanics of an interior fault: jump-corrected right-hand side ("lift")

There is no interface, multi-block or jump support anywhere in the code today. Two routes exist.

### Route B (chosen for stages 1–3): single grid, jump-corrected right-hand side

The momentum-balance operator is

    A = A_y + A_z + A_L + A_R + A_T + A_B
    A_y = H z_r (D_q μ q_y D_q − H_y⁻¹ R_y)       (y part, built in constructDyymu)
    A_z = H y_q (D_r μ r_z D_r − H_z⁻¹ R_z)       (z part, built in constructDzzmu)

plus SAT boundary terms. For fault k with slip `δ_k(z)` define the step field

    U_k = δ_k(z) · step_k(y),   step_k = 1 for grid rows iy > i_f(k), 0 otherwise
    Ũ_k = δ_k(z) for every row (its smooth extension)

and solve

    A u = rhs_BC + Σ_k J_k,       J_k = A U_k − diag(step_k) A Ũ_k  (= A_y U_k)

The identity `J_k = A_y U_k` is exact because `D_q` and the undivided `R_y` stencils annihilate
constants (so `A_y Ũ_k = 0`), the z-operators and the top/bottom SATs act within one y-row (so they
commute with `diag(step_k)`), and the left/right SATs see a constant field when the fault is at least
6 nodes away. `J_k` costs two `MatMult`s with the existing `A`; **`A` and its factorization are
unchanged**, and with `δ_k ≡ 0` the code path is bit-identical to today's.

Then:
- The physical displacement is `u`. The continuous part is `w = u − Σ U_k`.
- Every **y-strain** must be computed from `w`, not `u`: `σ_xy = μ D_y w`
  (`LinearElastic::computeStresses`, `PowerLaw::computeTotalStrains`). The z-strain keeps `u`
  (the jump `∂z δ_k` is physical).
- Fault traction: average `μ D_y w` on rows `i_f` and `i_f + 1` (two row scatters).
- Accuracy: `w` is C¹ but not C² across the fault when slip varies with depth
  (`[w_yy] = −(μ δ')'/μ`). Displacement converges at O(h²); traction at O(h), with error
  ≈ 0.24·h·|[w_yy]|. With `δ` linear in depth the traction is O(h²) (a useful test).
- If O(h) traction is not good enough, lift the kink too ("B+"): add `½ c(z) (y − y_f)²` with
  `c = −(μ δ')'/μ` to `U_k` for rows past the fault; traction becomes O(h²). Implement only if
  the convergence test in stage 2 requires it.
- Rejected alternative: representing fault slip as an eigenstrain in one column of cells through the
  power-law viscous source term. It omits the `R_y` contribution and leaves an O(δ/h) error in the
  fault traction.

Limits of Route B:
- **μ must be continuous across each fault** (no bimaterial contrast).
- **Not for fully dynamic phases** (`_fd` and `_qd_fd` mediators); those need characteristic
  interface conditions.

### Route A (deferred): split-node multi-block SBP-SAT

Blocks concatenated in y keep the global index layout `Ii = iy·Nz + iz`. Needs per-block 1D
operators, interface selectors, symmetric energy-stable interface SATs for `[u] = δ` and
`[μ u_y] = 0`, a fifth data vector in `setRhs`, and a per-block `R` term. References:
Kozdon, Dunham & Nordström (2012); Erickson & Day (2016); Kozdon, Erickson & Wilcox (2021,
hybridized SBP-SAT for multiple faults). Needed only for a material contrast across a fault, fully
dynamic ruptures, or if Route B traction accuracy proves insufficient. Encapsulate Route B behind a
small interface (`InteriorFaultLift`: `setSlip`, `addToRhs`, `strainField`, `traction`) so Route A
can replace it without touching the mediators.

## 4. Code changes by stage

### Stage 0: regression baseline (done in branch `audit/fixes-2026-10`)
- `tools/regress.sh baseline|compare <dir>` runs `examples/ex1.in` and `examples/ex2.in` and
  compares every HDF5 dataset with `h5diff`. Every later stage must keep these bit-identical in
  single-fault mode. A baseline from the end of the audit is in `data/regress-baseline/`.
- `tools/mms.in`: manufactured-solution convergence of the elastic quasi-dynamic solver on
  Ny = Nz = 21, 41, 81 (order 4: u converges at about 3.5, sxy at about 2.5; order 2: about 2).
- `SEAS_benchmarks/BP1`: runnable (`createICs.py` writes the grid and initial conditions; 25 m
  spacing gives 401 x 2401 nodes). This is the reference for gate 2d.
- The partition-dependent coefficient of the 4th-order `R` operator is fixed (`b65b5e0`): 1 and
  2 ranks agree to about 1e-10, so convergence tests may run in parallel.
- `tools/checkkeys.py` lists input keys that no component reads; run it on every new input.
- 94 findings were fixed first (`docs/AUDIT.md`), including ones on the paths stages 2 to 4
  rely on: the IMEX coupling of temperature and pressure, permeability stability through events,
  creeping nodes in fully dynamic phases, viscous dissipation in the heat and grain-size equations,
  and atomic checkpoints for long runs.

### Stage 1: generalize the fault class, still one fault, bit-identical
- `Fault` / `Fault_qd` already take the body-to-fault scatter as a constructor argument.
  Add a name (default `"fault"`) that sets:
  - integrand keys: `slip`, `psi` for the default fault, `<name>_slip`, `<name>_psi` otherwise
    (hard-coded today in `Fault_qd::initiateIntegrand`, `updateFields`, `d_dt`),
  - HDF5 groups `/<name>` and `/<name>_qd` (hard-coded `/fault`, `/fault_qd`),
  - per-fault radiation-damping scale (replaces the global `_faultTypeScale` in the `Fault_qd`
    constructor; `η = √(μρ)/2` for an interior fault between identical materials),
  - per-fault creep rate for `lockedVals < −0.5` (today `_D->_vL`, in both `ComputeVel_qd` and
    `Fault_fd::d_dt`).
- Mediator holds `std::vector<Fault_qd*>` and one scatter per fault; `HeatEquation::be` and
  `PressureEq::dk_dt` take the slip-rate key instead of `"slip"`.
- `Domain`: add `makeRowScatter(iy)` next to the `body2L` construction in `setScatters`.

### Stage 2: one interior fault
- New `InteriorFaultLift` (Route B). `LinearElastic::setRHS` adds `J`; `computeStresses`
  differentiates `w`.
- Full-domain loading (`bcL = −vL t/2`), grid from file.
- `solveSS` / `solveSSb` (steady-state initial guess) impose Neumann `bcL = τ_ss` and only work
  for a boundary fault. Disable them for interior faults; initialise ψ at steady state for an
  assumed slip rate instead.

### Stage 3: two interior faults, elastic
- Two lifts, two `Fault_qd` objects, two scatters; qd step-size control already picks up any new
  integrand key when `timeIntInds` is empty (`OdeSolver` uses all explicit keys). Inputs that list
  `timeIntInds` must add `<name>_slip`, `<name>_psi` themselves (`PressureEq::addErrorControl`
  shows the pattern for automatic additions to a non-empty list).

### Stage 4: viscoelastic and reversible-strength physics
- `PowerLaw::computeTotalStrains` uses `w` for `gTxy`; the viscous source term `B·γV` is unchanged.
- Heat: frictional source becomes a sum of per-fault Gaussians centred at `y_f(k)`
  (`HeatEquation` builds one centred at `y = 0`, normalised for a half-space); reject the
  `w = 0` boundary-flux option for interior faults.
- Grain size: no fault-specific code.
- Pore pressure / valving: one `PressureEq` per fault. Use `hydraulicTimeIntType = implicit`:
  permeability is then relaxed exactly over each step (`PressureEq::relaxPermeability`); the
  explicit form is unstable during events when `kL_p` is small.

## 5. Verification gates

| Stage | Test | Pass criterion |
|---|---|---|
| 1 | ex1, ex2 via `tools/regress.sh` | bit-identical |
| 2a | Prescribed smooth `δ(z)`, μ uniform, fault at `y = 0` of a symmetric full domain, vs the existing `symmFault` half-space run with `bcL = δ/2` | `u` O(h²), traction O(h) over three grids |
| 2b | Same with `δ` linear in depth | traction O(h²) |
| 2c | Analytic antiplane screw dislocation with free-surface image | `u` O(h²) |
| 2d | Rate-and-state on the interior fault, other boundary far | recurrence interval and peak V within 1% of the half-space BP1-style run |
| 3a | Fault 2 locked (`lockedVals`) | reproduces 2d |
| 3b | Fault 2 velocity-strengthening | long-term `V1 + V2 = vL` from cumulative slip |
| 3c | Coseismic stress change on fault 2 from a fault-1 event | matches the 2D screw-dislocation formula at distance `y_f` |
| 4 | Power-law single-fault limit | matches the existing `StrikeSlip_PowerLaw_qd` run; heat budget `∫Q = Σ τ_k V_k` |

Keep μ uniform across both faults (Route B does not handle a modulus jump at a fault). The
partition-dependent `R` coefficient that once made results depend on the rank count is fixed.

## 6. Output for the research goal

- Per fault, per 1D stride: slip, V, τ, τ_QS, ψ (already written by `Fault::writeStep`), plus scalar
  series for max V, moment rate, cumulative slip at a reference depth, and the partition fraction
  `Δslip_k / Σ Δslip` over a trailing window.
- Surface velocity: write `surfVel` directly. Because the right-hand side is linear in
  `(bcL, bcR, δ1, δ2)`, one extra back-substitution with the factored `A`,
  `A v = rhsL(−vL/2) + rhsR(+vL/2) + J(V1) + J(V2)`, gives the instantaneous velocity field; scatter
  its top row. This avoids differencing displacement across adaptive time steps.
- Context file: fault positions `y_f`, grid rows `i_f`, local grid spacing, so interseismic profiles
  can be selected in post-processing by a max-V threshold.

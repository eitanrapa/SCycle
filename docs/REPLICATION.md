# Replicating Allison & Dunham (2018, 2021)

The two papers this code was written for, run again with the current code (the fork's `master`):

- **AD18**: Allison, K. L., & Dunham, E. M. (2018). Earthquake cycle simulations with rate-and-state
  friction and power-law viscoelasticity. *Tectonophysics* 733, 232-256.
- **AD21**: Allison, K. L., & Dunham, E. M. (2021). Influence of shear heating and thermomechanical
  coupling on earthquake sequences and the brittle-ductile transition. *JGR Solid Earth* 126,
  e2020JB021394.

The runs are in `data/replication/runs` (a `tools/batch.sh` batch; `data/` is ignored by git) and the
inputs come from `examples/replication/make_ad2021.py` and `make_ad2018.py`. Started 2026-10-08.

## 1. Sources

**AD21.** The paper's inputs, grids, steady-state initial conditions, results and a spreadsheet of
summary numbers are public at osf.io/nzmw2 (downloaded to `data/replication/ad2021/osf`: inputFiles
1.05 GB, simulationData 3.8 GB). The code it used is very likely upstream `e9b1aa9` (2020-11-06,
"updating files for heat paper submission", 13 days before submission); the results' metadata
(`/med`, `/he`) record the settings actually used, which differ from the posted inputs in the error
norm, the tolerance and the step bounds.

**AD18.** No code, inputs or data were published. Upstream's history has inputs of the time:
`max_cremeBrulee_cycle.in` of June 2017 (`dd0c35c`) matches the paper in domain (120 x 40 km), dc
(16 mm), friction (Blanpied et al. 1991 and 1995, set with 30 K/km) and effective normal stress
(22.3 MPa/km) but has the layered rheology the paper rejected; `in/test.in` of January 2018
(`28e9f45`, "compare fixed point iteration method with spin up method used in Allison and Dunham
Tectonophysics") has the paper's quartz-diorite flow law. The rest comes from the paper: the mixed
quartz-diorite and wet-olivine flow law (Ji et al. 2003, volume fractions of its Fig. 5), the
geotherms, the y mapping and the depth spacing of its Fig. 21.

## 2. Inputs for the current code

AD21's inputs need these changes (`make_ad2021.py`):

| Paper's input | Now | Why |
|---|---|---|
| `normType = max_relative`, `timeStepTol = 1e-8` | `L2_absolute`, `1e-7` | the values the results record; the implicit-explicit integrator takes only L2 norms |
| transient heating implied by `heatEquationType` | `evolveTemperature = 1` | required since 2022 |
| `guessSteadyStateICs = 1` | and `computeSSMomBal = 1` | what it did in 2020 (A-88) |
| ssT: coupled heat with the explicit RK43, so T stayed at its initial value | `thermalCoupling = no` | the power law keeps the T file (ambient plus steady-state anomaly); the current code refuses the old pairing |
| `momentumBalanceType = steadyStateIts` (steady-state runs) | `systemEvolutionType = steadyStateIts` | A-01 |
| `maxDeltaT = 2.5e5` s (elastic run) | `1e11` s | the cap only adds output between earthquakes and would make an elastic cycle take 18 hours here |

AD18 is reconstructed (`make_ad2018.py`): grid 201 x 492 (y = 120 km sinh(10q)/sinh(10); depth
spacing of Fig. 21, 30 m at 12-18 km), mu 30 GPa, rho 3 g/cm^3, dc 16 mm, f0 0.6, V0 1e-6 m/s, a and b
and the normal stress of the June 2017 input, quartz-diorite (A 1.3e-3 MPa^-n/s, n 2.4, Q 219 kJ/mol)
over wet olivine (A 3.6e5, n 3.5, Q 480) mixed per node, T = 293 K + G z with G = 20, 25, 30 K/km held
fixed, RK32 with P control at 1e-8 on slip and state, minDeltaT 1e-3 s, maxDeltaT 2e6 s. Unknowns
inferred: the surface temperature (293 K), Ny (from the June 2017 input), the volume fraction (a tanh
fit to Fig. 5).

**Spin-up.** AD21's published time series all start after long spin-ups: about 950 m (no shear
heating), 1240 m (ssT) and 2250 m (elastic) of slip had accumulated at 5 km, some 110 to 270 cycles,
with the clock reset on the way. From the code's steady state (`guessSteadyStateICs`,
`computeSSMomBal`) our runs start with much shorter cycles (ssT: 234, 243, 252 yr; no heating: 77, 95
yr) and would need weeks per case here. So AD21's checks start from the paper's own states instead
(`examples/replication/seed_ad2021.py`): the heating cases from the first 2D snapshot of its w = 0.3 m
run (158 s after its start, at the end of an earthquake; its w = 0.3 and 1 m runs start from the same
state, slip within 3 mm) and the elastic case from the first output of its elastic run. Started from
these states, the current code reproduces the paper's fault state: shear stress within 0.15 MPa
(heating) and 0.02 MPa (elastic) at every depth, slip rates within 1-3% and 0.3% (our first output at
1 ms, theirs at 147 and 921 s). Each check then runs one or two cycles. The case without shear
heating has no published 2D state and is spun up from the steady state (to 30000 yr). AD18 starts
from the steady state.

## 3. Runs

Two batches, one job each (the stage 5 matrix has the other two slots):

| Batch | Run | Case | Length | Compare with |
|---|---|---|---|---|
| `data/replication/seeded` | `q1_ad21_seed_ssT` | AD21 reference (LAB 50 km, lambda 0.37), anomaly at steady state | 800 yr | the paper's ssT cycle (354.2 yr) |
| | `q2_ad21_seed_le` | AD21 elastic | 600 yr | the paper's next earthquake at 266.81 yr |
| | `q3_ad21_seed_w1` | AD21 transient shear heating, w = 1 m | 370 yr | the paper's w = 1 m cycle (earthquake at 352.84 yr, rise 321 K) |
| | `q4a_ad18_le`, then `q5a_ad18_g30`, `q5b_ad18_g25`, `q5c_ad18_g20` | AD18 elastic, then 30, 25, 20 K/km | until the recurrence converges (at most 15000 yr) | the paper's text and figures |
| `data/replication/runs` | `r3_ad21_nsh_zlab50_l0.37` | AD21 reference without shear heating, from the steady state | 30000 yr | the paper's no-heating cycle (271.85 yr) |

**Stopping rule for AD18** (agreed with the user 2026-10-17): each run stops once its recurrence
interval changed by less than 0.5 yr over two consecutive cycles (`tools/converge_stop.py`, running
detached, log `data/replication/seeded/converge_stop.log`), which marks it finished so the queue moves
on; when the stage 5 matrix has finished it restarts the seeded batch with three jobs, so the
viscoelastic cases run in parallel. Its earthquakes run at the 1 ms step floor of the 2017 input, about
60,000 steps (7 hours) per cycle. The elastic run's recurrence by 5475 yr: 260.5, 293.7, 304.1, ...,
333.3, 334.6, 335.8 yr, changing 1.2 yr per cycle (the paper: 366 yr).

Set aside in `data/replication/runs` (inputs renamed, checkpoints kept): the ssT spin-up from the
steady state (`r1`, stopped at 821 yr), the elastic spin-up with the paper's 2.5e5 s step cap (`r2`) and
without it (`r5`).

## Results so far (2026-10-11)

Measured with `examples/replication/events_ad2021.py`, which applies the paper's definitions to our
runs and to its own time series alike (earthquake: maximum slip rate above 1 mm/s; down-dip limit: 98%
of the coseismic potency):

| | Recurrence (yr) | Nucleation (km) | Max slip (m) | Down-dip limit (km) | log10 M |
|---|---|---|---|---|---|
| Elastic, seeded (`q2`), 2nd event | 269.08 | 14.54 | 8.418 | 17.32 | 15.6554 |
| Elastic, paper (data / spreadsheet) | 266.21 | 14.52 / 14.30 | 8.420 / 8.42 | 17.33 / 17.31 | 15.6548 / 15.654 |
| ssT, seeded (`q1`), events 1 and 2 | 361.86 | 12.46, 12.49 | 11.288, 11.346 | 14.40, 14.39 | 15.7034, 15.7053 |
| ssT, paper (data / spreadsheet) | 355.47 | 12.46 / 12.66 | 11.276 / 11.28 | 14.39 / 14.40 | 15.7031 / 15.703 |

From the paper's own states, the current code gives the same ruptures: maximum slip, down-dip limit and
moment within 0.1-0.6%, nucleation at the same depth. The recurrence intervals come out 1.1% (elastic)
and 1.8% (ssT) longer; the first elastic earthquake after the seed comes 0.87 yr after the paper's
(267.68 against 266.81 yr; its first event is also 1.3% smaller, the second matches). The elastic run
did not cap its steps at 2.5e5 s as the paper's did, but the ssT run used the paper's step bounds and
differs as much, so the cap is not the cause; candidates are the code changes since 2020 (section 4).

Spin-up without shear heating (`r3`), from the steady state: recurrence 94.6, 136.8, 152.3, 161.2,
167.2, 172.1, 176.3 yr over its first eight cycles (the paper's limit cycle: 271.85 yr), about 12,000
steps per cycle.

## 4. What should differ, and why

Fixes since the papers that touch these configurations (`docs/AUDIT.md`):

- **A-31, viscous dissipation: no difference.** Upstream's code had 0.5 x the dissipated work from July 2019
  and sqrt(2) x from the submission commit `e9b1aa9` on, but the paper's runs used neither: its published
  results balance a source of exactly the dissipation (factor 0.998 in the viscous-only steady state, from
  -k lap(dT) against the stored source; 0.999 in all 56 interseismic snapshot pairs of the transient run
  with w = 0.3 m, from rho c dT/dt - k lap(dT) against sdev^2/effVisc). The fixed code does the same, so
  no run with the old factor is needed (a build with sqrt(2) restored is in data/replication/build_sqrt2).
- **A-50 and the step floor.** The paper's earthquakes ran at its 1 ms floor (the results' smallest
  steps), where steps are accepted whatever their error; A-50 made such steps consistent (time and
  state advance with the same step). Same floor here.
- **A-55, implicit-explicit coupling.** It lagged one step in the transient-heating runs.
- **A-33, A-34**, the 4th-order R coefficient (4e-8 relative in ex2); **A-23, A-24**, the root finder
  (converged velocities unchanged); **A-45**, the damping of the steady-state viscosity iteration
  (initial conditions only).
- Not applicable: A-27 (variable grids here), A-28 (the heat solver was CG), A-29 (w > 0), A-30
  (temperatures from files), A-95 (minDeltaT given), A-97 (w resolved by 6 cm cells at the fault).

Upstream's other changes between the papers' code and `74a132f` were not audited for their effect on
results; the comparison will show them if they matter.

## 5. Targets

**AD21**, reference case, from the paper's `summaryData.xlsx` (the heating values hold for every w and
for ssT):

| | Recurrence (yr) | Nucleation (km) | Down-dip limit (km) | BDT (km) | Max slip (m) | log10 M |
|---|---|---|---|---|---|---|
| Elastic | 266.21 | 14.30 | 17.31 | | 8.42 | 15.654 |
| No heating | 271.85 | 14.26 | 17.34 | 19.49-21.00 | 8.76 | 15.672 |
| Shear heating | 355.47 | 12.66 | 14.40 | 13.26-15.09 | 11.28 | 15.703 |

Maximum transient temperature rise (LAB 50, lambda 0.37): 657, 321, 226 and 193 K for w = 0.3, 1, 3
and 10 m. The results' time series (`data_eqs/*.mat`) allow comparing slip, slip rate, stress and
state at every depth.

**AD18** (text and figures): recurrence 366 yr elastic, 366, 375 and 365 yr at 20, 25, 30 K/km; about
10 m of slip; nucleation at 10 km and coseismic slip to 14 km in every run; brittle-ductile
transition zone 14-17 km at 30 K/km and 16-20 km at 25 K/km; fault creep to the bottom at 20 K/km; the
surface displacement 20 km from the fault 10 yr after an earthquake 0.188, 0.148 and 0.138 m at 30, 25
and 20 K/km (Fig. 18).

## 6. Cost on the Mac mini

AD21's grid (301 x 1295, 390,000 nodes) costs 4.6 s per step with transient heating and 1.9 s per step
elastic, with the stage 5 matrix running beside it; its steady-state viscosity iteration takes about
an hour. The transient-heating cycle (steps capped at 2.5e5 s, as in the paper) is about 47,000 steps,
2.5 days. AD18's grid (201 x 492) is a quarter the size.

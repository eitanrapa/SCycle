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

Both start from the code's steady state (`guessSteadyStateICs`, `computeSSMomBal`) and run until the
earthquake cycle repeats (AD21 needed 50-100 cycles from its steady state; AD18 used its own
four-step spin-up).

## 3. Runs

| Run | Case | Length | Notes |
|---|---|---|---|
| `r1_ad21_ssT_zlab50_l0.37` | AD21 reference (LAB 50 km, lambda 0.37), anomaly held at steady state | 7500 yr | spin-up for the transient-heating cycle; the paper's ssT data are at 4843-5197 yr |
| `r3_ad21_nsh_zlab50_l0.37` | AD21 reference without shear heating | 7500 yr | the paper's data: 270-540 yr |
| `r4a_ad18_le` ... `r4d_ad18_g30` | AD18 elastic, 20, 25, 30 K/km | 15000 yr | |
| `r5_ad21_le_l0.37` | AD21 elastic | 2000 yr | the paper's data: the first cycle, event at 266.8 yr |
| (to come) `ad21_w1` | AD21 transient heating, w = 1 m | one cycle | from the end of r1 just after an earthquake (`tools/ckpt_to_ics.py`) |

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

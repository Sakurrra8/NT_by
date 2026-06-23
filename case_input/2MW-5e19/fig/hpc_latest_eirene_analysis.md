# Latest HPC NT versus SOLPS/EIRENE

The comparison uses HPC jobs 10087195-10087203 and the triangle-grid
SOLPS/EIRENE references in each `case_input/2MW-*e19` directory. NT and
reference panels in the PDFs use the same colorbar.

## Density scan

| Case | D relative L1 | D total bias | D2 relative L1 | D2 total bias |
|---|---:|---:|---:|---:|
| 3.2e19 | 0.520 | -24.2% | 0.415 | +12.1% |
| 3.6e19 | 0.334 | -18.8% | 0.657 | +50.6% |
| 4.0e19 | 0.291 | -15.8% | 0.800 | +68.9% |
| 4.6e19 | 0.230 | -15.7% | 1.017 | +97.0% |
| 5.0e19 | 0.218 | -18.5% | 1.031 | +100.3% |
| 5.4e19 | 0.214 | -19.6% | 1.013 | +100.4% |
| 5.5e19 | 0.221 | -20.4% | 0.996 | +98.9% |
| 5.6e19 | 0.216 | -19.9% | 0.976 | +97.2% |
| 5.8e19 | 0.219 | -20.7% | 0.949 | +94.6% |

The atomic D distribution approaches the EIRENE distribution as density
increases, but its integrated inventory remains about 16-24% lower. The
molecular D2 inventory changes in the opposite direction and reaches about
twice the EIRENE value from 4.6e19 upward. This anti-correlated D/D2 shift
points first to the target/wall recycling partition, molecular survival, and
reflection launch model rather than a common normalization error.

## Detailed 5.0e19 comparison

| Quantity | Relative L1 | Total bias | Median NT/EIRENE |
|---|---:|---:|---:|
| D density | 0.218 | -18.5% | 0.883 |
| D2 density | 1.031 | +100.3% | 2.118 |
| D2+ density | 1.324 | +61.3% | 1.600 |
| Electron-impact D+ source | 0.496 | +0.25% | 0.787 |
| Raw `Diss1+Diss2` versus `mol_ds_2D.data` | 188.644 | +18857% | 189.918 |

The electron-impact D+ production totals agree to 0.25%, while the relative
L1 error is about 0.50. Its main discrepancy is therefore spatial
redistribution, not global source normalization.

The raw D2 dissociation comparison is not presently a valid like-for-like
benchmark. NT gives:

- sampled `Diss1`: 3.9051e22 s^-1;
- sampled `Diss2`: 5.8168e20 s^-1;
- total: 3.9632e22 s^-1.

An independent reconstruction using the NT D2 density, local `ne`, `Te`,
triangle volume, and the same AMJUEL H.4 reaction 2.2.5g gives
3.8018e22 s^-1. The sampled `Diss1` total is only 2.7% above this
reconstruction. In contrast, `mol_ds_2D.data` sums to 2.0906e20 s^-1.
Even reconstruction with the ground-state H.4 reaction 2.2.5 and the EIRENE
D2 density gives 1.3242e22 s^-1. Thus `mol_ds_2D.data` likely contains a
reaction/stratum subset or a different normalization and should not be used
as the total molecular dissociation source until its EIRENE tally definition
is identified.

## Interpretation of the latest switches

- The HPC run used `K_DWTrimReflection=1`, specular reflection, and a total of
  50000 target-recycling histories. These directly affect the D/D2 launch
  partition and wall-adjacent spatial distribution.
- `K_D2Flight=0`. The reported D2+ density is the local-balance result; the
  D2+ flight diagnostics contain zero tracked steps. This run therefore does
  not test the new H.3 drift-corrected D2+ orbit velocity model.
- Variance reduction was disabled (`K_Roulette=0`, `K_Splitting=0`), so the
  observed scan trend is not caused by regional splitting.
- Disabling plasma collisions outside the plasma region is physically
  consistent with the intended vacuum treatment. It can change residence
  near the mesh boundary relative to older NT runs, but it does not explain
  the internally consistent D2 dissociation rate.

The next decisive comparison is to export the exact EIRENE reaction-resolved
H.4 2.2.5/2.2.5g and 2.2.10 tallies, summed over the same strata, and to run
one matched case with the legacy and database wall-reflection switches while
keeping the target source weights identical.

## EIRENE surface-model matching test

The EIRENE manual and `input.dat` identify `EWALL=-0.1 eV` as a
Maxwellian-flux thermal emission model with mean emitted energy `2*Tw=0.2 eV`.
They also specify `ERMIN=1 eV`, TRIM D-on-W reflection, and total recycling
coefficient one. NT was changed to reproduce these settings, including thermal
reemission of incident D2 and the local target incidence angle.

HPC ablation jobs were:

- 10143866: thermal emission plus ERMIN, fixed 60 degree target angle;
- 10143868: thermal emission plus ERMIN and local target angle;
- 10143124: the above plus target incident-energy profiles (`K_Ei=1`);
- 10143976: thermal emission without ERMIN.

The final merged-code validation is job 10144428 (`K_Ei=2`, local target
angle). Relative to the previous baseline, atomic D improves from relative L1
0.218 to 0.191 and total bias from -18.5% to -14.9%. The D+ electron-impact
source remains close in total (bias +1.8%) and has relative L1 0.522.
Molecular D2, however, worsens from total bias +100.3% to +148.8%.

The source audit shows that the manual-consistent surface treatment shifts
recycling from fast D to thermal D2 and thermalizes repeatedly incident D2.
This increases D2 residence time and source weight. Therefore the remaining
D2 discrepancy is upstream of the flight estimator: the meaning and
normalization of `recycled_neutral_flux_D_l/r.data`, and its conversion into
separate D and D2 source channels, must be matched to the EIRENE surface
tallies before further tuning the collision model.

# NT_by

NT_by is a neutral and molecular test-particle transport code for SOLPS/EIRENE
edge-plasma cases.  The current working focus is EAST triangular-mesh
deuterium transport, including neutral D/D2 recycling and optional D2+ test-ion
tracking.

The code is still primarily a Linux/HPC codebase.  Local development works best
through WSL or another Linux environment.

## Repository Layout

```text
src/                      C++ implementation
include/                  shared declarations and global state
Inputfile/settingfile/    run settings
Inputfile/database/       ADAS/EIRENE/AMJUEL-style rate data
case_input/               SOLPS/EIRENE case data and triangular mesh files
triangle/EAST/            comparison and plotting utilities
Outputfile/               generated run output
Makefile                  HPC-oriented build
Makefile.local            local/WSL build
crslum.slurm              Slurm run script
```

`Outputfile/` is generated data.  Keep source, settings, and case input under
version control; avoid committing new run products unless they are deliberately
curated reference results.

## Build And Run

Local WSL/Linux debug build:

```bash
make -f Makefile.local BUILD=debug -j1
```

Default local run:

```bash
./bin/edit Inputfile/settingfile/setting_Trimesh_D_5.log
```

HPC run:

```bash
sbatch crslum.slurm /data/leuven/379/vsc37950/nt/Inputfile/settingfile/setting_Trimesh_D_5.log
```

For an EAST density scan, use the prepared files:

```text
Inputfile/settingfile/setting_Trimesh_D_1.log
...
Inputfile/settingfile/setting_Trimesh_D_9.log
```

These are no-variance-reduction triangular-mesh settings.  `Tri_sbatch.sh`
submits the sweep.

## Current EAST Triangular Inputs

The EAST triangular settings are intentionally kept simple for manual testing:

- `MeshMode=3`
- `K_D2Flight=0`
- `D2ElasticModel=2`
- `K_NNCs=1`
- `K_database_Pra=2`
- `K_Maxwell=1`
- `K_Puff=0`
- `Num_D2_pump=0`
- `numPar_flight_D2=0`
- `K_Roulette=0`
- `K_Splitting=0`
- `Databasepath Inputfile/database/`

This matches the reference EIRENE species card (`NFOL=-1`). D2+ is created by
D2 electron-impact ionization or D2-D+ charge exchange, destroyed immediately
through the local DS branching ratios, and reported with the local `P/L`
density estimator. It is neither injected nor transported as a finite-orbit
test ion in the validation inputs.

## Input Setting Reference

This reference follows
`Inputfile/settingfile/setting_Trimesh_D_5.log`, which is the template for the
current EAST density scan. The settings through `coeff_puff` are read
**positionally**: keep those lines in the same order and do not remove a line
when disabling a feature. In particular, the historical spellings
`K_Rflect`, `K_RflectDirection`, and `coeff_ercyc_wall` are intentional input
labels; their values are stored in the code variables `K_Reflect`,
`K_ReflectDirection`, and the wall recycling coefficient. The later optional
settings are searched by name.

Paths are interpreted from the run working directory unless stated otherwise.
`DBoundaryFluxFile` is relative to `Casepath`, while `DWTrimDatabase` names a
directory below `Inputpath/database/`.

### Paths And Mesh

| Setting | Current value | Effect |
| --- | ---: | --- |
| `Inputpath` | `Inputfile/` | Root for code input tables and geometry support files. |
| `Casepath` | `case_input/2MW-5e19/` | SOLPS/EIRENE case root containing `2D_data/`, `profiles_data/`, and mesh data. This is the only path changed between the density-scan setting files. |
| `Outputpath` | `Outputfile/test/data/` | Destination for generated densities, rates, audits, and optional HDF5 files. |
| `Databasepath` | `Inputfile/database/` | Root of the EIRENE/AMJUEL/HYDHEL/AMMONX reaction database. It must be valid when `K_database=2`; otherwise the particles can be transported without the intended collisions. |
| `N_poloidal`, `N_radial` | `98`, `38` | B2 grid dimensions, including guard cells. They must match the exported SOLPS case. |
| `MeshMode` | `3` | `1`: structured orthogonal grid; `2`: legacy full/external grid; `3`: triangular neutral mesh coupled to B2 cells. The EAST validation uses mode 3. |

### Species And Sources

| Setting | Current value | Effect |
| --- | ---: | --- |
| `K_H`, `K_D`, `K_T` | `0`, `1`, `0` | Enable H, D, and T test-particle systems. The current case transports only deuterium neutrals. |
| `K_back` | `2` | Select the hydrogenic plasma background isotope: `1=H`, `2=D`, `3=T`. |
| `K_DT` | `0` | Enable mixed D-T simulation branches. Keep off for the pure-D validation case. |
| `Ratio_T` | `0` | Tritium fraction of the hydrogenic background used by mixed-isotope branches. |
| `K_CX_DT` | `0` | Enable cross-isotope D-T/T-D charge exchange. It is inactive when `K_DT=0`. |
| `K_Recyc` | `1` | Launch target-recycling D/D2 histories from the prepared target source. |
| `K_Rec` | `1` | Launch neutral histories produced by plasma recombination. This is a volume source, separate from target recycling. |
| `K_Puff[0..4]` | all `0` | Five hard-coded gas-puff/pumping-location toggles. All zero means that the validation case has no puff source. |
| `K_Methane`, `K_C`, `K_Ar` | all `0` | Enable CD4, carbon, or argon source/transport branches. They are excluded from the D/D2 comparison. |
| `numPar_flight` | `50000` | Monte Carlo history budget used for volume/grid sources such as recombination. |
| `numPar_flight_Target` | `50000` | Total history budget shared by target-recycled atomic and molecular sources. |
| `numPar_flight_D2`, `Num_D2_pump` | `0`, `0` | D2 pump history budget and represented source rate. Both must be positive, together with `K_D=1`, to create an artificial D2 pump source. They do not control D2 made by recycling or reactions. |
| `numPar_flight_T2`, `Num_T2_pump` | `0`, `0` | Equivalent T2 pump controls. |
| `numPar_flight_CD4`, `Num_CD4_pump` | `0`, `0` | Equivalent methane-source controls; `K_Methane` must also be enabled. |

### Collision And Transport

| Setting | Current value | Effect |
| --- | ---: | --- |
| `K_database` | `2` | Main hydrogenic collision database: `1=ADAS`, `2=EIRENE/AMJUEL/HYDHEL`. Use 2 for direct EIRENE validation. |
| `K_database_Pra` | `2` | Product-energy/momentum treatment: `1=legacy ADAS-style`, `2=EIRENE H.10/H.8 energy-weighted rates`. |
| `K_Maxwell` | `1` | Background/product velocity sampler: `1=full Maxwellian`, `2=legacy fixed RMS-speed shortcut`. |
| `K_H2_elastic` | `1` | Master switch for D2-D+ molecular elastic collisions. Setting it to zero forces `D2ElasticModel=0`. |
| `D2ElasticModel` | `2` | `0`: off; `1`: H.3 first-moment approximation; `2`: H.0/H.1/H.3 two-body kinematics with a drifting Maxwellian ion; `3`: mode 2 plus H.1 `sigma*v` rejection sampling of the ion partner. Mode 3 is the closest current kinetic-parity option, while mode 2 is the validation baseline. |
| `K_NNCs` | `1` | Enable AMMONX neutral-neutral D-D, D-D2, D2-D, and D2-D2 collisions against the fixed neutral background. |
| `K_MarColl` | `1` | Molecular-ion loss representation: `1=separate DS1/DS2/DS3 branches`; `2=legacy summed MAR channel`. Separate branches are required for product-aware D tracking and local D2+ balance. |
| `K_CX_impurity` | `1` | Enable charge exchange for enabled impurity species such as C and Ar. It has no effect while those species switches are zero. |
| `K_Prob` | `1` | `1`: sum local reaction frequencies first and sample from the total; `2`: legacy sum of individual finite-step probabilities. Mode 1 is used for physical collision competition. |
| `K_flight` | `1` | Neutral flight algorithm: `1`: sample a mean free time/path from the local total rate; `2`: fixed-`dt` collision checks; `3`: legacy null-collision/delta tracking using a global minimum mean free time. |
| `dt` | `1e-8 s` | Fixed neutral step for `K_flight=2` and requested charged-particle step when optional D2+ flight is enabled. It does not replace the sampled neutral free-flight time in mode 1. |
| `t_max` | `10` | Legacy maximum-time input. It is parsed but is not currently used as an active termination condition. |
| `K_D2Flight` | `0` | `0`: EIRENE `NFOL=-1` local D2+ production/loss balance, with no D2+ orbit; `1`: experimental fixed-time D2+ test-ion transport and secondary D-product tracking. This does not enable a D2 gas source. |
| `K_EcrossBDrift` | `0` | Add E-cross-B drift to charged test-particle pushes. It does not affect neutral D/D2. |
| `K_abnormal_transport` | `1` | Enable anomalous cross-field diffusion for charged tracked particles. Neutral particles return before this operator, and the switch has no D2+ effect while `K_D2Flight=0`. |
| `K_dn` | `1` | Anomalous diffusivity source: `1`: fixed 0.3 m^2/s; `2`: read the SOLPS radial profile from `Casepath/dn3da.data` in the main chamber. Used only when anomalous charged transport is active. |

### Wall, Target, And Interface Sources

| Setting | Current value | Effect |
| --- | ---: | --- |
| `K_Wallelement` | `1` | Global legacy wall material selector: `1=W`, `2=C`. The current comparison is an all-W wall case. |
| `K_WallRefl` | `1` | Legacy positional compatibility flag (`1=TRIM`, `2=fixed` in the old input convention). The active current reflection backend is selected by `K_Rflect` and, for D-on-W, `K_DWTrimReflection`. |
| `K_Rflect` | `2` | Reflection coefficient source stored as `K_Reflect`: `1`: Eckstein empirical formula; `2`: tabulated TRIM `Rn/RE` files; `3`: alternate generated TRIM tables. |
| `K_RflectDirection` | `3` | Reflected direction model: `1=cosine`; `2=forward/half-cosine`; `3=specular`. The D-on-W distribution model may provide its own sampled direction. |
| `K_backScatter` | `1` | Enable stochastic reflected-atom versus thermal-molecule recycling in the legacy structured-grid wall path. MeshMode 3 performs this branching directly from the reflection and recycling coefficients. |
| `backGridBoundry` | `1` | Legacy structured-grid charged-particle boundary recovery switch. MeshMode 3 uses triangle-boundary and additional-cell logic instead. |
| `coeff_recyc_target` | `1` | Fraction of incident target nuclei returned as recycled neutrals. The reflected atomic part is bounded by this value; the remainder is emitted as molecules. |
| `coeff_ercyc_wall` | `1` | Equivalent recycling coefficient for vessel-wall impacts. The misspelling is part of the positional input format. |
| `coeff_puff` | `0.99` | Local recycling/pumping coefficient at any enabled `K_Puff` location; inactive when all puff switches are zero. |
| `K_DWTrimReflection` | `1` | `1`: use the selected differential D-on-W TRIM distribution for reflection probability, energy, and angle; `0`: use the legacy `K_Rflect` D-W coefficient model. |
| `DWTrimDatabase` | `D_on_W` | D-on-W database directory below `Inputpath/database/`. |
| `DWTrimERMIN_eV` | `1 eV` | Minimum incident energy for fast D-on-W reflection. Below it, fast reflection is disabled and recycled nuclei enter the thermal branch. |
| `SurfaceTemperature_eV` | `0.1 eV` | Wall temperature used for thermal Maxwellian-flux re-emission. |
| `K_DWTargetActualAngle` | `1` | Legacy `DTargetIncidentModel=0` angle: `1`: local B-to-actual-target-normal angle; `0`: fixed 60 degrees. Mode 1 derives both energy and angle from each sampled post-sheath ion velocity. |
| `K_Ei` | `2` | Legacy target incident-energy input for `DTargetIncidentModel=0`: `1`: read `ei_Dion_l/r.data`; `2`: use `2*Ti+3*Te`; `3`: legacy Ti-only path. It does not set the mode-1 sampled ion energy. |
| `DTargetIncidentModel` | `1` | `0`: legacy mean energy/angle; `1`: EIRENE `NEMODS=7` drifting Maxwellian ion flux plus sheath; `2`: normal monoenergetic sensitivity model. |
| `DTargetSheathFactor` | `-1` | For target incident modes, values `<=0` use the EIRENE-style dynamic sheath from local D+ flow; positive values use a fixed multiple of `Te`. |
| `DTargetIncidentSamples` | `4096` | Low-discrepancy samples per target element for incident/reflection source averages. Values below 256 are raised to 256. |
| `K_DTargetSourceMode` | `1` | Target source strength: `1`: read `recycled_neutral_flux_D_l/r.data`; `2`: reconstruct from local `ni*abs(ua)*abs(b.normal)*target_area`. |
| `K_DBoundarySource` | `1` | Add EIRENE interface strata 3-5 from outward D+ `FNIY`. Requires `K_D=1`. |
| `DBoundaryFluxFile` | `2D_data/fnay_Dplus.dat` | Species-resolved face-integrated D+ flux file, relative to `Casepath`. |
| `DBoundaryLaunchModel` | `0` | `0`: EIRENE-style surface recycling into fast D plus thermal D2; `1`: direct outward-D sensitivity source. |
| `numPar_flight_DBoundary` | `6000` | Combined history budget for the two PFR-side and one outer-side interface strata. |

### Variance Reduction, Logs, And Output

| Setting | Current value | Effect |
| --- | ---: | --- |
| `K_Roulette` | `0` | Enable low-weight Russian roulette. It is disabled for the baseline comparison. |
| `W_RouletteMin`, `W_RouletteTarget` | `0.05`, `0.2` | Trigger weight and surviving target weight used only when roulette is enabled. |
| `K_Splitting` | `0` | Enable regional importance splitting. It is disabled for the baseline comparison. |
| `W_SplitMax`, `W_SplitTarget`, `W_SplitMin` | `2`, `1`, `0.05` | Split trigger, desired child weight, and minimum child weight. They are inactive when `K_Splitting=0`. |
| `MaxSplit`, `MaxSplitDepth` | `4`, `1` | Maximum children per split and maximum split generations per source history. |
| `ImportanceOutside`, `ImportanceDivertor`, `ImportanceMainChamber` | `1`, `2`, `1` | Regional importance values used to decide splitting/roulette at region changes. |
| `ImportanceMainPoloidalBegin`, `ImportanceMainPoloidalEnd` | `25`, `72` | B2 poloidal-index interval classified as the main chamber; the remainder is classified with the divertor logic. |
| `K_log` | `0` | Write particle trajectory/velocity logs for supported non-hydrogenic species. Potentially large. |
| `StepLog` | `0` | Enable detailed per-step diagnostic logging. Use only for short debugging runs. |
| `K_H5Output` | `0` | Write HDF5 particle/wall statistics, including D-on-W wall and target data. CSV/text diagnostics remain available when this is off. |

## Physical Model

### Neutral D And D2

Neutral particles are followed with the existing mean-free-path transport
logic.  The local collision set is built from plasma background fields and
database reaction rates.  For the EAST D2 case, the main molecular channels are:

- electron-impact D2 ionization to D2+
- molecular dissociation channels
- D2-D+ charge exchange
- molecular elastic collisions where enabled
- AMMONX D-D, D-D2, D2-D, and D2-D2 collisions against the fixed EIRENE
  neutral background profiles
- target and wall recycling/re-emission

The D2-D+ elastic model follows the three data records requested by the EIRENE
input. AMJUEL H.3 `0.3T` supplies the Maxwell-averaged event frequency, H.1
`0.3T` supplies the relative-energy total cross section and sampled impact
parameter, and H.0 `0.3T` supplies the Morse interaction potential used for the
classical center-of-mass deflection angle. The outgoing relative velocity is
rotated around the incoming relative velocity before the exact two-body
center-of-mass transformation. Heavy-particle arguments are converted from D
or T velocities to the H database mass scale; electron temperature is not
rescaled.

`D2ElasticModel` controls comparison runs:

- `0`: disable molecular ion elastic collisions
- `1`: legacy H.3 `0.3T/0.3D` first-angular-moment approximation
- `2`: H.0/H.1/H.3 kinetics with a drifting Maxwellian ion sample
- `3`: model 2 plus H.1-weighted `sigma*g*f` rejection sampling of the ion
  partner; this is the closest current option to EIRENE conditional partner
  sampling

Run `make -f Makefile.local check-eirene-elastic` to verify the H.0/H.1
database offsets, units, potential minimum, and cross-section extrapolation.

AMMONX collisions are real kinematic events for both atomic and molecular test
particles. The collision samples the matching local atom or molecule background
velocity, updates all three velocity components, and recalculates the test
particle energy after the event. The validation inputs use the full Maxwell
speed distribution rather than the fixed RMS-speed shortcut.

The target model uses the D-on-W TRIM database selected by the EIRENE input.
`DTargetIncidentModel=1` matches the automatically generated EIRENE recycling
source (`NEMODS=7` in input block 14). It samples the local drifting Maxwellian
ion flux directed into the surface, then adds the dynamic sheath acceleration
along the surface normal. The energy and incidence angle derived from that same
post-sheath velocity are passed to the D-on-W TRIM model. The D/D2 source split is a fixed low-discrepancy average,
and fast-D histories sample the corresponding reflection-conditioned incident
distribution. For this input's target normals (`NINCT=-1/+1`), the plasma state
comes from the left/right B2 guard cells. `DTargetIncidentModel=2` is the
normal-monoenergetic sensitivity model using `3*Ti + 0.5*Te + sheath`; mode 0
retains the configured mean-energy and B-to-target-normal angle path. The
optional `K_Ei=1` reader selects the `D@S@+1@N@` section of
`ei_Dion_l/r.data` for mode 0 only. Run
`make -f Makefile.local check-target-incident` to verify the analytic flux
moments, sheath-energy increment, target-surface sampling, and thermal D2
emission moments. Target-recycling histories are sampled over each full
axisymmetric target element with probability proportional to local major
radius, then nudged into the element's mapped plasma triangle. The actual D and
D2 launch position, speed, energy, inward cosine, and mesh-index checks are
written to `target_launch_D.csv` and `target_launch_D2.csv`.
Subsequent real-wall and target re-emissions are audited separately in
`surface_reemission_D.csv` and `surface_reemission_D2.csv`, including the
surface-coordinate residual and the emitted velocity's inward-normal cosine.
The matching `surface_reemission_by_source_D*.csv` files retain the primary
IT, OT, PFR-side, outer-side, or recombination source label through every
later wall cycle. They can therefore be compared directly with EIRENE's
per-stratum `wldra(N)` and `wldrm(N)` surface fluxes.

The optional `K_DBoundarySource=1` path reconstructs EIRENE interface strata
3--5 from the species-resolved SOLPS `fnay_Dplus.dat` field. Negative flux on
the lower radial boundary supplies the two PFR-side strata; positive flux on
the upper radial boundary supplies the outer-side stratum. The file values are
face-integrated D+ rates in s^-1. Mode 1 samples `NEMODS=7` from the matching
radial guard cell and local face orientation. These transparent PFR/outer
surfaces use EIRENE's default D-on-Fe modified Behrisch model (`ILREF=2`,
`EWALL=0.0388 eV`), not the target's D-on-W TRIM model. Fast products launch
cosine-distributed D atoms and the remaining recycled nuclei launch thermal D2
molecules. Source points are sampled on the exact plasma/non-plasma triangle
edge and emitted into the plasma triangle. `D_boundary_source.csv` records both
geometry and guard-cell indices, FNIY strength, incident moments, D/D2 split,
and the source-nuclei closure. `numPar_flight_DBoundary` controls the combined
history budget. Mode 2 retains the former normal-monoenergetic boundary source
as a sensitivity case.

### Transparent Boundaries And Additional Cells

In MeshMode 3, triangle edges tagged `-1` or `-2` are EIRENE transparent
standard surfaces, not reflecting walls. A neutral crossing one of these edges
enters a collision-free additional-cell flight. The code ray-traces to the
nearest exact event: a real vessel-wall segment or a transparent edge through
which the particle re-enters the triangle mesh. No triangle density or plasma
optical depth is accumulated outside the mesh. A particle reflected from a real
wall remains in the additional-cell state until it re-enters the mesh.

The per-species `*_transport_loss_summary.csv` records additional-cell exits,
re-entries, real-wall hits, no-hit losses, H.0 samples, and H.0 numerical
fallbacks. `wall_nearest_fallback` is retained as a regression diagnostic and
should remain zero. MeshMode 3 wall and target incidence angles use the real
surface index, not the triangle-edge number.

### D2+ Local Balance Density

The formal triangular-grid D2+ density is the local balance estimate:

```text
n_D2+ = (P_Ion_D2 + P_CX_D2 + P_CXDT_D2)
        / (ne * (k_2.2.11 + k_2.2.12 + k_2.2.14))
```

In code this is stored in the D2 charged density slot and exported through
`D2p_balance_Tri.csv` and `D2p_balance_B2.csv`.

The D2+ loss channels are AMJUEL H.4:

- `2.2.11`
- `2.2.12`
- `2.2.14`

The local balance output is useful as a collision-rate closure check, but it is
not a finite-orbit residence-time estimator.

### Optional D2+ Test-Ion Transport

`K_D2Flight=1` retains the experimental finite-orbit D2+ implementation for
separate studies, but it is intentionally disabled in the EIRENE validation
inputs because the reference case uses `NFOL=-1`. Its transport principle is:

1. D2 ionization or D2-D+ charge exchange converts the current D2 test particle
   into a charged D2+ test ion.
2. The new D2+ initial charged velocity is obtained by projecting the parent
   neutral D2 velocity onto the local magnetic field.
3. Charged particles are pushed with fixed-time substeps, not neutral
   mean-free-path steps.
4. Each charged substep is limited by the current triangle length scale so that
   a D2+ cannot skip over many mesh cells in one update.
5. Triangle crossing updates `Tri_Index_` and the corresponding B2 cell index
   before the next force/collision evaluation.
6. D2+ undergoes DS1/DS2/DS3 loss sampling after finite residence time.
7. Leaving the valid plasma triangle mesh, reaching non-finite state, or
   exceeding the numerical speed guard terminates the current D2+ as a boundary
   or numerical loss.

The charged push includes parallel ion friction/thermal relaxation and optional
electric-field/drift terms when the corresponding switches are enabled.

### Speed Guard

`D2pMaxAllowedSpeed` is a numerical sanity guard, currently `1.0e6 m/s`.  It is
not intended as a physical clamp in normal simulations.  For D2+ this is far
above the expected edge-ion velocity scale; normal speeds in the EAST tests are
well below the guard.  If the guard triggers, inspect mesh indexing, field
lookup, and charged substep logic.

## Important Output Files

For D2+ studies, check these files first:

```text
D2p_transport_status.txt
D2p_mesh3_flight_audit.csv
D2p_tri_density_interpretation.txt
D2p_balance_Tri.csv
D2p_balance_B2.csv
D2p_track_length_Tri.csv
D2p_track_length_summary.csv
D2p_mesh3_tracking_fate_summary.csv
D2p_initial_velocity_audit.csv
D2p_production_capability_Tri.csv
D2p_production_capability_summary.csv
D2p_physics_decomposition_B2.csv
D2p_physics_decomposition_summary.csv
```

`D2p_track_length_Tri.csv` contains the finite-orbit track-length density
diagnostic.  Compare it against the local-balance density to separate transport
effects from collision-rate closure.

Useful summary quantities:

- total D2 inventory
- `P_Ion_D2`
- `P_CX_D2`
- D2+ local-balance inventory
- D2+ track-length inventory
- DS1/DS2/DS3 fate weights
- boundary-loss weight
- maximum D2+ speed relative to the guard

## SOLPS/EIRENE Comparison

The comparison workflow is:

1. Run NT on one of the EAST triangular settings.
2. Use `triangle/EAST/compare_tri_programs.py` or
   `triangle/EAST/benchmark_tri.py` to map NT and SOLPS/EIRENE quantities onto
   the same triangle grid.
3. For D2+, compare both:
   - local-balance density
   - track-length density
4. Also compare D2 density, because D2+ depends strongly on the parent D2
   inventory and spatial distribution.

`benchmark_tri.py` reports both the direct `dab2/dmb2` B2-field comparison and
the sum of the six primary-source track-length inventories against EIRENE's
`pdena_int_b2(0)`/`pdenm_int_b2(0)`. It also reports a reference
self-consistency row for `d?b2 * vol_2D.dat` versus those EIRENE integrals;
do not attribute that residual to NT. Use `check_b2_alignment.py` when a narrow
target region disagrees strongly, to rule out an axis flip or integer-cell
offset before changing transport physics.

For the 5e19 no-pump/no-puff test, D2+ disagreement should be debugged as a
chain:

```text
recycling source -> D2 density -> D2 ionization/CX -> D2+ loss and residence
```

Do not tune D2+ transport in isolation before confirming the D2 source and D2
density are consistent.

### Reaction-Coefficient Export

Build and run the B2-grid exporter with an explicit physical projectile energy
for the H.3 diagnostic slice:

```bash
make -f Makefile.local export-eirene-reactions
./bin/export_eirene_reactions Inputfile/settingfile/setting_Trimesh_D_5.log 1.0
```

An optional third argument overrides `Outputpath`. Use it for density scans so
each case keeps its own coefficient maps:

```bash
./bin/export_eirene_reactions \
  Inputfile/settingfile/setting_Trimesh_D_1.log 1.0 \
  Outputfile/10292310/data/reaction_coefficients
python3 triangle/EAST/check_collision_closure.py \
  --case case_input/2MW-3.2e19 \
  --output Outputfile/10292310 \
  --rates Outputfile/10292310/data/reaction_coefficients/eirene_reaction_b2
```

The exporter reads `ne`, `Te`, `Ti`, atom temperature, and molecule temperature
from the setting's `Casepath`. Electron reactions use `ne,Te`; H.2 3.2.3 uses
local `Ti`; AMMONX reactions use the matching neutral temperature; hydrogen H.3
uses local `Ti` and the physical projectile energy supplied as the second
argument. D/T arguments are automatically mapped to the H-database mass scale.
The generated `value_m3_s` is generally a rate coefficient `<sigma v>`, not a
raw cross section. During particle flight, H.3 is still evaluated from each
test particle's instantaneous energy; the exported fixed-energy map is only a
controlled comparison slice.

## Program Implementation Notes

The high-level flow is:

1. `Initialize.cpp` reads the setting file.
2. `Read.cpp` loads SOLPS/EIRENE case data and plasma fields.
3. `Prepare.cpp` builds collision rates and mesh-indexed coefficients.
4. `Moncar.cpp` creates source particles and drives histories.
5. `Particle.cpp` handles particle state, collision sampling, neutral tracking,
   D2+ creation, charged pushing, and diagnostic accumulation.
6. `Output.cpp` writes densities, source terms, D2+ audits, and comparison
   diagnostics.

Mesh mode 3 uses the external triangular mesh.  When editing charged transport,
keep the triangle index and B2 index synchronized; many D2+ failures come from
using a stale cell index after crossing an edge.

## Optimization And Numerical Stability

Recent D2+ transport changes are intentionally conservative:

- Charged D2+ substeps are triangle-limited before force updates.
- Parallel ion relaxation uses an exact exponential relaxation form rather than
  explicit Euler for stiff friction terms.
- Invalid triangle/cell transitions are treated as boundary/numerical loss
  instead of continuing with undefined plasma fields.
- Target hits with a stale or missing B2 radial index are resolved from the
  current triangle and nearest target-segment midpoint only at the boundary
  event; target geometry is not searched on every transport step.
- A neutral step whose endpoint remains inside the current triangle follows the
  normal collision-probability path instead of being counted as a lost particle.
- Diagnostic summaries report maximum speed, low-speed residence fractions,
  fate weights, and source channel weights.

When optimizing further:

- preserve local-balance closure before changing transport;
- keep neutral mean-free-path logic separate from charged fixed-time stepping;
- avoid silently reusing stale `Tri_Index_` or `XY_`;
- add diagnostics before changing collision tables or source normalization;
- compare D2 first, then D2+.

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
- `2`: EIRENE-default H.0/H.1/H.3 kinetics with a drifting Maxwellian ion sample
- `3`: model 2 plus H.1-weighted `sigma*g*f` rejection sampling of the ion partner

Run `make -f Makefile.local check-eirene-elastic` to verify the H.0/H.1
database offsets, units, potential minimum, and cross-section extrapolation.

AMMONX collisions are real kinematic events for both atomic and molecular test
particles. The collision samples the matching local atom or molecule background
velocity, updates all three velocity components, and recalculates the test
particle energy after the event. The validation inputs use the full Maxwell
speed distribution rather than the fixed RMS-speed shortcut.

The target model uses the D-on-W TRIM database selected by the EIRENE input.
`DTargetIncidentModel=1` matches the automatically generated EIRENE recycling
source (`NEMODS=-3`): every source ion is monoenergetic and normally incident.
For the main target-face flux, the incident energy is reconstructed from
EIRENE's `ELSTEP/FLSTEP` expression, `3*Ti + 0.5*Te`, and the dynamic sheath
energy is then added. For this input's target normals (`NINCT=-1/+1`),
`INFCOP` obtains those plasma values from the left/right B2 guard cells. The
optional `K_Ei=1`
reader still selects the `D@S@+1@N@` section of `ei_Dion_l/r.data` for legacy
mode 0, but that diagnostic is not used by mode 1. `DTargetIncidentModel=2`
retains the drifting Maxwellian-flux plus sheath model as a sensitivity case,
while mode 0 retains the configured mean-energy and B-to-target-normal angle
path. Run
`make -f Makefile.local check-target-incident` to verify the analytic flux
moments and sheath-energy increment.

The optional `K_DBoundarySource=1` path reconstructs EIRENE interface strata
3--5 from the species-resolved SOLPS `fnay_Dplus.dat` field. Negative flux on
the lower radial boundary supplies the two PFR-side strata; positive flux on
the upper radial boundary supplies the outer-side stratum. The file values are
face-integrated D+ rates in s^-1. For EIRENE mode 1, the source state comes from
the radial guard cell, has normal incidence, and uses the radial-boundary energy
`2*Ti + sheath`; the optional inclined-face flux term is not present in the
available export. These transparent PFR/outer surfaces use EIRENE's default
D-on-Fe modified Behrisch model (`ILREF=2`, `EWALL=0.0388 eV`), not the target's
D-on-W TRIM model. Fast products launch cosine-distributed D atoms and the
remaining recycled nuclei launch thermal D2 molecules. Source points are
sampled on the exact plasma/non-plasma triangle edge and emitted into the
plasma triangle. `D_boundary_source.csv` records both geometry and guard-cell
indices, FNIY strength, incident moments, D/D2 split, and the source-nuclei
closure. `numPar_flight_DBoundary` controls the combined history budget.

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

For the 5e19 no-pump/no-puff test, the key diagnostic lesson was that D2 itself
was already low relative to EIRENE, and D2+/D2 was also lower.  Therefore D2+
disagreement should be debugged as a chain:

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

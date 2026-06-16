# Unused Code Candidates

Generated: 2026-06-16
Scope: src/ and include/

This is a static-scan candidate list. Before deleting, remove in small batches and run make -j2.
Backup files such as *.bak_memopt_20260615 were treated as noise.

## High-Confidence Unused Functions

These appear to have only declaration + definition, with no real call sites.

| Symbol | Definition | Declaration |
|---|---|---|
| ADAS::waveLength(int) | src/ADAS.cpp:193 | include/ADAS.h:45 |
| GRID::Rectangular_Step(int) | src/Grid.cpp:710 | include/Grid.h:98 |
| GRID::Target1_Index(int) | src/Grid.cpp:733 | include/Grid.h:101 |
| GRID::Target2_Index(int) | src/Grid.cpp:738 | include/Grid.h:102 |
| andomSign() | src/utils.cpp:119 | include/utils.h:43 |

## Likely Unused Functions, Confirm Before Deleting

These look unused in the current main flow, but may be old debug/output interfaces or kept for manual experiments.

| Symbol | Definition / note |
|---|---|
| Particle::RecombinCal(...) | src/Particle.cpp:224, include/Particle.h:490 |
| Particle::lambda_now() | src/Particle.cpp:1206, include/Particle.h:558 |
| Particle::Setdt(double) | src/Particle.cpp:6708, include/Particle.h:509 |
| Particle::FlightTime() | src/Particle.cpp:6713, include/Particle.h:510 |
| Particle::CoreTotalcs(int) | src/Particle.cpp:6779, include/Particle.h:591 |
| Particle::SetCharge(int) | src/Particle.cpp:7358, include/Particle.h:584 |
| Particle::T_n(int,int,int) | src/Particle.cpp:7425, include/Particle.h:562 |
| Particle::Weight_Target() | src/Particle.cpp:7435, include/Particle.h:564 |
| Particle::SetTn(double) | src/Particle.cpp:7453, include/Particle.h:513 |
| Particle::setname(string) | src/Particle.cpp:7466, include/Particle.h:580 |
| Particle::SetChargeTag(int) | src/Particle.cpp:7468, include/Particle.h:578 |
| Particle::setPar(string,double,int) | src/Particle.cpp:7470, include/Particle.h:581 |
| Particle::regression() | src/Particle.cpp:8395, include/Particle.h:508 |
| Particle::setmass(double) | src/Particle.cpp:8895, include/Particle.h:582 |
| ParCollCar::EIRENEInput(EIRENE*) | src/Particle.cpp:9262, include/Particle.h:144 |
| ParCollCar::Set_K(int) | src/Particle.cpp:9438, include/Particle.h:89; implementation currently looks suspicious: K_ - i; |
| ParCollCar::Mu_Add1(int[],double) | src/Particle.cpp:9765, include/Particle.h:124 |
| ParCollCar::Mu_Add2(int[],double) | src/Particle.cpp:9775, include/Particle.h:125 |
| Xrecycling(double,int,int) | src/Particle.cpp:9876, include/Global.h:59 |
| CX_DT_Fix() | src/Prepare.cpp:1515; call in Prepare() is commented out near src/Prepare.cpp:15 |

## High-Confidence Unused Globals

These appear only as definitions, or otherwise have no real use sites.

| Symbol | Location |
|---|---|
| K_CD4 | src/Global.cpp:69 |
| 
_H_1_Tri | src/Global.cpp:117 |
| 
_D_1_Tri | src/Global.cpp:117 |
| 
_T_1_Tri | src/Global.cpp:117 |
| HHmass | src/Global.cpp:166 |

## Likely Unused Globals, Confirm Before Deleting

Usually declaration + definition, but no observed read/write in current code.

| Symbol | Location / note |
|---|---|
| K_PartoPar | include/Global.h:66, src/Global.cpp:73 |
| Te_vacuum | include/Global.h:82, src/Global.cpp:87 |
| 
e_vacuum | include/Global.h:82, src/Global.cpp:87 |
| nax_e | include/Global.h:149, src/Global.cpp:233 |
| Ar_W | include/Global.h:243, src/Global.cpp:319 |
| He_W | include/Global.h:244, src/Global.cpp:320 |
| R2_1_7_H2 | check nearby H2 reaction declarations/definitions |
| R7_2_3a_H4 | include/Global.h:287 |
| R7_2_3b_H4 | include/Global.h:288 |
| RM_2_23 | include/Global.h:305 |
| RM_3_2 | include/Global.h:310 |

## Likely Unused Member Variables

These were found with very low reference counts, usually declaration plus initialization/assignment only.

| Symbol | Location / note |
|---|---|
| ParCollCar::cc | include/Particle.h:30 |
| ParCollCar::Tri_Smu_0_ | include/Particle.h:64 |
| ParCollCar::Tri_Smu_1_ | include/Particle.h:65 |
| ParCollCar::Tri_Smu_2_ | include/Particle.h:66 |
| ParCollCar::Tri_Smu1_ | include/Particle.h:67 |
| ParCollCar::Tri_Smu2_ | include/Particle.h:68 |
| ParCollCar::Tri_Num_mu1_ | include/Particle.h:75 |
| ParCollCar::Tri_Num_mu2_ | include/Particle.h:76 |
| ParCollCar::Tri_crossSection | include/Particle.h:79 |
| Particle::T_array_ | include/Particle.h:403 |
| PartoPar::ifChange_ | include/Particle.h:596 |
| Grid_extern::dnx_ | include/Grid.h:119; assigned in src/Grid.cpp:947 |
| Grid_extern::dny_ | include/Grid.h:120; assigned in src/Grid.cpp:948 |
| Grid_extern::Nbc_ | include/Grid.h:121; assigned in src/Grid.cpp:1008 |
| Grid_extern::Nwall_ | include/Grid.h:124; assigned in src/Grid.cpp:1107 |
| Grid_extern::cvX_ | include/Grid.h:131; read from file in src/Grid.cpp:845 |
| Grid_extern::cvY_ | include/Grid.h:132; read from file in src/Grid.cpp:846 |
| Grid_extern::cvVol_ | include/Grid.h:142; read from file in src/Grid.cpp:855 |
| Grid_extern::cvReg_ | include/Grid.h:143; read from file in src/Grid.cpp:856 |
| Grid_extern::erx_e_ | include/Grid.h:151; filled in src/Grid.cpp:860 |
| Grid_extern::ery_e_ | include/Grid.h:152; filled in src/Grid.cpp:861 |

## Suggested First Deletion Batch

Safest first batch based on the scan:

- ADAS::waveLength(int)
- GRID::Rectangular_Step(int)
- GRID::Target1_Index(int)
- GRID::Target2_Index(int)
- andomSign()
- K_CD4
- 
_H_1_Tri, 
_D_1_Tri, 
_T_1_Tri
- HHmass
- ParCollCar::cc

After deleting each batch:

`ash
make -j2
`

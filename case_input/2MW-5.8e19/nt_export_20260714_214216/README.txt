SOLPS-ITER / EIRENE comparison export
======================================
Case directory: /home/bianyu/SOLPS_3.0.8/solps-iter/runs/EAST_62421/2MW-5.8e19-new-test2-wo-n-n
Created: 2026-07-14T21:42:19+0800

Important files
---------------
parsed/case_metadata.txt
  Parsed NFILE, NGSTAL and NSTRAI.
parsed/availability_report.txt
  States what this run can and cannot provide, especially individual ISTRA data.
parsed/reaction_map.tsv
  EIRENE internal reaction number to AMJUEL/HYDHEL label.
parsed/custom_collision_tallies.tsv
  Existing COLV reaction tallies from input block 10b.
parsed/interface_strata_raw.tsv
  Raw interface records for ISTRA 1..5 from input block 14.
b2plot/fnax_all_species.dat and b2plot/fnay_all_species.dat
  Face-integrated B2 particle flows for every B2 species.
parsed/eirene_tally_text_hits.txt
  Grep extract of FLUX/ISTRA/density/source/surface/energy tally lines.
raw/
  Unmodified input and raw output files for reproducibility.

D+ zsel export
---------------
No B2 species index is guessed by this script. Set it explicitly, e.g.:
  DPLUS_INDEX=0 bash export_solps_eirene_compare.sh
The EIRENE bulk-ion index '1 D+' in input.dat is NOT proof that B2PLOT zsel is 1.

Per-ISTRA limitation
--------------------
Read parsed/availability_report.txt. If NFILE does not end in 1 or 2, total
fort.44/fort.46 cannot be decomposed after the run into PDENA/PDENM/PPAT/PPML
for ISTRA 1..5. A new EIRENE run must save individual-stratum FT10/FT11 data.

#!/bin/bash
set -eu

# Submit a matched OFF/ON pair. Repeat this pair at least five times before
# comparing median Moncar time and statistical uncertainty.
sbatch crslum.slurm Inputfile/settingfile/setting_Trimesh_D_5_vr_off.log
sbatch crslum.slurm Inputfile/settingfile/setting_Trimesh_D_5.log

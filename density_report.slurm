#!/bin/bash
#SBATCH -M mindwell
#SBATCH -p batch_graniterapids
#SBATCH -N 1
#SBATCH -J NT_density_report
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=log/DensityReport-%j.out
#SBATCH --error=log/DensityReport-%j.err
#SBATCH -t 01:00:00
#SBATCH -A lp_opt_plasma_edge_conf

set -e

python_executable=${NT_PYTHON:-/data/leuven/379/vsc37950/opt/miniconda/envs/dbg/bin/python3}
output_dir=${10:-Outputfile/density_scan_report}

mkdir -p log "${output_dir}"

"${python_executable}" triangle/EAST/build_density_scan_report.py \
  --case "3.2e19=Outputfile/${1}/fig/tri_benchmark" \
  --case "3.6e19=Outputfile/${2}/fig/tri_benchmark" \
  --case "4.0e19=Outputfile/${3}/fig/tri_benchmark" \
  --case "4.6e19=Outputfile/${4}/fig/tri_benchmark" \
  --case "5.0e19=Outputfile/${5}/fig/tri_benchmark" \
  --case "5.4e19=Outputfile/${6}/fig/tri_benchmark" \
  --case "5.5e19=Outputfile/${7}/fig/tri_benchmark" \
  --case "5.6e19=Outputfile/${8}/fig/tri_benchmark" \
  --case "5.8e19=Outputfile/${9}/fig/tri_benchmark" \
  --summary-output "${output_dir}/density_scan_summary.pdf" \
  --full-output "${output_dir}/eirene_density_scan_full.pdf"

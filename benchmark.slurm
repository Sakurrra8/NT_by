#!/bin/bash
#SBATCH -M mindwell
#SBATCH -p batch_graniterapids
#SBATCH -N 1
#SBATCH -J NT_benchmark
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=log/Benchmark-%j.out
#SBATCH --error=log/Benchmark-%j.err
#SBATCH -t 00:30:00
#SBATCH -A lp_opt_plasma_edge_conf

set -e

case_path=$1
simulation_job=$2
case_label=$3
setting_file=$4
benchmark_dir="Outputfile/${simulation_job}/fig/tri_benchmark"
python_executable=${NT_PYTHON:-/data/leuven/379/vsc37950/opt/miniconda/envs/dbg/bin/python3}

mkdir -p log "${benchmark_dir}" "Outputfile/${simulation_job}/.matplotlib"
export MPLCONFIGDIR="$PWD/Outputfile/${simulation_job}/.matplotlib"

"${python_executable}" triangle/EAST/benchmark_tri.py \
  --case "${case_path}" \
  --output "Outputfile/${simulation_job}" \
  --setting "${setting_file}"

"${python_executable}" triangle/EAST/assemble_eirene_case_report.py \
  --benchmark-dir "${benchmark_dir}" \
  --output "${benchmark_dir}/eirene_comparison_full.pdf" \
  --case-label "${case_label}" \
  --job-id "${simulation_job}"

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
eirene_output=$5
eirene_input=$6
benchmark_dir="Outputfile/${simulation_job}/fig/tri_benchmark"
if [[ -n "${eirene_output}" ]]; then
  benchmark_dir="Outputfile/${simulation_job}/fig/tri_benchmark_$(basename "$(dirname "${eirene_output}")")"
fi
python_executable=${NT_PYTHON:-/data/leuven/379/vsc37950/opt/miniconda/envs/dbg/bin/python3}

mkdir -p log "${benchmark_dir}" "Outputfile/${simulation_job}/.matplotlib"
export MPLCONFIGDIR="$PWD/Outputfile/${simulation_job}/.matplotlib"

benchmark_args=(
  --case "${case_path}"
  --output "Outputfile/${simulation_job}"
  --setting "${setting_file}"
  --outdir "${benchmark_dir}"
)
if [[ -n "${eirene_output}" ]]; then
  benchmark_args+=(--eirene-output "${eirene_output}")
fi
if [[ -n "${eirene_input}" ]]; then
  benchmark_args+=(--eirene-input "${eirene_input}")
fi
"${python_executable}" triangle/EAST/benchmark_tri.py "${benchmark_args[@]}"

"${python_executable}" triangle/EAST/assemble_eirene_case_report.py \
  --benchmark-dir "${benchmark_dir}" \
  --output "${benchmark_dir}/eirene_comparison_full.pdf" \
  --case-label "${case_label}" \
  --job-id "${simulation_job}"

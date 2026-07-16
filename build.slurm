#!/bin/bash
#SBATCH -M mindwell
#SBATCH -p batch_graniterapids
#SBATCH -N 1
#SBATCH -J NT_build
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=log/Build-%j.out
#SBATCH --error=log/Build-%j.err
#SBATCH -t 00:20:00
#SBATCH -A lp_opt_plasma_edge_conf

set -e
mkdir -p log obj bin
make -j "${SLURM_CPUS_PER_TASK:-1}"

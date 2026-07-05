#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

make -f Makefile.local BUILD=release -j1

for case_id in {1..9}; do
  setting="Inputfile/settingfile/setting_wallscan_recyc500k_D_${case_id}.log"
  outdir="Outputfile/wallscan_recyc500k_D_${case_id}"
  mkdir -p "${outdir}/data" "${outdir}/log"
  echo "=== wallscan_recyc500k_D_${case_id}: ${setting} ==="
  ./bin/edit "${setting}" >"${outdir}/log/run.out" 2>"${outdir}/log/run.err"
done

python3 triangle/EAST/plot_wallscan_erosion.py \
  --root . \
  --run-prefix wallscan_recyc500k_D \
  --summary-name wallscan_recyc500k_summary

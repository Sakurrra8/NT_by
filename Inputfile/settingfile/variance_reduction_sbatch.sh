#!/bin/bash
set -euo pipefail

# Usage:
#   bash Inputfile/settingfile/variance_reduction_sbatch.sh pair
#   bash Inputfile/settingfile/variance_reduction_sbatch.sh on 5
#   bash Inputfile/settingfile/variance_reduction_sbatch.sh off 1 3 5
#
# mode:
#   off  submit original setting_Trimesh_D_N.log files (VR disabled)
#   on   submit setting_Trimesh_D_N_vr_on.log files (VR enabled)
#   pair submit matched OFF and ON jobs for direct comparison
#
# If no case numbers are supplied, cases 1 through 9 are submitted.

mode=${1:-pair}
if (($# > 0)); then
    shift
fi

if (($# > 0)); then
    cases=("$@")
else
    cases=(1 2 3 4 5 6 7 8 9)
fi

submit_case()
{
    local case_number=$1
    local variant=$2
    local setting

    if [[ $variant == "off" ]]; then
        setting="Inputfile/settingfile/setting_Trimesh_D_${case_number}.log"
    else
        setting="Inputfile/settingfile/setting_Trimesh_D_${case_number}_vr_on.log"
    fi

    [[ -f $setting ]] || {
        echo "Missing setting file: $setting" >&2
        return 1
    }

    echo "Submitting case ${case_number}, variance reduction ${variant}: ${setting}"
    sbatch crslum.slurm "$setting"
}

mkdir -p log Outputfile

case "$mode" in
    off|on)
        for case_number in "${cases[@]}"; do
            submit_case "$case_number" "$mode"
        done
        ;;
    pair)
        for case_number in "${cases[@]}"; do
            submit_case "$case_number" off
            submit_case "$case_number" on
        done
        ;;
    *)
        echo "Unknown mode '$mode'. Use: off, on, or pair." >&2
        exit 2
        ;;
esac

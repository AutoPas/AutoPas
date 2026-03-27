#!/bin/bash
# Local dry-check for rerun status.
# Mirrors the iteration logic of job_RerunFailed.sh without Slurm and without executing srun.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_INPUT_DIR="${SCRIPT_DIR}/generated_inputs"

# Mirror the Slurm array range from job_RerunFailed.sh by default.
START_INDEX="${START_INDEX:-40}"
END_INDEX="${END_INDEX:-49}"

if [[ ! "${START_INDEX}" =~ ^[0-9]+$ ]] || [[ ! "${END_INDEX}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: START_INDEX and END_INDEX must be non-negative integers." >&2
    exit 2
fi

if [ "${START_INDEX}" -gt "${END_INDEX}" ]; then
    echo "ERROR: START_INDEX (${START_INDEX}) must be <= END_INDEX (${END_INDEX})." >&2
    exit 2
fi

if [ ! -d "${BASE_INPUT_DIR}" ]; then
    echo "ERROR: Input base directory not found: ${BASE_INPUT_DIR}" >&2
    exit 2
fi

# Generate arrays of parameters indexed by job ID.
declare -a container
declare -a traversal
declare -a sigma_ratio
declare -a cell_size
index=0

for container_iter in HierarchicalGridMatching HierarchicalGrid LinkedCells
 do
    case "${container_iter}" in
        HierarchicalGridMatching)
            traversals=(hgrid_matching)
            ;;
        HierarchicalGrid)
            traversals=(hgrid_block4 hgrid_block8)
            ;;
        LinkedCells)
            traversals=(lc_c08 lc_c04_HCP)
            ;;
        *)
            echo "ERROR: Unknown container '${container_iter}'" >&2
            exit 2
            ;;
    esac

    for traversal_iter in "${traversals[@]}"
    do
        for sigma_ratio_iter in 0p10 0p20 0p30 0p40 0p48
        do
            for cell_size_iter in 0p50 1p00
            do
                container[$index]="${container_iter}"
                traversal[$index]="${traversal_iter}"
                sigma_ratio[$index]="${sigma_ratio_iter}"
                cell_size[$index]="${cell_size_iter}"
                index=$((index + 1))
            done
        done
    done
 done

max_index=$((index - 1))
if [ "${START_INDEX}" -gt "${max_index}" ] || [ "${END_INDEX}" -gt "${max_index}" ]; then
    echo "ERROR: Requested index range ${START_INDEX}-${END_INDEX} exceeds valid range 0-${max_index}." >&2
    exit 2
fi

echo "Checking indices ${START_INDEX}-${END_INDEX} (valid: 0-${max_index})"
echo "Base input dir: ${BASE_INPUT_DIR}"

skipped=0
would_run=0
missing_run_dirs=0

for task_id in $(seq "${START_INDEX}" "${END_INDEX}")
do
    target_base_dir="${BASE_INPUT_DIR}/container_${container[${task_id}]}/traversal_${traversal[${task_id}]}/sigmaRatio_${sigma_ratio[${task_id}]}/cellSize_${cell_size[${task_id}]}"

    if [ ! -d "${target_base_dir}" ]; then
        echo "[MISSING BASE] task=${task_id} path=${target_base_dir}"
        continue
    fi

    for count_ratio_iter in 0p50 1p00 2p00 4p00 8p00 16p00
    do
        count_dir="${target_base_dir}/countRatio_${count_ratio_iter}"
        if [ ! -d "${count_dir}" ]; then
            echo "[MISSING COUNT] task=${task_id} path=${count_dir}"
            continue
        fi

        for run in $(seq 0 2)
        do
            run_dir="${count_dir}/run_${run}"
            if [ ! -d "${run_dir}" ]; then
                echo "[MISSING RUN] task=${task_id} path=${run_dir}"
                missing_run_dirs=$((missing_run_dirs + 1))
                continue
            fi

            if compgen -G "${run_dir}/MDFlex_end*" > /dev/null; then
                echo "[SKIP] task=${task_id} run_dir=${run_dir}"
                skipped=$((skipped + 1))
            else
                echo "[WOULD RUN] task=${task_id} run_dir=${run_dir}"
                would_run=$((would_run + 1))
            fi
        done
    done
done

echo ""
echo "Summary:"
echo "  Skipped (MDFlex_end* exists): ${skipped}"
echo "  Would run (MDFlex_end* missing): ${would_run}"
echo "  Missing run directories: ${missing_run_dirs}"

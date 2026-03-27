#!/bin/bash
#SBATCH --job-name=Uni_2_single_10_17
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --cpus-per-task=56
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
#SBATCH --mail-user=alexander.glas@tum.de
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load slurm_setup

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
MD_FLEXIBLE_BIN="${MD_FLEXIBLE_BIN:-${SCRIPT_DIR}/../../build/examples/md-flexible/md-flexible}"

if [ ! -x "${MD_FLEXIBLE_BIN}" ]; then
    echo "ERROR: md-flexible executable not found or not executable: ${MD_FLEXIBLE_BIN}" >&2
    echo "Hint: export MD_FLEXIBLE_BIN=/absolute/path/to/md-flexible before sbatch." >&2
    exit 2
fi

RUN_DIR="${SCRIPT_DIR}/generated_inputs/container_HierarchicalGrid/traversal_hgrid_block4/sigmaRatio_0p10/cellSize_0p50/countRatio_16p00/run_2"
if [ ! -d "${RUN_DIR}" ]; then
    echo "ERROR: Run directory not found: ${RUN_DIR}" >&2
    exit 2
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=true

unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

cd "${RUN_DIR}"
srun -n 1 "${MD_FLEXIBLE_BIN}" --yaml-filename input.yaml > "logOutput_${SLURM_JOB_ID}_10.out" 2>&1

echo "Completed single run successfully."

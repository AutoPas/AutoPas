#!/bin/bash
#SBATCH --job-name=Uni_2_Job_Array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny  # ! This probably needs to change !
#SBATCH --cpus-per-task=56 # ! This should be the maximum number of CPUs you wish to use per single MPI rank !
#SBATCH --time=03:00:00
#SBATCH --output=/dev/null
# 3 containers with traversals: 1 + 2 + 2 = 5 traversal combinations
# Index only container/traversal, sigma ratio and cell size.
# The count-ratio sweep and 3 repeat runs are executed inside each array job.
# 5 (container/traversal) * 5 sigma ratios * 2 cell sizes = 50 jobs
#SBATCH --array=49-49
#SBATCH --mail-user=alexander.glas@tum.de
#SBATCH --mail-type=END,FAIL  # Send email on end and failure

# !! Modify the above as relevant for your computer !!
# !! The ones most likely needing changes are highlighted, but check all of them !!

echo "#==================================================#"
echo " num nodes: " $SLURM_JOB_NUM_NODES
echo " num tasks: " $SLURM_NTASKS
echo " cpus per task: " $SLURM_CPUS_PER_TASK
echo " nodes used: " $SLURM_JOB_NODELIST
echo " job cpus used: " $SLURM_JOB_CPUS_PER_NODE
echo "#==================================================#"

# Change this to whatever is needed to run an MPI simulation. Can be removed if not running MPI simulations.
# The data collection only needs one rank, so this is not strictly necessary even if running MPI simulations
# with the trained model.
#module load openmp
module load slurm_setup

# Below is a common method for assigning the SLURM job array's ID to a particular directory, as created by
# the input_generator.py.

# In batch mode, SLURM_SUBMIT_DIR points to the directory from which sbatch was called.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
MD_FLEXIBLE_BIN="${MD_FLEXIBLE_BIN:-${SCRIPT_DIR}/../../build/examples/md-flexible/md-flexible}"

if [ ! -x "${MD_FLEXIBLE_BIN}" ]; then
    echo "ERROR: md-flexible executable not found or not executable: ${MD_FLEXIBLE_BIN}" >&2
    echo "Hint: export MD_FLEXIBLE_BIN=/absolute/path/to/md-flexible before sbatch." >&2
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
                container[$index]="$container_iter"
                traversal[$index]="$traversal_iter"
                sigma_ratio[$index]="$sigma_ratio_iter"
                cell_size[$index]="$cell_size_iter"
                index=$((index + 1))
            done
        done
    done
done

if [ "${SLURM_ARRAY_TASK_ID}" -ge "${index}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range [0,$((index - 1))]." >&2
    exit 2
fi
     
# Go to base directory of experiment for this job ID.
TARGET_BASE_DIR="${SCRIPT_DIR}/generated_inputs/container_${container[${SLURM_ARRAY_TASK_ID}]}/traversal_${traversal[${SLURM_ARRAY_TASK_ID}]}/sigmaRatio_${sigma_ratio[${SLURM_ARRAY_TASK_ID}]}/cellSize_${cell_size[${SLURM_ARRAY_TASK_ID}]}"
if [ ! -d "${TARGET_BASE_DIR}" ]; then
    echo "ERROR: Input base directory not found: ${TARGET_BASE_DIR}" >&2
    echo "Hint: Generate inputs with input_generatorCompare2.py before submitting." >&2
    exit 2
fi

# Set up environment variables
# Check your machine to see what is recommended here
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=true

unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

# Sequentially loop over innermost sweep (count ratio) and over each repeat run.
for count_ratio_iter in 0p50 1p00 2p00 4p00 8p00 16p00
do
    COUNT_DIR="${TARGET_BASE_DIR}/countRatio_${count_ratio_iter}"
    if [ ! -d "${COUNT_DIR}" ]; then
        echo "ERROR: Input count-ratio directory not found: ${COUNT_DIR}" >&2
        exit 2
    fi

    cd "${COUNT_DIR}"
    for run in $(seq 0 2)
    do
        cd "run_${run}"

        # Skip completed runs that already produced MDFlex_end*.
        if compgen -G "MDFlex_end*" > /dev/null; then
            echo "Skipping ${COUNT_DIR}/run_${run}: found MDFlex_end*."
            cd ..
            continue
        fi

        # Clean stale output logs before rerun.
        rm -f -- *.out

        # If not running an mpi executable, change srun to whatever is recommended for your machine.
        # Override MD_FLEXIBLE_BIN if your executable is in a different location.
        srun -n 1 "${MD_FLEXIBLE_BIN}" --yaml-filename input.yaml > logOutput_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1
        cd ..
    done
    cd "${TARGET_BASE_DIR}"
done

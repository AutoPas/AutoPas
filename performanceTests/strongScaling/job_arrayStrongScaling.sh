#!/bin/bash
#SBATCH --job-name=Uni_2_Job_Array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny  # ! This probably needs to change !
#SBATCH --cpus-per-task=17 # ! ${in 1 2 4 8 14 28 56 112}
#SBATCH --hint=nomultithread
#SBATCH --time=03:00:00
#SBATCH --output=/dev/null
# 3 containers with traversals: 1 + 2 + 1 = 4 traversal combinations
# Index only container/traversal and cell size.
# The thread sweep and 3 repeat runs are executed inside each array job.
# 4 (container/traversal) * 2 cell sizes = 8 jobs
#SBATCH --array=2-5 # 0-7
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
            traversals=(lc_c08)
            ;;
        *)
            echo "ERROR: Unknown container '${container_iter}'" >&2
            exit 2
            ;;
    esac

    for traversal_iter in "${traversals[@]}"
    do
        for cell_size_iter in 0p50 1p00
        do
            container[$index]="$container_iter"
            traversal[$index]="$traversal_iter"
            cell_size[$index]="$cell_size_iter"
            index=$((index + 1))
        done
    done
done

if [ "${SLURM_ARRAY_TASK_ID}" -ge "${index}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range [0,$((index - 1))]." >&2
    exit 2
fi
     
# Resolve base directory of experiment for this job ID.
TRAVERSAL_DIR="${SCRIPT_DIR}/generated_inputs/container_${container[${SLURM_ARRAY_TASK_ID}]}/traversal_${traversal[${SLURM_ARRAY_TASK_ID}]}"
if [ ! -d "${TRAVERSAL_DIR}" ]; then
    echo "ERROR: Input traversal directory not found: ${TRAVERSAL_DIR}" >&2
    echo "Hint: Generate inputs with input_generatorStrongScaling.py before submitting." >&2
    exit 2
fi

shopt -s nullglob
sigma_dirs=("${TRAVERSAL_DIR}"/sigmaRatio_*)
if [ "${#sigma_dirs[@]}" -ne 1 ]; then
    echo "ERROR: Expected exactly one sigmaRatio_* directory in ${TRAVERSAL_DIR}, found ${#sigma_dirs[@]}." >&2
    exit 2
fi

CELL_DIR="${sigma_dirs[0]}/cellSize_${cell_size[${SLURM_ARRAY_TASK_ID}]}"
if [ ! -d "${CELL_DIR}" ]; then
    echo "ERROR: Input cell-size directory not found: ${CELL_DIR}" >&2
    exit 2
fi

count_dirs=("${CELL_DIR}"/countRatio_*)
shopt -u nullglob
if [ "${#count_dirs[@]}" -ne 1 ]; then
    echo "ERROR: Expected exactly one countRatio_* directory in ${CELL_DIR}, found ${#count_dirs[@]}." >&2
    exit 2
fi

TARGET_BASE_DIR="${count_dirs[0]}"

# Set up environment variables
# Check your machine to see what is recommended here
export OMP_PLACES=cores
export OMP_PROC_BIND=true

unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

# Sequentially loop over thread counts and over each repeat run.
for num_threads in 1 2 4 8 14
do
    if [ "${num_threads}" -gt "${SLURM_CPUS_PER_TASK}" ]; then
        echo "ERROR: Requested ${num_threads} threads but only ${SLURM_CPUS_PER_TASK} CPUs allocated." >&2
        exit 2
    fi

    export OMP_NUM_THREADS=${num_threads}

    cd "${TARGET_BASE_DIR}"
    for run in `seq 0 2`
    do
        cd "run_${run}"
        # If not running an mpi executable, change srun to whatever is recommended for your machine.
        # Override MD_FLEXIBLE_BIN if your executable is in a different location.
        srun -n 1 -c "${num_threads}" "${MD_FLEXIBLE_BIN}" --yaml-filename input.yaml > logOutput_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_threads_${num_threads}.out 2>&1
        cd ..
    done
done

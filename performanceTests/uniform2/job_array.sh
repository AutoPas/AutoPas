#!/bin/bash
#SBATCH --job-name=Uni_2_Job_Array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny  # ! This probably needs to change !
#SBATCH --cpus-per-task=56 # ! This should be the maximum number of CPUs you wish to use per single MPI rank !
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
#SBATCH --array=0-35
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
declare -a sigma_ratio
declare -a count_ratio
declare -a data_layout
index=0
for sigma_ratio_iter in 0p05 0p15 0p25 0p35 0p45 0p55
do
    for count_ratio_iter in 0p50 1p00 2p00
    do
        for data_layout_iter in AoS SoA
        do
            sigma_ratio[$index]="$sigma_ratio_iter"
            count_ratio[$index]="$count_ratio_iter"
            data_layout[$index]="$data_layout_iter"
            index=$((index + 1))
        done
    done
done
     
# Go to directory of experiment for this job ID
TARGET_DIR="${SCRIPT_DIR}/generated_inputs/sigmaRatio_${sigma_ratio[${SLURM_ARRAY_TASK_ID}]}/countRatio_${count_ratio[${SLURM_ARRAY_TASK_ID}]}/dataLayout_${data_layout[${SLURM_ARRAY_TASK_ID}]}"
if [ ! -d "${TARGET_DIR}" ]; then
    echo "ERROR: Input directory not found: ${TARGET_DIR}" >&2
    echo "Hint: Generate inputs with input_generator.py before submitting." >&2
    exit 2
fi
cd "${TARGET_DIR}"

# Set up environment variables
# Check your machine to see what is recommended here
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=true

unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

#sequentially loop over each run
for run in `seq 0 2`
do
    cd run_${run}
    # If not running an mpi executable, change mpirun to whatever is recommended for your machine (e.g. srun)
	# Override MD_FLEXIBLE_BIN if your executable is in a different location.
    srun -n 1 "${MD_FLEXIBLE_BIN}" --yaml-filename input.yaml > logOutput_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1
    cd ..
done
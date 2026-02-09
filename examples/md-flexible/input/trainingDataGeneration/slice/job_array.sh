#!/bin/bash
#SBATCH --job-name=Slice_Job_Array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=small # ! This probably needs to change !
#SBATCH --cpus-per-task=36 # ! This should be the maximum number of CPUs you wish to use per single MPI rank !
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
#SBATCH --array=0-119

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
# with the trained model, but it is recommended to use the same binary for both data generation and actual
# use, so you may still want to use an MPI executable.
module load intel-oneapi-mpi

# Below is a common method for assigning the SLURM job array's ID to a particular directory, as created by
# the input_generator.py.

# Generate arrays of parameters indexed by job ID
declare -a slice_size
declare -a thread_count
declare -a rebuild_frequency
declare -a verlet_skin
declare -a run
index=0
for slice_size_iter in 05 10 15 20
do
    for thread_count_iter in 06 12 18 24 30 36
    do
	for rebuild_frequency_iter in 10
	do
	    for verlet_skin_iter in 0.1 0.2 0.3 0.4 0.5
	    do
     		slice_size[$index]="$slice_size_iter"
		thread_count[$index]="$thread_count_iter"
		rebuild_frequency[$index]="$rebuild_frequency_iter"
		verlet_skin[$index]="$verlet_skin_iter"
		index=$((index + 1))
	    done
	done
    done
done
    
# Go to directory of experiment for this job ID
cd slice_size_${slice_size[${SLURM_ARRAY_TASK_ID}]}/${thread_count[${SLURM_ARRAY_TASK_ID}]}_threads/skin_${verlet_skin[${SLURM_ARRAY_TASK_ID}]}/

# Set up environment variables
# Check your machine to see what is recommended here
export OMP_NUM_THREADS=${thread_count[${SLURM_ARRAY_TASK_ID}]}
export OMP_PLACES=cores
export OMP_PROC_BIND=true

unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0


#sequentially loop over each run
# If not running an mpi executable, change mpirun to whatever is recommended for your machine (e.g. srun)
for run in `seq 0 4`
do
    cd run_${run}
    
    mpirun -np 1 $HOME/AutoPas/build-MPI/examples/md-flexible/md-flexible --yaml-filename input.yaml   > logOutput_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1
    cd ..
done

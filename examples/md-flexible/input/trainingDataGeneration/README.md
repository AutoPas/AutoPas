# Training Data Generation

This directory contains scripts for generating training data to be used with machine learning tuning strategies such as random forest tuning.

It is based on that used for dataset generation by [Newcome et al. 2025](https://doi.org/10.1007/978-3-031-97635-3_35).

## WARNING

This subdirectory is merely meant to provide a starting point for generating your own data and may require modification for your computer or use case. Within the job_array(_big).sh, we have highlighted that which probably needs to be modified, but this is not an exclusive list. You must ensure that these jobscripts conform to your computer's SLURM setup.

The distribution of experiments between jobs has been roughly balanced for the computer used in this paper, and may not be appropriate for your computer and its SLURM setup.

The experiments may not be suitable as training data for the simulations you run. In particular, they do not consider the impact of different rebuild frequencies which is particularly important with dynamic containers enabled (the default option).

## How to generate data

There are five subdirectories, each with a different type of simulation:
- empty: An empty simulation (relevant generally in some large MPI simulations, as models were found to perform poorly if this data was not provided). Simulation domain: 40x40x40.
- gaussian: a number of clusters of 100 particles each's center uniformly distributed within a 40x40x40 domain, with the particles distributed about the centre in a gaussian distribution. Recreates "fuzzy" heterogeneous simulations. Different numbers of clusters are considered.
- slice: a large slice of particles in an otherwise sparsely packed 40x40x40 domain. Mimics simulations like exploding liquid. Different slice sizes are considered.
- sphere: a sphere in an otherwise empty 200x200x200 domain. Mimics simulations like heating sphere. Different sphere sizes and densities are considered.
- uniform: a number of particles uniformly distributed in a 20x20x20 domain. 100 to 2048000 particles are considered (highest density roughly corresponds to that realistic with cells of length ~6 times sigma, for Lennard-Jones molecules.)

All simulations have a cutoff of 3 and a rebuild frequency of 10.

For each simulation, 
- verlet skins between 0.1 and 0.5 are considered, in intervals of 0.1;
- thread counts between 6 and 36 are considered, in intervals of 6; and
- five runs are taken.

In each subdirectory, there is:
- `template_input.yaml`: A template input file.
- `input_generator.py`: A python script that generates a directory structure where each root directory contains an input file - one for each experiment.
- `job_array.sh`: A SLURM job script that runs a number of SLURM jobs, each running one or more experiments. There may be several to handle differently-sized simulations differently.

The job array likely needs modification:
- The script must be changed to suit your computer.
- The thread count may need to be changed - 36 was chosen as it was the number of hardware threads per NUMA core on the computer used by Newcome et al. (This may also require changing the input generator).
- Jobs may need to be merged if you are finding that the current setup produces too many too small simulations.
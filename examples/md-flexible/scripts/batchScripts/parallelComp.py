
import subprocess
import sys
import os
import time as t


def get_according_cluster_and_nodes(num_cpus: int) -> [str, str, int, int]:
    if num_cpus <= 56:
        return "serial", "serial_std", 1, 1
    if num_cpus <= 112:
        return "cm4", "cm4_tiny", 1, 1
    if num_cpus > 448:
        raise ValueError("Too many cpus")
    return "cm4", "cm4_std", 2, 4


def get_run_script(name: str, num_tasks: int, num_threads: int, executable: str, input_file: str, memory_per_task: int, time: str) -> str:
    clusters, partition, min_nodes, max_nodes = get_according_cluster_and_nodes(
        num_tasks * num_threads)
    num_nodes = min_nodes
    num_tasks_per_node = int(num_tasks / num_nodes)
    mem_per_node = memory_per_task * num_tasks_per_node
    # create output directory
    basename = name.split('.')[0] + t.strftime("%Y%m%d-%H%M%S")
    try:
        os.mkdir(basename)
    except OSError:
        print("Could not create the output directory: " + basename)
        exit(2)
    SCRIPT = f"""#!/bin/bash
#SBATCH -J {name}
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D {basename}
#SBATCH --mail-type=ALL
#SBATCH --get-user-env
#SBATCH --clusters={clusters}
#SBATCH --partition={partition}
#SBATCH --qos={partition}
#SBATCH --mem={mem_per_node}mb
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={num_tasks_per_node}
#SBATCH --cpus-per-task={num_threads}
#SBATCH --export=NONE
#SBATCH --mail-user=ge42joq@mytum.de
#SBATCH --time={time}
module load gcc/13.2.0
module load openmpi/4.1.5-gcc11
OMP_NUM_THREADS={num_threads} mpirun -np {num_tasks} {executable} --yaml-filename {input_file}
"""
    return SCRIPT


if __name__ == "__main__":

    # read args
    if len(sys.argv) < 7 or len(sys.argv) > 8:
        print(
            "Usage: python parallelComp.py <num_tasks> <num_threads> <executable> <input_file> <mem> <time> [true|false]")
        print(" - num_tasks  : number of MPI ranks to use")
        print(" - num_threads: number of threads to use per rank")
        print(" - executable : executable to use")
        print(" - input_file : input file")
        print(" - mem        : memory to use per task in MB")
        print(" - time       : time to use in HH:MM:SS")
        print(" - true/false : Step up the number of threads up until num_tasks")
        print("                (if true num_tasks must be power of 2, default: false)")
        sys.exit(1)

    stepUp = False
    if (len(sys.argv) == 8 and sys.argv[7] == "true"):
        stepUp = True

    num_tasks = int(sys.argv[1])
    num_threads = int(sys.argv[2])
    if (stepUp):
        # check if num_threads is power of 2
        if (num_tasks & (num_tasks - 1)):
            print("num_tasks must be power of 2")
            sys.exit(1)

    executable = os.path.abspath(sys.argv[3])
    input_file = os.path.abspath(sys.argv[4])
    mem = int(sys.argv[5])
    time = sys.argv[6]
    # check time format
    if (len(time.split(":")) != 3):
        print("time format must be HH:MM:SS")
        sys.exit(1)

    current_tasks = 1
    if (not stepUp):
        current_tasks = num_tasks

    base_name = input_file.split("/")[-1].split(".")[0]

    while current_tasks <= num_tasks:
        name = f"{base_name}_{current_tasks}"
        run_script = get_run_script(
            name, current_tasks, num_threads, executable, input_file, mem, time)

        with open(f"{name}.sh", "w") as f:
            f.write(run_script)

        print(f"Running {name}.sh")
        subprocess.call(["sbatch", f"{name}.sh"])

        current_tasks *= 2


import subprocess
import sys
import os


def get_according_cluster_and_nodes(num_cpus: int) -> [str, str, int, int]:
    if num_cpus <= 56:
        return "serial", "serial_std", 1, 1
    if num_cpus <= 112:
        return "cm4", "cm4_tiny", 1, 1
    if num_cpus > 448:
        raise ValueError("Too many cpus")
    return "cm4", "cm4_std", 2, 4


def get_run_script(name: str, num_tasks: int, num_threads: int, executable: str, input_file: str, memory: int, time: str) -> str:
    clusters, partition, min_nodes, max_nodes = get_according_cluster_and_nodes(
        num_tasks * num_threads)
    num_nodes = min_nodes
    num_tasks_per_node = int(num_tasks / num_nodes)
    SCRIPT = f"""#!/bin/bash
#SBATCH -J {name}
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D .
#SBATCH --mail-type=ALL
#SBATCH --get-user-env
#SBATCH --clusters={clusters}
#SBATCH --partition={partition}
#SBATCH --qos={partition}
#SBATCH --mem={memory}mb
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={num_tasks_per_node}
#SBATCH --cpus-per-task={num_threads}
#SBATCH --export=NONE
#SBATCH --mail-user=ge42joq@mytum.de
#SBATCH --time={time}
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
        print(" - mem        : memory to use in MB")
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

        subprocess.call(["sbatch", f"{name}.sh"])

        current_tasks *= 2

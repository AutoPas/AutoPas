
import subprocess
import sys
import os
import time as t


def get_according_cluster_and_nodes(num_cpus: int) -> [str, str, int, int]:
    if num_cpus <= 16:
        return "serial", "serial_std", 1, 1
    if num_cpus <= 224:
        return "cm4", "cm4_tiny", 1, 1
    if num_cpus > 896:
        raise ValueError("Too many cpus")
    return "cm4", "cm4_std", 2, 4


# pessimistic parallel portion
optimized_portion = 0.3


def get_run_script(name: str, num_tasks: int, num_threads: int,
                   executable: str, input_file: str, memory_per_task: int,
                   time: int, adapt_time: bool, cm4: bool) -> str:
    num_cpus = num_tasks * num_threads
    num_assigned_threads = num_threads
    if cm4 and num_cpus <= 16:
        num_assigned_threads = (int(17 / num_tasks) + 1)
        num_cpus = num_tasks * num_assigned_threads
    clusters, partition, min_nodes, max_nodes = get_according_cluster_and_nodes(
        num_cpus)
    num_nodes = min_nodes
    num_tasks_per_node = int(num_tasks / num_nodes)
    mem_per_node = memory_per_task * num_tasks_per_node
    # create output directory
    basename = name.split('.')[0] + t.strftime("-%Y%m%d-%H%M%S")
    try:
        os.mkdir(basename)
    except OSError:
        print("Could not create the output directory: " + basename)
        exit(2)
    qos = f"\n#SBATCH --qos={partition}\n"
    if clusters == "serial":
        qos = ""
    # adapt time if needed with Amdahl's law
    if adapt_time:
        timeOptimized = time * optimized_portion
        timeSerial = time - timeOptimized
        time = timeSerial + timeOptimized / num_tasks
    # put time into HH:MM:SS
    hours = int(time / 60)
    minutes = int(time % 60)
    time_str = f"{hours:02d}:{minutes:02d}:00"
    SCRIPT = f"""#!/bin/bash
#SBATCH -J {name}
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D {basename}
#SBATCH --mail-type=ALL
#SBATCH --get-user-env
#SBATCH --clusters={clusters}
#SBATCH --partition={partition}{qos}
#SBATCH --mem={mem_per_node}mb
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={num_tasks_per_node}
#SBATCH --cpus-per-task={num_assigned_threads}
#SBATCH --export=NONE
#SBATCH --mail-user=ge42joq@mytum.de
#SBATCH --time={time_str}
module load gcc/13.2.0
module load openmpi/4.1.5-gcc11
OMP_NUM_THREADS={num_threads} mpirun -np {num_tasks} {executable} --yaml-filename {input_file}
"""
    return SCRIPT


if __name__ == "__main__":

    # read args
    if len(sys.argv) < 7 or len(sys.argv) > 9:
        print(
            "Usage: python parallelComp.py <num_tasks> <num_threads> <executable> <input_file> <mem> <time> [true|false] [cm4]")
        print(" - num_tasks  : number of MPI ranks to use")
        print(" - num_threads: number of threads to use per rank")
        print(" - executable : executable to use")
        print(" - input_file : input file")
        print(" - mem        : memory to use per task in MB")
        print(" - time       : expected time in minutes (longest time if with stepup)")
        print(" - true/false : Step up the number of threads up until num_tasks")
        print("                (if true num_tasks must be power of 2, default: false)")
        print("                This will also adapt the expected time heuristically.")
        sys.exit(1)

    cm4 = False
    stepUp = False
    if (len(sys.argv) >= 8 and sys.argv[7] == "true"):
        stepUp = True
    if (len(sys.argv) == 9 and sys.argv[8] == "cm4"):
        cm4 = True

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
    time = int(sys.argv[6])

    current_tasks = 1
    if (not stepUp):
        current_tasks = num_tasks

    base_name = input_file.split("/")[-1].split(".")[0]

    adaptTime = False

    while current_tasks <= num_tasks:
        name = f"{base_name}_{current_tasks}"
        run_script = get_run_script(
            name, current_tasks, num_threads, executable, input_file, mem, time, adaptTime, cm4)

        with open(f"{name}.sh", "w") as f:
            f.write(run_script)

        print(f"Running {name}.sh")
        subprocess.call(["sbatch", f"{name}.sh"])

        current_tasks *= 2
        adaptTime = True

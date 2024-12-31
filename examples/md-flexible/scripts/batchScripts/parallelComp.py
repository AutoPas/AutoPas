
import subprocess
import sys
import os
import time

def get_run_script(name: str, num_threads: str, input_file: str, memory:int, time: str,
                   clusters: str = "cm4", partition: str = "cm4_tiny") -> str:
    SCRIPT = f"""#!/bin/bash
#SBATCH -J {name}
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D .
#SBATCH --get-user-env
#SBATCH --clusters={clusters}
#SBATCH --partition={partition}
#SBATCH --mem={memory}mb
#SBATCH --cpus-per-task={num_threads}
#SBATCH --export=NONE
#SBATCH --mail-user=ge42joq@mytum.de
#SBATCH --time={time}
OMP_NUM_THREADS=mpirun -np {num_threads} ./examples/md-flexible/md-flexible --yaml-filename {input_file}
"""
    return SCRIPT


if __name__ == "__main__":

    # read args
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print("Usage: python parallelComp.py <num_threads> <input_file> <mem> <time> [true|false]")
        print(" - num_threads: number of threads to use")
        print(" - input_file : input file")
        print(" - mem        : memory to use in MB")
        print(" - time       : time to use in HH:MM:SS")
        print(" - true/false : Step up the number of threads up until num_threads")
        print("                (if true num_threads must be power of 2, default: false)")
        sys.exit(1)

    stepUp = False
    if (len(sys.argv) == 6 and sys.argv[5] == "true"):
        stepUp = True


    num_threads = int(sys.argv[1])
    if(stepUp):
        # check if num_threads is power of 2
        if (num_threads & (num_threads - 1)):
            print("num_threads must be power of 2")
            sys.exit(1)

    input_file = os.path.abspath(sys.argv[2])
    mem = int(sys.argv[3])
    time = sys.argv[4]
    #check time format
    if (len(time.split(":")) != 3):
        print("time format must be HH:MM:SS")
        sys.exit(1)


    current_threads = 1
    if (not stepUp):
        current_threads = num_threads

    base_name = input_file.split("/")[-1]
    base_name = base_name.split(".")[0]

    while current_threads <= num_threads:
        name = f"{base_name}_{current_threads}"
        run_script = ""
        if (num_threads < 16):
            run_script = get_run_script(
                name, current_threads, input_file, mem, time)
        else:
            run_script = get_run_script(
                name, current_threads, input_file, mem, time, "cm4", "cm4_std")

        with open(f"{name}.sh", "w") as f:
            f.write(run_script)

        subprocess.call(["sbatch", f"{name}.sh"])

        current_threads *= 2

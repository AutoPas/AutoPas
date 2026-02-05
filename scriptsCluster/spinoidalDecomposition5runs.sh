#!/bin/bash
#SBATCH -J NoEarlyStoppingVCLC06
#SBATCH -o %x.%j.%N.out
#SBATCH -e %x.%j.%N.err
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=200G
#SBATCH --time=24:00:00
set -euo pipefail

module load slurm_setup
module load gcc/15.2.0

PROJECT_DIR="$SLURM_SUBMIT_DIR"
BIN="$PROJECT_DIR/../build/examples/md-flexible/md-flexible"
YAML="$PROJECT_DIR/../build/examples/md-flexible/SpinodalDecomposition_equilibration.yaml"

RUNBASE="$PROJECT_DIR/runs/virialSpin_$SLURM_JOB_ID"
mkdir -p "$RUNBASE"

echo "Running in $RUNBASE"

for i in 1 2 3; do
    RUNDIR="$RUNBASE/run_$i"
    mkdir -p "$RUNDIR"
    cd "$RUNDIR"

    echo "Normal run $i"
    "$BIN" --yaml-filename "$YAML"
done


# ---------- gprof ----------
#RUNDIR="$RUNBASE/gprof"
#mkdir -p "$RUNDIR"
#cd "$RUNDIR"

#"$BIN" --yaml-filename "$YAML"

#gprof "$BIN" gmon.out > profile_report_gprof.txt
#echo "gprof profiling completed."
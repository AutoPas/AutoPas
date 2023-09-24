#!/bin/bash
export LC_NUMERIC=en_US.UTF-8
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    echo "Usage: $0 PATH_TO_MD-Flexible"
    exit -1
fi

MDFlex=$1

configPrinted=false

for t in 1 2 4 8
do
    export OMP_NUM_THREADS=${t}
    output=$(${MDFlex} \
    --box-length                        40 \
    --container                         lc \
    --cutoff                            1 \
    --distribution-mean                 10 \
    --distribution-stddeviation         6 \
    --data-layout                       aos \
    --functor                           lj \
    --iterations                        10 \
    --particles-generator               grid \
    --particles-per-dimension           80 \
    --particle-spacing                  .5 \
    --traversal                         sliced \
    --verlet-rebuild-frequency          20 \
    --verlet-skin-radius-per-timestep   .01 \
    )

    if [ "${configPrinted}" = false ] ; then
        configPrinted=true
        # print all output lines until, excluding, "Using" (this is the whole config part)
        sed '/Using/Q' <<< "${output}"
        echo
        printf "%6s%15s%15s%15s%24s\n" "Threads" "GFLOPs/s" "MFUPs/s" "Time[micros]" "SingleIteration[micros]"
    fi

    gflops=$(echo "$output" | sed --quiet -e 's|GFLOPs/sec.*: \(.*\)|\1|gp')
    mfups=$(echo "$output" | sed --quiet -e 's|MFUPs/sec.*: \(.*\)|\1|gp')
    timeTotal=$(echo "$output" | sed --quiet -e 's|Time total.*: \(.*\) .*s (.*|\1|gp')
    timeOne=$(echo "$output" | sed --quiet -e 's|One iteration.*: \(.*\) .*s (.*|\1|gp')

    printf "%6d%15.2f%15.2f%15.2f%24.2f\n" "$t" "$gflops" "$mfups" "$timeTotal" "$timeOne"
done
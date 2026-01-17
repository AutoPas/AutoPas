#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

BINARY="./build-remote-release/examples/md-flexible/md-flexible"
PARTICLES=64000
TYPE="GPU"
CSV_FILE="results.csv"
BASE_LENGTH=40

echo "num_particles,type,scale,side_length,mfups_sec" > "$CSV_FILE"

# Iterating from 0.1 up to 3.0
for i in $(seq 1 1 30); do
    VALS=$(awk "BEGIN { scale=$i/10; len=$BASE_LENGTH*scale; printf \"%.4f %.4f\", scale, len }")
    SCALE=$(echo $VALS | cut -d' ' -f1)
    LENGTH_FMT=$(echo $VALS | cut -d' ' -f2)

    echo "--- Scale: $SCALE | Length: $LENGTH_FMT ---"

    # tee /dev/tty allows live viewing while RAW_OUT captures the stream for awk
    RAW_OUT=$($BINARY \
      --container KokkosVerletClusterLists \
      --traversal kk_vcl \
      --functor lennard-jones-kokkos \
      --data-layout SoA \
      --boundary-type reflective reflective reflective \
      --box-length "$LENGTH_FMT" "$LENGTH_FMT" "$LENGTH_FMT" \
      --newton3 enabled \
      --particle-generator uniform \
      --particles-total $PARTICLES \
      --deltaT 0.0 \
      --iterations 1000 \
      --cutoff 3 \
      --no-end-config \
      --log-level info | tee /dev/tty)

    MFUPS=$(echo "$RAW_OUT" | awk '/MFUPs\/sec/ {print $3}')

    if [ -n "$MFUPS" ]; then
      echo "$PARTICLES,$TYPE,$SCALE,$LENGTH_FMT,$MFUPS" >> "$CSV_FILE"
    fi
done
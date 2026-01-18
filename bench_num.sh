#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

BINARY="./build-remote-release/examples/md-flexible/md-flexible"
TYPE="GPU"
CSV_FILE="results.csv"

echo "num_particles,type,side_length,mfups_sec" > "$CSV_FILE"

# Iterating side lengths to vary particle count at density = 1
# N = V * density = L^3 * 1
for L in $(seq 10 5 100); do
    PARTICLES=$(awk "BEGIN { print $L^3 }")

    echo "--- Side Length: $L | Particles: $PARTICLES ---"
#      --verlet-rebuild-frequency 1 \

    RAW_OUT=$($BINARY \
      --container KokkosVerletClusterLists \
      --traversal kk_vcl \
      --functor lennard-jones-kokkos \
      --data-layout SoA \
      --boundary-type reflective reflective reflective \
      --box-length "$L" "$L" "$L" \
      --newton3 enabled \
      --particle-generator uniform \
      --particles-total "$PARTICLES" \
      --deltaT 0.0 \
      --iterations 100 \
      --cutoff 3 \
      --no-end-config \
      --log-level info | tee /dev/tty)

    MFUPS=$(echo "$RAW_OUT" | awk '/MFUPs\/sec/ {print $3}')

    if [ -n "$MFUPS" ]; then
      echo "$PARTICLES,$TYPE,$L,$MFUPS" >> "$CSV_FILE"
    fi
done
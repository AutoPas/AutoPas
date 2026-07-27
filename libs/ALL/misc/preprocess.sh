#!/bin/bash

# Usage:
# $1: path to include directory
# $2: output file
# $3..: additional defines (without -D)

# Either run this script in the include directory, or pass the include
# directory as an additional argument ($1).

set -euo pipefail

set -x

HEADERS="ALL_CustomExceptions.hpp
ALL_Defines.h
ALL_ForceBased.hpp
ALL_Functions.hpp
ALL_Histogram.hpp
ALL.hpp
ALL_LB.hpp
ALL_Point.hpp
ALL_Staggered.hpp
ALL_Tensor.hpp
ALL_Voronoi.hpp"

PROCESSED_DIR=`mktemp -d`

CWD=`pwd`
cd "$1"
if [[ ${2:0:1} == / ]]
then
	OUTFILE="$2"
else
	OUTFILE="$CWD/$2"
fi
shift
shift

for f in $HEADERS
do
	cat $f | sed -e 's!#\s*include\s*<!// PPIGNORE <!' > "$PROCESSED_DIR/$f"
done

cd "$PROCESSED_DIR"

DEFINES=
for d in "$@"
do
	DEFINES="$DEFINES -D$d"
done

cpp -undef -fdirectives-only $DEFINES -o ALL.ii ALL.hpp

mkdir -p "${OUTFILE%/*}"
cat ALL.ii | sed -e 's!// PPIGNORE <!#include <!' > "$OUTFILE"

rm -rf $PROCESSED_DIR

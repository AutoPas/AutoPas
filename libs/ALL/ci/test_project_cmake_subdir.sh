#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")"; pwd -P)"
source $SCRIPT_DIR/ci_funcs.sh
set -euxo pipefail

find_ALL #sets ALL_ROOT_DIR
find_VTK #sets VTK_DIR

TEMP=$(mktemp -d)

cp -rL "$CI_SCRIPT_PATH/../example/CMakeProjectSubdir/." "$TEMP"
cd "$TEMP"
ln -s "$ALL_ROOT_DIR" all
ln -s "$VTK_DIR" vtk_bin
./build_all.sh


rm -rf "$TEMP"

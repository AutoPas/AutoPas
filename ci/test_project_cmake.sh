#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")"; pwd -P)"
source $SCRIPT_DIR/ci_funcs.sh
set -euxo pipefail

find_ALL #sets ALL_ROOT_DIR
find_VTK #sets VTK_DIR

TEMP=$(mktemp -d)

cp -rL "$CI_SCRIPT_PATH/../example/CMakeProject/." "$TEMP"
cd "$TEMP"
sed -i\
	-e "/^ALL_ROOT_DIR=/c\\\nALL_ROOT_DIR=\"$ALL_ROOT_DIR\"\n"\
	-e "/^VTK_DIR=/c\\\nVTK_DIR=\"$VTK_DIR\"\n"\
	build_all.sh
cat build_all.sh
./build_all.sh


rm -rf "$TEMP"

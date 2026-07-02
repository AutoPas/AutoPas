#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")"; pwd -P)"
source $SCRIPT_DIR/ci_funcs.sh
set -euxo pipefail

find_ALL #sets ALL_ROOT_DIR
find_VTK #sets VTK_DIR

TEMP=$(mktemp -d)

cp -rL "$CI_SCRIPT_PATH/../example/MakefileProject/." "$TEMP"
cd "$TEMP"
sed -i\
	-e "/^ALL_SOURCE=/c\\\nALL_SOURCE=\"$ALL_ROOT_DIR\"\n"\
	-e "/^VTK_DIR=/c\\\nVTK_DIR=\"$VTK_DIR\"\n"\
	Makefile
cat Makefile
make


rm -rf "$TEMP"

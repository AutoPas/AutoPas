#!/bin/bash

export VTK_DIR="`pwd`/vtk_bin/lib/cmake/vtk-7.1"

prepare_example () {
	cat $1 | sed \
		-e 's/ALL_VTK_OUTPUT/TEST_VTK_OUTPUT/'\
		-e 's/ALL_VORONOI_ACTIVE/TEST_VORONOI/'\
		> $2
}

prepare_example ALL_test_src.cpp ALL_test.cpp
rm -rf build && CC=gcc CXX=g++ cmake -S . -B build && VERBOSE=1 cmake --build build

#!/bin/bash

# all commands must execute successfully
set -e
set -o pipefail
set -u
set -x

# $ALL_INSTALL_DIR must be an absolute path!
ALL_ROOT_DIR=../..
ALL_BUILD_DIR=all_build
ALL_INSTALL_DIR=`pwd`/all_bin #where ALL installs itself to
ALL_PACKAGE=`pwd`/all_package #where our project expects the installed ALL
VTK_DIR=`pwd`/../../../vtk_bin
BUILD_DIR=build

# We only move ALL after installation from $ALL_INSTALL_DIR to $ALL_PACKAGE to
# test for errors in relocatability of the library. A typical user does not
# need to do that and can just set $ALL_PACKAGE to the directory ALL installs
# itself to.

export CC=gcc
export CXX=g++

build_all () {
	rm -rf "$ALL_BUILD_DIR"
	mkdir "$ALL_BUILD_DIR"
	cmake -S "$ALL_ROOT_DIR" -B "$ALL_BUILD_DIR"\
		-DCMAKE_INSTALL_PREFIX="$ALL_INSTALL_DIR"\
		-DCM_ALL_FORTRAN=ON\
		-DCM_ALL_VTK_OUTPUT=ON\
		-DCM_ALL_VORONOI=ON\
		-DCMAKE_BUILD_TYPE=Release\
		-DVTK_DIR="$VTK_DIR"/lib/cmake/vtk-7.1

	cmake --build "$ALL_BUILD_DIR"
	rm -rf "ALL_INSTALL_DIR"
	cmake --install "$ALL_BUILD_DIR"
	if [[ $ALL_INSTALL_DIR != $ALL_PACKAGE ]]
	then
		rm -rf "$ALL_PACKAGE"
		mv "$ALL_INSTALL_DIR" "$ALL_PACKAGE"
	fi
}

build_self () {
	rm -rf "$BUILD_DIR"
	mkdir "$BUILD_DIR"
	cmake -S . -B "$BUILD_DIR"\
		-DALL_DIR="$ALL_PACKAGE"/lib/cmake/ALL\
		-DVTK_DIR="$VTK_DIR"/lib/cmake/vtk-7.1

	VERBOSE=1 cmake --build "$BUILD_DIR"
}

prepare_example () {
	cat $1 | sed \
		-e 's/ALL_VTK_OUTPUT/TEST_VTK_OUTPUT/'\
		-e 's/ALL_VORONOI_ACTIVE/TEST_VORONOI/'\
		> $2
}

prepare_example ALL_test_src.cpp ALL_test.cpp

build_all
build_self

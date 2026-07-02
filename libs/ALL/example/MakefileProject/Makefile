# Set our compilers
CC=gcc
CXX=g++
FC=gfortran
# as well as MPI wrapper
MCC=mpicc
MCXX=mpic++
MFC=mpifort

# This Makefile will build ALL automatically, but not VTK. That is assumed to
# be available in $(VTK_DIR). It can be built automatically by reproducing the
# way ALL is included. In general, changes to the library do not trigger a
# rebuild, so in those cases, you need to run `make clean-ext-liball`
# explicitly to trigger the rebuild.
#
# This Makefile also enables the Voronoi and VTK features of ALL. If these are
# not necessary, you can remove the corresponding lines.

ALL_DIR=all_bin
VTK_DIR=vtk_bin

# Source tree and temporary build directory
ALL_BUILD=all_build
ALL_SOURCE=all
ALL_SOURCE=../..

# Which features should be enabled for ALL
LIBALL_CONFIGURE := -DCM_ALL_FORTRAN=ON -DCM_ALL_VORONOI=ON -DCM_ALL_VTK_OUTPUT=ON -DCMAKE_BUILD_TYPE=Release


LIBALL_CONFIG_FILE := $(ALL_DIR)/lib/make/Makefile.libALL.in
# This will provide the variables LIBALL_DEFINES, LIBALL_INCLUDE and
# LIBALL_LIB, which need to be passed to the compiler and linker.
-include $(LIBALL_CONFIG_FILE)

# Since gfortran does not link the standard C++ library by default, we need to
# explicitly state it.
FFLAGS += -lstdc++

FLAGS += -DTEST_VTK_OUTPUT -DTEST_VORONOI
FLAGS += $(LIBALL_DEFINES)

# The linking order of VTK can be gleaned from the CMake invocation of any
# executable of the library. Just compile the library using `VERBOSE=1` and the
# command lines are output as well.
VTK_INCLUDE := -I$(VTK_DIR)/include/vtk-7.1
VTK_LIB := $(subst lib/,$(VTK_DIR)/lib/, lib/libvtkFiltersProgrammable-7.1.a lib/libvtkIOParallelXML-7.1.a lib/libvtkIOXML-7.1.a lib/libvtkIOXMLParser-7.1.a lib/libvtkexpat-7.1.a lib/libvtkParallelMPI-7.1.a lib/libvtkParallelCore-7.1.a lib/libvtkIOLegacy-7.1.a lib/libvtkIOCore-7.1.a lib/libvtkCommonExecutionModel-7.1.a lib/libvtkCommonDataModel-7.1.a lib/libvtkCommonTransforms-7.1.a lib/libvtkCommonMisc-7.1.a lib/libvtkCommonMath-7.1.a lib/libvtkCommonSystem-7.1.a lib/libvtkCommonCore-7.1.a lib/libvtksys-7.1.a -ldl -lpthread lib/libvtkzlib-7.1.a)

INCLUDES += $(LIBALL_INCLUDE)
INCLUDES += $(VTK_INCLUDE)
LIBS += $(LIBALL_LIB) $(VTK_LIB)

.PHONY: all clean clean-ext clean-ext-liball

all: ALL_Staggered ALL_Staggered_f ALL_test

clean:
	rm -rf ALL_Staggered ALL_Staggered_f ALL_test.cpp ALL_test

clean-ext: clean-ext-liball

clean-ext-liball:
	rm -rf $(ALL_BUILD) $(ALL_DIR)

# Our projects executables:
ALL_Staggered: ALL_Staggered.cpp $(LIBALL_CONFIG_FILE)
	$(MCXX) $(FLAGS) $(CXXFLAGS) $(INCLUDES) $< $(LIBS) -o $@

ALL_Staggered_f: ALL_Staggered_f.F90 $(LIBALL_CONFIG_FILE)
	$(MFC) $(FLAGS) $(FFLAGS) $(INCLUDES) $< $(LIBS) -o $@

ALL_test.cpp: ALL_test_src.cpp
	$(shell cat $< | sed -e 's/ALL_VTK_OUTPUT/TEST_VTK_OUTPUT/' -e 's/ALL_VORONOI_ACTIVE/TEST_VORONOI/' > $@)

ALL_test: ALL_test.cpp $(LIBALL_CONFIG_FILE)
	$(MCXX) $(FLAGS) $(CXXFLAGS) $(INCLUDES) $< $(LIBS) -o $@

# The library part follows:
$(LIBALL_CONFIG_FILE):
	rm -rf $(ALL_BUILD)
	mkdir -p $(ALL_BUILD)
	CC=$(CC) CXX=$(CXX) FC=$(FC) VTK_DIR=$(VTK_DIR)/lib/cmake/vtk-7.1 cmake -S $(ALL_SOURCE) -B $(ALL_BUILD)  -DCMAKE_INSTALL_PREFIX=`pwd`/$(ALL_DIR) $(LIBALL_CONFIGURE)
	cmake --build $(ALL_BUILD)
	rm -rf $(ALL_DIR)
	cmake --install $(ALL_BUILD)

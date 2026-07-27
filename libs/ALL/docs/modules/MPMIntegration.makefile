MJOBS ?= $(shell getconf _NPROCESSORS_ONLN)
JUELICH_ALL_INCLUDE := external/juelich_all_build/include/modules
JUELICH_ALL_LIB := external/juelich_all_build/lib/libALL_fortran.a external/juelich_all_build/lib/libALL.a
VTK_LIB := $(subst lib/,external/vtk_build/lib/, lib/libvtkFiltersProgrammable-7.1.a lib/libvtkIOParallelXML-7.1.a lib/libvtkIOXML-7.1.a lib/libvtkIOXMLParser-7.1.a lib/libvtkexpat-7.1.a lib/libvtkParallelMPI-7.1.a lib/libvtkParallelCore-7.1.a lib/libvtkIOLegacy-7.1.a lib/libvtkIOCore-7.1.a lib/libvtkCommonExecutionModel-7.1.a lib/libvtkCommonDataModel-7.1.a lib/libvtkCommonTransforms-7.1.a lib/libvtkCommonMisc-7.1.a lib/libvtkCommonMath-7.1.a lib/libvtkCommonSystem-7.1.a lib/libvtkCommonCore-7.1.a lib/libvtksys-7.1.a -ldl -lpthread lib/libvtkzlib-7.1.a)

# ...

# VTK
VTKCONFIG_FILE := external/vtk_build/lib/cmake/vtk-7.1/VTKConfig.cmake
$(VTKCONFIG_FILE):
	mkdir -p external/vtk_build
	cd external/vtk_build && CC=$(CC) CXX=$(CXX) cmake ../vtk -DCMAKE_INSTALL_PREFIX=`pwd` $(EXT_LTO_CMFLAGS) -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DVTK_Group_MPI=OFF -DVTK_Group_Rendering=OFF -DVTK_Group_StandAlone=OFF -DVTK_RENDERING_BACKEND=None -DVTK_USE_CXX11_FEATURES=ON -DModule_vtkCommonDataModel=ON -DModule_vtkFiltersProgrammable=ON -DModule_vtkIOParallelXML=ON -DModule_vtkParallelMPI=ON
	$(MAKE) -j $(MJOBS) -C external/vtk_build
	$(MAKE) -C external/vtk_build install

# juelich_all
$(JUELICH_ALL_LIB): $(VTKCONFIG_FILE)
	mkdir -p external/juelich_all_build
	cd external/juelich_all_build && CC=$(CC) CXX=$(CXX) cmake ../juelich_all $(EXT_LTO_CMFLAGS) -DCMAKE_Fortran_FLAGS="-Wall -Wextra -fbacktrace $(EXT_LTO)" -DCMAKE_INSTALL_PREFIX=`pwd` -DCM_ALL_FORTRAN=ON -DCM_ALL_USE_F08=$(ALL_USE_F08) -DCMAKE_BUILD_TYPE=Release -DCM_ALL_DEBUG=OFF -DCM_ALL_VTK_OUTPUT=ON -DVTK_DIR=`pwd`/../vtk_build/lib/cmake/vtk-7.1
	$(MAKE) -C external/juelich_all_build
	$(MAKE) -C external/juelich_all_build install


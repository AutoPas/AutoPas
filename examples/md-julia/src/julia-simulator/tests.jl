include("testBoundaryCondition.jl")

executeTest()

CMake Error in /mnt/c/Users/laura/Documents/BA_SH/testB/AutoPas/build/CMakeFiles/CMakeTmp/CMakeLists.txt:
  Imported target "MPI::MPI_C" includes non-existent path

    "/usr/lib/x86_64-linux-gnu/openmpi/include"

  in its INTERFACE_INCLUDE_DIRECTORIES.  Possible reasons include:

  * The path was deleted, renamed, or moved to another location.

  * An install or uninstall procedure did not complete successfully.

  * The installation package was faulty and references files it does not
  provide.



CMake Error in /mnt/c/Users/laura/Documents/BA_SH/testB/AutoPas/build/CMakeFiles/CMakeTmp/CMakeLists.txt:
  Imported target "MPI::MPI_C" includes non-existent path

    "/usr/lib/x86_64-linux-gnu/openmpi/include"

  in its INTERFACE_INCLUDE_DIRECTORIES.  Possible reasons include:

  * The path was deleted, renamed, or moved to another location.

  * An install or uninstall procedure did not complete successfully.

  * The installation package was faulty and references files it does not
  provide.



CMake Error at /usr/share/cmake-3.22/Modules/FindMPI.cmake:1264 (try_compile):
  Failed to generate test project build system.
Call Stack (most recent call first):
  /usr/share/cmake-3.22/Modules/FindMPI.cmake:1315 (_MPI_try_staged_settings)
  /usr/share/cmake-3.22/Modules/FindMPI.cmake:1638 (_MPI_check_lang_works)
  cmake/modules/autopas_mpi.cmake:16 (find_package)
  CMakeLists.txt:83 (include)

  /usr/lib/x86_64-linux-gnu/openmpi/include"

  ./configure --prefix=/usr/lib/x86_64-linux-gnu/openmpi 2>&1 | tee config.out
  make all 2>&1 | tee make.out
  make install 2>&1 | tee install.out
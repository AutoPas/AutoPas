# REMEMBER TO SET CC AND CXX!
# CC=gcc CXX=g++ cmake -S . -B build -DALL_DIR=`pwd`/../../jall_bin/lib/cmake/ALL

cmake_minimum_required(VERSION 3.14)
project(ExampleProject)

enable_language(CXX)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# For cleaner code we should also include MPI here, since this project's source
# also needs MPI. Otherwise, we are free to leave it out, since the parts
# required by ALL are automatically included.
find_package(MPI REQUIRED COMPONENTS CXX Fortran)

# Some variants for testing VTK output and the voronoi method. These are only
# required, so the code used in .this. project has the corresponding features
# enabled and calls the respective functions from the library. This is not
# required for the library itself. Although the library needs have been built
# with the respective features enabled.
add_compile_definitions(TEST_VTK_OUTPUT)
add_compile_definitions(TEST_VORONOI)

find_package(ALL 0.9.1)

add_executable(ALL_Staggered ALL_Staggered.cpp)
target_link_libraries(ALL_Staggered PUBLIC ALL::ALL)
# This is in this case not .necessary., since it is already inherited from ALL,
# however, for readability and clear intent, we should still link against MPI
# ourselves, since our code also uses MPI directly.
target_link_libraries(ALL_Staggered PUBLIC ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
target_include_directories(ALL_Staggered PRIVATE ${MPI_CXX_INCLUDE_PATH})

add_executable(ALL_Staggered_f ALL_Staggered_f.F90)
set_property(TARGET ALL_Staggered_f PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(ALL_Staggered_f PUBLIC ALL::ALL_fortran)
target_link_libraries(ALL_Staggered_f PUBLIC ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
target_include_directories(ALL_Staggered_f PRIVATE ${MPI_Fortran_INCLUDE_PATH})

add_executable(ALL_test ALL_test.cpp)
target_link_libraries(ALL_test PUBLIC ALL::ALL)
target_link_libraries(ALL_test PUBLIC ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
target_include_directories(ALL_test PRIVATE ${MPI_CXX_INCLUDE_PATH})

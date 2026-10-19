# LB4OMP (https://github.com/unibas-dmi-hpc/LB4OMP) replaces LLVM's OpenMP runtime libomp, adding further scheduling
# techniques and Auto4OMP's selection methods.
# Kept as a git subtree in libs/lb4omp, with some fixes added to the upstream version.

set(AUTOPAS_LB4OMP_DOC "Uses LB4OMP as the OpenMP runtime, adding its scheduling techniques. Requires Clang.")
option(AUTOPAS_LB4OMP "${AUTOPAS_LB4OMP_DOC}" OFF)

if (NOT AUTOPAS_LB4OMP)
    return()
endif ()
if (NOT AUTOPAS_OPENMP)
    message(WARNING "OpenMP must be enabled to use LB4OMP. LB4OMP disabled.")
    set(AUTOPAS_LB4OMP OFF CACHE BOOL "${AUTOPAS_LB4OMP_DOC}" FORCE)
    return()
endif ()
if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "LB4OMP requires Clang. Use Clang, or disable LB4OMP with -DAUTOPAS_LB4OMP=OFF.")
endif ()
message(STATUS "LB4OMP - using bundled version (commit b9d25eca) with AutoPas fixes")

# Build LB4OMP for OpenMP 4.5, the version AutoPas is compiled for with LB4OMP (see -fopenmp-version below).
set(LIBOMP_OMP_VERSION 45 CACHE INTERNAL "")
# Offloading is not needed, and its build fails if CUDA is installed.
set(OPENMP_ENABLE_LIBOMPTARGET OFF CACHE INTERNAL "")
# Don't copy omp.h and libomp.so into libs/lb4omp/runtime/exports; AutoPas uses the copies in the build tree.
set(LIBOMP_COPY_EXPORTS OFF CACHE INTERNAL "")
# LB4OMP's scheduling code needs a high resolution timer, but LB4OMP only looks for one if LIBOMP_STATS is on.
if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86|AMD64")
    set(LIBOMP_HAVE___RDTSC ON CACHE INTERNAL "")
    set(LIBOMP_HAVE_X86INTRIN_H ON CACHE INTERNAL "")
else ()
    set(LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER ON CACHE INTERNAL "")
endif ()
# Build LB4OMP on its own, not as part of an LLVM build.
set(OPENMP_STANDALONE_BUILD ON)

include(FetchContent)
FetchContent_Declare(lb4omp SOURCE_DIR ${AUTOPAS_SOURCE_DIR}/libs/lb4omp)
# LB4OMP's test setup uses the FindPythonInterp module. Allow it without a warning.
set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
FetchContent_MakeAvailable(lb4omp)
unset(CMAKE_POLICY_DEFAULT_CMP0148)

# Hide LB4OMP's cache variables from ccmake / cmake-gui.
get_cmake_property(cacheVariables CACHE_VARIABLES)
list(FILTER cacheVariables INCLUDE REGEX "^(LIBOMP|OPENMP)")
mark_as_advanced(${cacheVariables})

# OpenMP::OpenMP_CXX is not linked with LB4OMP, so: -fopenmp enables pragma recognition and __kmpc_* code generation
# (required for LB4OMP dispatch). It must be a compile option (not a link option) so the compiler emits
# __kmpc_dispatch_* calls without also implicitly linking the system libomp; only LB4OMP's libomp.so is linked.
# -fopenmp-version=45: LB4OMP is not compatible with OpenMP 5.0 or newer. 
# The include directory holds LB4OMP's generated omp.h.
target_compile_options(omp INTERFACE -fopenmp -fopenmp-version=45)
target_include_directories(omp INTERFACE ${lb4omp_BINARY_DIR}/runtime/src)
add_library(lb4omp::omp ALIAS omp)
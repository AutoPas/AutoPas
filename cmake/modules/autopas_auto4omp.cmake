# Version settings for convenience. If Auto4OMP supports new packages, update here.
set(AUTOPAS_NVCC_MAX_GCC 12)

# Auto4OMP option description.
set(AUTOPAS_AUTO4OMP_DOC "Activates Auto4OMP to automatically select OpenMP scheduling algorithms \
from the LB4OMP portfolio during runtime.")

# NVCC GNUC path option description.
set(AUTOPAS_NVCC_GNUC_PATH_DOC "Path to gcc ${AUTOPAS_NVCC_MAX_GCC} or less, required by NVCC.")

# AutoPas Auto4OMP CMake option.
option(AUTOPAS_AUTO4OMP ${AUTOPAS_AUTO4OMP_DOC} ON)

# CMake option to specify the gcc 12 or less path for CUDA's NVCC, in case the default gcc is newer.
option(AUTOPAS_NVCC_GNUC_PATH ${AUTOPAS_NVCC_GNUC_PATH_DOC} "")

if (NOT AUTOPAS_AUTO4OMP)
    # If Auto4OMP disabled, notify.
    message(STATUS "Auto4OMP disabled.")

elseif (NOT AUTOPAS_OPENMP)
    # If OpenMP disabled, warn.
    message(WARNING "OpenMP must be enabled to use Auto4. Auto4OMP disabled.")

else ()
    # Notify.
    message(STATUS "Auto4OMP enabled.")

    # Used the following old policies to suppress Auto4OMP's CMake warnings.
    # TODO: didn't supress the warnings. Why?
    cmake_policy(SET CMP0146 OLD) # FindCUDA.
    cmake_policy(SET CMP0148 OLD) # FindPythonInterp and FindPythonLibs.

    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        # If Clang isn't used, warn. Auto4OMP is tested with up to clang 7 and gcc 4. Modern gcc fails.
        message(WARNING "Auto4OMP needs clang, but ${CMAKE_CXX_COMPILER_ID} is used. Building may produce errors.")
    endif ()

    # Mark as standalone build to fix CMake's missing declaration errors.
    set(OPENMP_STANDALONE_BUILD ON)

    # Enable the FetchContent CMake module.
    include(FetchContent)

    ## Option 1: download Auto4OMP from GIT and make the CMake targets available.
    #FetchContent_Declare(
    #        auto4omp
    #        GIT_REPOSITORY https://github.com/unibas-dmi-hpc/LB4OMP
    #        GIT_TAG v0.1 # The Auto4OMP release.
    #)

    # Option 2: build the pre-packaged Auto4OMP and make the CMake targets available.
    FetchContent_Declare(
            auto4omp
            URL ${AUTOPAS_SOURCE_DIR}/libs/LB4OMP-0.1.zip
            URL_HASH MD5=7bfcf0f896ed99945515a33cd5734cf0 # Calculated with the md5sum command.
    )

    # Auto4OMP's CMake options:

    # Auto4OMP needs a high resolution counter to make measurements; either a cycle counter or time-stamp counter.
    # TODO: processes not executing correctly, outputs are empty.

    # __builtin_readcyclecounter() is available since clang 4 (possibly prior). [1]
    # On systems lacking a cycle counter register or similar, the function defaults to 0.
    # Compile and run CycleCounterQuery.cpp to check if a cycle counter's available, inspired from [2, 3].
    add_executable(has_builtin_readcyclecounter
            EXCLUDE_FROM_ALL "${CMAKE_SOURCE_DIR}/src/autopas/utils/CycleCounterQuery.cpp")
    set_target_properties(has_builtin_readcyclecounter
            PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/scripts")
    execute_process(COMMAND "${CMAKE_BINARY_DIR}/scripts/has_builtin_readcyclecounter"
            OUTPUT_VARIABLE HAS_BUILTIN_READCYCLECOUNTER_OUTPUT)

    if (HAS_BUILTIN_READCYCLECOUNTER_OUTPUT)
        # If cycle counter available, enable.
        message(STATUS "Cycle counter found, enabling.")
        option(LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER
                "Indicates whether the system supports __builtin_readcyclecounter()." ON)
    else ()
        message(STATUS "No cycle counter found. Will look for time-stamp counter next.")
    endif ()

    # __rdtscp() reads from the time-stamp counter, supported by x86 and x86_64.
    execute_process(COMMAND lscpu | grep rdtsc OUTPUT_VARIABLE LSCPU_GREP_RDTSC_OUTPUT)
    if (LSCPU_GREP_RDTSC_OUTPUT)
        # If time-stamp counter available, enable.
        message(STATUS "Time-stamp counter found, enabling.")
        option(LIBOMP_HAVE___RDTSC "Indicates whether the system supports __rdtscp()." ON)
    elseif (NOT HAS_BUILTIN_READCYCLECOUNTER_OUTPUT)
        message(WARNING "No high resolution counter found, Auto4OMP may complain.")
    endif ()

    # CUDA:
    # CMake's FindCUDA is deprecated, but necessary (I think) as Auto4OMP looks for CUDA in the system.
    # For now, use old policy. TODO: for future-proofing, is there a way to use the new policy?
    find_package(CUDA QUIET)
    if (${CUDA_FOUND})
        # If CUDA's installed, Auto4OMP will use it by default. The following fixes build errors.
        #message(STATUS "CUDA found.")

        # GPU compute capability:
        if (${CUDA_VERSION_MAJOR} GREATER 10)
            # CUDA 11+ dropped support for compute 3.5. Use 5.3 instead.
            # List of Nvidia GPUs and their compute capabilities: [4]
            option(LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES "Nvidia GPU's compute capability." 53)
        endif ()

        # GCC: CUDA requires gcc 12 or less. If installed, pass its path to NVCC.
        # According to LB4OMP's docs, it won't be used to produce binaries, only to pass NVCC's compiler checks.
        if (NOT ${AUTOPAS_NVCC_GNUC_PATH} STREQUAL "")
            # If alternate gcc specified, pass it to Auto4OMP.
            option(LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER ${AUTOPAS_NVCC_GNUC_PATH_DOC} ${AUTOPAS_NVCC_GNUC_PATH})
        elseif (DEFINED ENV{__GNUC__})
            if ($ENV{__GNUC__} GREATER AUTOPAS_NVCC_MAX_GCC)
                # If system's gcc is incompatible with NVCC, warn.
                message(WARNING "CUDA enabled, but system's gcc is newer than ${AUTOPAS_NVCC_MAX_GCC}. \
                NVCC may complain. To fix, disable CUDA or install e.g. gcc-${AUTOPAS_NVCC_MAX_GCC} \
                and and pass its path (e.g. \"/usr/bin/gcc-12\") \
                with the CMake option LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER.")
            endif ()
            # If system's gcc is compatible, do nothing.
        else ()
            # If no gcc found, warn.
            message(WARNING "CUDA enabled, but no gcc found. NVCC may complain. \
            To fix, install gcc ${AUTOPAS_NVCC_MAX_GCC} or less.")
        endif ()
    endif ()

    # Mark Auto4OMP's options advanced.
    mark_as_advanced(
            OPENMP_ENABLE_WERROR
            OPENMP_LIBDIR_SUFFIX
            OPENMP_TEST_C_COMPILER
            OPENMP_TEST_CXX_COMPILER
            OPENMP_LLVM_TOOLS_DIR
            OPENMP_LLVM_LIT_EXECUTABLE
            OPENMP_FILECHECK_EXECUTABLE
            LIBOMP_ARCH
            LIBOMP_MIC_ARCH
            LIBOMP_OMP_VERSION
            LIBOMP_LIB_TYPE
            LIBOMP_USE_VERSION_SYMBOLS
            LIBOMP_ENABLE_SHARED
            LIBOMP_FORTRAN_MODULES
            LIBOMP_USE_ADAPTIVE_LOCKS
            LIBOMP_USE_INTERNODE_ALIGNMENT
            LIBOMP_OMPT_SUPPORT
            LIBOMP_OMPT_OPTIONAL
            LIBOMP_STATS
            LIBOMP_USE_DEBUGGER
            LIBOMP_USE_HWLOC
            LIBOMP_HWLOC_INSTALL_DIR
            LIBOMP_CPPFLAGS
            LIBOMP_CFLAGS
            LIBOMP_CXXFLAGS
            LIBOMP_ASMFLAGS
            LIBOMP_LDFLAGS
            LIBOMP_LIBFLAGS
            LIBOMP_FFLAGS
            LIBOMPTARGET_OPENMP_HEADER_FOLDER
            LIBOMPTARGET_OPENMP_HOST_RTL_FOLDER
            LIBOMPTARGET_NVPTX_ENABLE_BCLIB
            LIBOMPTARGET_NVPTX_CUDA_COMPILER
            LIBOMPTARGET_NVPTX_BC_LINKER
            LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER
            LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES
            LIBOMPTARGET_NVPTX_DEBUG
    )

    # Integrate Auto4OMP into the project.
    FetchContent_MakeAvailable(auto4omp)
    include_directories(SYSTEM "${auto4omp_SOURCE_DIR}/runtime/src")
    link_directories(SYSTEM "${auto4omp_SOURCE_DIR}/runtime/src")

    # Reset Policies.
    cmake_policy(SET CMP0146 NEW) # FindCUDA.
    cmake_policy(SET CMP0148 NEW) # FindPythonInterp and FindPythonLibs.

endif ()

# [1] https://releases.llvm.org/4.0.1/tools/clang/LanguageExtensions.html
# [2] https://stackoverflow.com/a/41380220
# [3] https://stackoverflow.com/a/50306091
# [4] https://developer.nvidia.com/cuda-gpus
# Some code was inspired from autopas_spdlog.cmake and other AutoPas CMake files.
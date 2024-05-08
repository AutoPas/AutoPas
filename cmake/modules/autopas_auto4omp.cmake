# Version settings for convenience. If Auto4OMP supports new packages, update here.
set(AUTOPAS_NVCC_MAX_GCC 12 CACHE STRING "Newest gcc version supported by NVCC.")

# Option descriptions:
set(AUTOPAS_AUTO4OMP_DOC
        "Activates Auto4OMP to automatically select OpenMP scheduling algorithms \
        from the LB4OMP portfolio during runtime.")

set(AUTOPAS_LLVM_LIT_EXECUTABLE_DOC
        "Path to the LLVM-lit executable, required by Auto4OMP. E.g., \"/usr/lib/llvm-18/build/utils/lit/lit.py\". \
        To find it, run the terminal command \"find /usr -name lit.py\". \
        By default, Auto4OMP looks for LLVM-lit with the PATH environment variable.")

set(AUTOPAS_FILECHECK_EXECUTABLE_DOC
        "Path to the FileCheck executable, required by Auto4OMP. E.g., \"/usr/lib/llvm-18/bin/FileCheck\". \
        To find it, run the terminal command \"find /usr -name lit.py\". \
        By default, Auto4OMP looks for FileCheck with the PATH environment variable.")

set(AUTOPAS_NVCC_GNUC_PATH_DOC
        "Path to gcc ${AUTOPAS_NVCC_MAX_GCC} or less, required by NVCC, in case system's gcc is newer.")

set(LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER_DOC
        "Indicates whether the system supports __builtin_readcyclecounter().")

set(LIBOMP_HAVE___RDTSC_DOC
        "Indicates whether the system supports __rdtscp().")

set(LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES_DOC
        "Nvidia GPU's compute capability.")

# AutoPas options:
option(AUTOPAS_AUTO4OMP ${AUTOPAS_AUTO4OMP_DOC} ON)
set(AUTOPAS_NVCC_GNUC_PATH "" CACHE FILEPATH ${AUTOPAS_NVCC_GNUC_PATH_DOC})
set(AUTOPAS_LLVM_LIT_EXECUTABLE "" CACHE FILEPATH ${AUTOPAS_LLVM_LIT_EXECUTABLE_DOC})
set(AUTOPAS_FILECHECK_EXECUTABLE "" CACHE FILEPATH ${AUTOPAS_FILECHECK_EXECUTABLE_DOC})

if (NOT AUTOPAS_AUTO4OMP)
    # If Auto4OMP disabled, notify.
    message(STATUS "Auto4OMP disabled.")

elseif (NOT AUTOPAS_OPENMP)
    # If OpenMP disabled, warn.
    message(WARNING "OpenMP must be enabled to use Auto4. Auto4OMP disabled.")

else ()
    # Notify.
    message(STATUS "Auto4OMP enabled.")

    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        # If Clang isn't used, warn. Auto4OMP is tested with up to clang 7 and gcc 4. Modern gcc fails.
        message(WARNING "Auto4OMP needs clang, but ${CMAKE_CXX_COMPILER_ID} is used. Building may produce errors.")
    endif ()

    # Auto4OMP's CMake options:

    ## High resolution timer:
    ### Auto4OMP needs a high resolution timer to make measurements; either a cycle counter or time-stamp counter.

    ### Cycle counter:
    #### __builtin_readcyclecounter() is available since clang 4 (possibly prior). [1]
    #### On systems lacking a cycle counter register or similar, the function defaults to 0.
    #### Compile and run CycleCounterQuery.cpp to check if a cycle counter's available, inspired from [2, 3].
    add_executable(has_builtin_readcyclecounter
            EXCLUDE_FROM_ALL "${CMAKE_SOURCE_DIR}/src/autopas/utils/CycleCounterQuery.cpp")
    set_target_properties(has_builtin_readcyclecounter
            PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/scripts")

    #### TODO: when building target with a cmake command, cmake fails to load its cache as it's not done generating.
    ####execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target has_builtin_readcyclecounter)

    #### For now, build manually.
    execute_process(
            COMMAND mkdir -p ${CMAKE_BINARY_DIR}/scripts
            COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_SOURCE_DIR}/src/autopas/utils/CycleCounterQuery.cpp
            -o ${CMAKE_BINARY_DIR}/scripts/has_builtin_readcyclecounter
    )
    execute_process(
            COMMAND "${CMAKE_BINARY_DIR}/scripts/has_builtin_readcyclecounter"
            OUTPUT_VARIABLE HAS_BUILTIN_READCYCLECOUNTER_OUTPUT
    )

    #### TODO: (extra) why doesn't if (string) return true when string's not empty?
    if (${HAS_BUILTIN_READCYCLECOUNTER_OUTPUT} GREATER 0)
        # __has_builtin() returns  a non-zero int if the feature's supported, else 0. [4, 5]
        # If cycle counter available, enable.
        message(STATUS "Cycle counter found, enabling.")
        option(LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER ${LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER_DOC} ON)
    else ()
        message(STATUS "No cycle counter found. Will look for time-stamp counter next.")
    endif ()

    ### Time-stamp counter:
    #### __rdtscp() reads from the time-stamp counter, supported by x86 and x86_64. Inspired from [6]
    execute_process(COMMAND lscpu COMMAND grep rdtsc OUTPUT_VARIABLE LSCPU_GREP_RDTSC_OUTPUT)
    if (NOT ${LSCPU_GREP_RDTSC_OUTPUT} STREQUAL "")
        # If time-stamp counter available, enable.
        message(STATUS "Time-stamp counter found, enabling.")
        option(LIBOMP_HAVE___RDTSC ${LIBOMP_HAVE___RDTSC_DOC} ON)
    elseif (${HAS_BUILTIN_READCYCLECOUNTER_OUTPUT} LESS_EQUAL 0)
        message(STATUS "No time-stamp counter found.")
        message(WARNING "No high resolution timer found, Auto4OMP may complain. \n
        Auto4OMP requires one of two high resolution timers: a cycle counter or time-stamp counter. \n
        The cycle counter register is read with the builtin function __builtin_readcyclecounter(). \
        To check if the system provides this timer, call __has_builtin(__builtin_readcyclecounter). \n
        The time-stamp counter is read with __rdtscp(). \
        To check if the system provides it, run the following command: lscpu | grep rdtsc \n
        Cmake attempted to query for the two timers but did not find them. \
        Auto4OMP will thus likely produce build errors.")
    endif ()

    ## CUDA:
    block(SCOPE_FOR POLICIES)
        # CMake's FindCUDA is deprecated, but necessary (I think) as Auto4OMP looks for CUDA in the system.
        # For now, use old policy. TODO: for future-proofing, is there a way to use the new policy?
        cmake_policy(SET CMP0146 OLD)

        find_package(CUDA QUIET)
        if (${CUDA_FOUND})
            # If CUDA's installed, Auto4OMP will use it by default. The following fixes build errors.
            message(STATUS "CUDA found.")

            # GPU compute capability:
            if (${CUDA_VERSION_MAJOR} GREATER 10)
                # CUDA 11+ dropped support for compute 3.5. Use 5.3 instead.
                # List of Nvidia GPUs and their compute capabilities: [7]
                set(LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES
                        53 CACHE STRING ${LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES_DOC} FORCE)
            endif ()

            # GCC: CUDA requires gcc 12 or less. If installed, pass its path to NVCC.
            # According to LB4OMP's docs, it won't be used to produce binaries, only to pass NVCC's compiler checks.
            if (NOT ${AUTOPAS_NVCC_GNUC_PATH} STREQUAL "")
                # If alternate gcc specified, pass it to Auto4OMP.
                set(LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER
                        ${AUTOPAS_NVCC_GNUC_PATH} CACHE STRING ${AUTOPAS_NVCC_GNUC_PATH_DOC} FORCE)
                message(STATUS "GCC: ${LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER}") # Debug.
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
    endblock()

    ## LLVM-lit:
    ### Look for FileCheck with the PATH environment variable.
    if (NOT AUTOPAS_LLVM_LIT_EXECUTABLE)
        message(DEBUG "Looking for LLVM-lit with the PATH environment variable.")
        find_program(AUTOPAS_LLVM_LIT_EXECUTABLE NAMES llvm-lit lit lit.py)
    endif ()

    ### Look for FileCheck with the find command.
    if (NOT AUTOPAS_LLVM_LIT_EXECUTABLE)
        message(DEBUG "Looking for LLVM-lit with the find command.")
        execute_process(
                COMMAND find /usr/lib -name lit.py # Lists paths containing lit.
                COMMAND sort -dfr # Puts paths through newer LLVM versions at the top.
                COMMAND grep -m1 "" #Keeps first path only, inspired from [8].
                COMMAND tr -d '\n' # Truncates newlines.
                OUTPUT_VARIABLE AUTOPAS_LLVM_LIT_EXECUTABLE
        )
    endif ()

    ### If a path to LLVM-lit executable was specified, pass it to Auto4OMP.
    if (AUTOPAS_LLVM_LIT_EXECUTABLE)
        message(STATUS "LLVM-lit executable at ${AUTOPAS_LLVM_LIT_EXECUTABLE}")
        set(
                OPENMP_LLVM_LIT_EXECUTABLE
                ${AUTOPAS_LLVM_LIT_EXECUTABLE}
                CACHE FILEPATH ${AUTOPAS_LLVM_LIT_EXECUTABLE_DOC} FORCE
        )
    else ()
        message(STATUS "No path to LLVM-lit was found, Auto4OMP may warn. \
        AutoPas attempted to look for lit with the PATH environment variable \
        and the command \"find /usr/lib -name lit.py\", but did not find it. \
        To fix, specify the path with -DAUTOPAS_LLVM_LIT_EXECUTABLE=\"path/to/lit.py\"")
    endif ()

    ## FileCheck:
    ### Look for FileCheck with the PATH environment variable.
    if (NOT AUTOPAS_FILECHECK_EXECUTABLE)
        message(DEBUG "Looking for FileCheck with the PATH environment variable.")
        find_program(AUTOPAS_FILECHECK_EXECUTABLE NAMES FileCheck)
    endif ()

    ### Look for FileCheck with the find command.
    if (NOT AUTOPAS_FILECHECK_EXECUTABLE)
        message(DEBUG "Looking for FileCheck with the find command.")
        execute_process(
                COMMAND find /usr/lib -name FileCheck # Lists paths containing FileCheck.
                COMMAND sort -dfr # Puts paths through newer LLVM versions at the top.
                COMMAND grep -m1 "" # Keeps first path only, inspired from [8].
                COMMAND tr -d '\n' # Truncates newlines.
                OUTPUT_VARIABLE AUTOPAS_FILECHECK_EXECUTABLE
        )
    endif ()

    ### If a path to FileCheck executable was found or specified, pass it to Auto4OMP.
    if (AUTOPAS_FILECHECK_EXECUTABLE)
        message(STATUS "FileCheck executable at ${AUTOPAS_FILECHECK_EXECUTABLE}")
        set(
                OPENMP_FILECHECK_EXECUTABLE
                ${AUTOPAS_FILECHECK_EXECUTABLE}
                CACHE FILEPATH ${AUTOPAS_FILECHECK_EXECUTABLE_DOC} FORCE
        )
    else ()
        message(STATUS "No path to FileCheck was found, Auto4OMP may warn. \
        AutoPas attempted to look for FileCheck with the PATH environment variable \
        and the command \"find /usr/lib -name FileCheck\", but did not find it. \
        To fix, specify the path with -DAUTOPAS_LLVM_LIT_EXECUTABLE=\"path/to/FileCheck\"")
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

    #set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
    #set(CMAKE_POLICY_DEFAULT_CMP0126 NEW)

    # Mark as standalone build to fix CMake's missing declaration errors.
    set(OPENMP_STANDALONE_BUILD ON)

    # Enable the FetchContent CMake module.
    include(FetchContent)

    ## Option 1: download Auto4OMP from GIT and make the CMake targets available. TODO: untested.
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
            CMAKE_ARGS "-DCMAKE_C_COMPILER=clang;-DCMAKE_CXX_COMPILER=clang++"
    )

    # Use the following old policies to suppress Auto4OMP's CMake warnings.
    # cmake_policy() doesn't work, likely because Auto4OMP isn't in a subdirectory of this module's directory.
    set(CMAKE_POLICY_DEFAULT_CMP0146 OLD) # The FindCUDA module.
    set(CMAKE_POLICY_DEFAULT_CMP0148 OLD) # The FindPythonInterp and FindPythonLibs modules.

    # Integrate Auto4OMP into the project.
    FetchContent_MakeAvailable(auto4omp)
    include_directories(SYSTEM "${auto4omp_SOURCE_DIR}/runtime/src")
    link_directories(SYSTEM "${auto4omp_SOURCE_DIR}/runtime/src")

    # Reset Policies.
    set(CMAKE_POLICY_DEFAULT_CMP0146 NEW) # The FindCUDA module.
    set(CMAKE_POLICY_DEFAULT_CMP0148 NEW) # The FindPythonInterp and FindPythonLibs modules.

endif ()

# [1] https://releases.llvm.org/4.0.1/tools/clang/LanguageExtensions.html
# [2] https://stackoverflow.com/a/41380220
# [3] https://stackoverflow.com/a/50306091
# [4] https://clang.llvm.org/docs/LanguageExtensions.html#has-builtin
# [5] https://gcc.gnu.org/onlinedocs/cpp/_005f_005fhas_005fbuiltin.html
# [6] https://stackoverflow.com/a/35693504
# [7] https://developer.nvidia.com/cuda-gpus
# [8] https://unix.stackexchange.com/a/294493
# Some code was inspired from autopas_spdlog.cmake and other AutoPas CMake files.
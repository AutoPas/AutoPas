#[=====================================================================================================================[
File: autopas_auto4omp.cmake
Author: MehdiHachicha
Date: 15.04.2024

This CMake module loads Auto4OMP into the AutoPas project.
It attempts to automatically handle its required CMake arguments, and warns if one must be manually entered.
Ideally, "cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++" should work without issues.
If cuda's installed, a gcc compatible with NVCC should also be installed. E.g., gcc-12.
#]=====================================================================================================================]

# Variable descriptions:
string(
        CONCAT AUTOPAS_AUTO4OMP_DOC
        "Activates Auto4OMP to automatically select OpenMP scheduling algorithms "
        "from the LB4OMP portfolio during runtime."
)

set(
        AUTOPAS_AUTO4OMP_GIT_DOC
        "Downloads Auto4OMP from Git instead of the bundled zip."
)

set(
        AUTOPAS_AUTO4OMP_GIT_FORCE_DOC
        "Downloads Auto4OMP from Git despite potential incompatibility with AutoPas."
)

set(
        AUTOPAS_AUTO4OMP_GIT_TAG_DOC
        "The Auto4OMP Git branch or commit to use."
)

string(
        CONCAT AUTOPAS_AUTO4OMP_DEBUG_DOC
        "Uses a custom bundled Auto4OMP that prints a message "
        "to confirm that Auto4OMP is indeed being used instead of standard OpenMP."
)

set(
        AUTOPAS_AUTO4OMP_INSTALL_DIR_DOC
        "Directory for Auto4OMP to install its libraries."
)

string(
        CONCAT AUTOPAS_LLVM_LIT_EXECUTABLE_DOC
        "Path to the LLVM-lit executable, required by Auto4OMP. E.g., \"/usr/lib/llvm-18/build/utils/lit/lit.py\". "
        "To find it, run the following command: find /usr -name lit.py. "
        "By default, Auto4OMP looks for LLVM-lit with the PATH environment variable."
)

string(
        CONCAT AUTOPAS_FILECHECK_EXECUTABLE_DOC
        "Path to the FileCheck executable, required by Auto4OMP. E.g., \"/usr/lib/llvm-18/bin/FileCheck\". "
        "To find it, run the following command: \"find /usr -name lit.py\". "
        "By default, Auto4OMP looks for FileCheck with the PATH environment variable."
)

set(
        AUTOPAS_NVCC_GNUC_PATH_DOC
        "Path to gcc ${AUTOPAS_NVCC_MAX_GCC} or less, required by NVCC, in case system's gcc is newer."
)

string(
        CONCAT AUTOPAS_pthread_LIBRARY
        "Path to pthread.a, required by OpenMP. E.g., \"/usr/lib/x86_64-linux-gnu/libpthread.a\". "
        "To find it, run the following command: \"find /usr -name libpthread.a\". "
        "By default, AutoPas looks for pthread with FindThreads and the command \"find /usr/lib -name FileCheck\"."
)

set(
        LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER_DOC
        "Indicates whether the system supports __builtin_readcyclecounter()."
)

set(
        LIBOMP_HAVE___RDTSC_DOC
        "Indicates whether the system supports __rdtscp()."
)

set(
        LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES_DOC
        "Nvidia GPU's compute capability."
)

set(
        AUTOPAS_NVCC_MAX_GCC_DOC
        "Newest gcc version supported by NVCC."
)

set(
        AUTOPAS_OMP_VERSION_DOC
        "Newest OpenMP version supported by AutoPas."
)

function(set_autopas_auto4omp_target_boolean_doc bool_var target)
    set(${bool_var} "Boolean that specifies whether Auto4OMP's target ${target} is defined." PARENT_SCOPE)
endfunction()
set_autopas_auto4omp_target_boolean_doc("AUTOPAS_TARGET_omp_DOC" "omp")
set_autopas_auto4omp_target_boolean_doc("AUTOPAS_TARGET_omptarget_DOC" "omptarget")
set_autopas_auto4omp_target_boolean_doc("AUTOPAS_TARGET_omptarget.rtl.cuda_DOC" "omptarget.rtl.cuda")
set_autopas_auto4omp_target_boolean_doc("AUTOPAS_TARGET_omptarget.rtl.x86_64_DOC" "omptarget.rtl.x86_64")
set_autopas_auto4omp_target_boolean_doc("AUTOPAS_TARGET_omptarget-nvptx_DOC" "omptarget-nvptx")

# AutoPas CMake variables for Auto4OMP.
option(AUTOPAS_AUTO4OMP "${AUTOPAS_AUTO4OMP_DOC}" ON)
option(AUTOPAS_AUTO4OMP_GIT "${AUTOPAS_AUTO4OMP_GIT_DOC}" OFF)
option(AUTOPAS_AUTO4OMP_GIT_FORCE "${AUTOPAS_AUTO4OMP_GIT_FORCE_DOC}" OFF)
set(AUTOPAS_AUTO4OMP_GIT_TAG "master" CACHE STRING "${AUTOPAS_AUTO4OMP_GIT_TAG_DOC}")
option(AUTOPAS_AUTO4OMP_DEBUG "${AUTOPAS_AUTO4OMP_DEBUG_DOC}" OFF)
set(AUTOPAS_AUTO4OMP_INSTALL_DIR "" CACHE PATH ${AUTOPAS_AUTO4OMP_INSTALL_DIR_DOC})
set(AUTOPAS_NVCC_GNUC_PATH OFF CACHE FILEPATH "${AUTOPAS_NVCC_GNUC_PATH_DOC}")
set(AUTOPAS_LLVM_LIT_EXECUTABLE OFF CACHE FILEPATH "${AUTOPAS_LLVM_LIT_EXECUTABLE_DOC}")
set(AUTOPAS_FILECHECK_EXECUTABLE OFF CACHE FILEPATH "${AUTOPAS_FILECHECK_EXECUTABLE_DOC}")
set(AUTOPAS_pthread_LIBRARY OFF CACHE FILEPATH "${AUTOPAS_FILECHECK_EXECUTABLE_DOC}")

# Version settings. If AutoPas or Auto4OMP support new packages, update here.
set(AUTOPAS_NVCC_MAX_GCC 12 CACHE STRING ${AUTOPAS_NVCC_MAX_GCC_DOC})
set(AUTOPAS_OMP_VERSION 45 CACHE STRING ${AUTOPAS_OMP_VERSION_DOC})

# Set boolean variables for Auto4OMP's targets, used in target_link_libraries's conditional generator expressions.
set(AUTOPAS_TARGET_omp OFF CACHE BOOL "${AUTOPAS_TARGET_omp_DOC}")
set(AUTOPAS_TARGET_omptarget OFF CACHE BOOL "${AUTOPAS_TARGET_omptarget_DOC}")
set(AUTOPAS_TARGET_omptarget.rtl.cuda OFF CACHE BOOL "${AUTOPAS_TARGET_omptarget.rtl.cuda_DOC}")
set(AUTOPAS_TARGET_omptarget.rtl.x86_64 OFF CACHE BOOL "${AUTOPAS_TARGET_omptarget.rtl.x86_64_DOC}")
set(AUTOPAS_TARGET_omptarget-nvptx OFF CACHE BOOL "${AUTOPAS_TARGET_omptarget-nvptx_DOC}")

if (NOT AUTOPAS_AUTO4OMP)
    # If Auto4OMP disabled, notify.
    message(STATUS "Auto4OMP disabled.")

elseif (NOT AUTOPAS_OPENMP)
    # If OpenMP disabled, warn.
    message(WARNING "OpenMP must be enabled to use Auto4. Auto4OMP disabled.")
    set(AUTOPAS_AUTO4OMP OFF CACHE BOOL ${AUTOPAS_AUTO4OMP_DOC} FORCE)

else ()
    # Notify.
    message(STATUS "Auto4OMP enabled.")

    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        # If Clang isn't used, warn. Auto4OMP is tested with up to clang 7 and gcc 4. Modern gcc fails.
        message(WARNING "Auto4OMP needs clang, but ${CMAKE_CXX_COMPILER_ID} is used. Building may produce errors.")
    endif ()

    # Auto4OMP's CMake variables:

    ## OpenMP version:
    ### AutoPas supports up to OpenMP 4.5. Pass this to Auto4OMP.
    set(LIBOMP_OMP_VERSION ${AUTOPAS_OMP_VERSION} CACHE INTERNAL "OpenMP version for Auto4OMP to build for.")

    ### OpenMP tools is only supported from OpenMP 5.
    if (${AUTOPAS_OMP_VERSION} GREATER_EQUAL 50)
        set(LIBOMP_OMPT_SUPPORT ON CACHE INTERNAL "OpenMP tools support.")
    else ()
        set(LIBOMP_OMPT_SUPPORT OFF CACHE INTERNAL "OpenMP tools support.")
    endif ()

    ## High resolution timer:
    ### Auto4OMP needs a high resolution timer to make measurements; either a cycle counter or time-stamp counter.

    ### Cycle counter:
    #### __builtin_readcyclecounter() is available since clang 4 (possibly prior). [1]
    #### On systems lacking a cycle counter register or similar, the function defaults to 0.
    #### Compile and run CycleCounterQuery.cpp to check if a cycle counter's available, inspired from [2, 3].
    add_executable(
            has_builtin_readcyclecounter
            EXCLUDE_FROM_ALL "${CMAKE_SOURCE_DIR}/src/autopas/utils/CycleCounterQuery.cpp"
    )
    set_target_properties(
            has_builtin_readcyclecounter
            PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/scripts"
    )

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

    if (${HAS_BUILTIN_READCYCLECOUNTER_OUTPUT} GREATER 0)
        # __has_builtin() returns  a non-zero int if the feature's supported, else 0. [4, 5]
        # If cycle counter available, enable.
        message(DEBUG "Cycle counter found, enabling.")
        option(LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER "${LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER_DOC}" ON)
    else ()
        message(STATUS "No cycle counter found. Will look for time-stamp counter next.")
    endif ()

    ### Time-stamp counter:
    #### __rdtscp() reads from the time-stamp counter, supported by x86 and x86_64. Inspired from [6]
    execute_process(COMMAND lscpu COMMAND grep rdtsc OUTPUT_VARIABLE LSCPU_GREP_RDTSC_OUTPUT)
    if (NOT "${LSCPU_GREP_RDTSC_OUTPUT}" STREQUAL "")
        # If time-stamp counter available, enable.
        message(DEBUG "Time-stamp counter found, enabling.")
        option(LIBOMP_HAVE___RDTSC ${LIBOMP_HAVE___RDTSC_DOC} ON)
    elseif (${HAS_BUILTIN_READCYCLECOUNTER_OUTPUT} LESS_EQUAL 0)
        message(STATUS "No time-stamp counter found.")
        message(
                WARNING
                "No high resolution timer found, Auto4OMP may complain.\n"
                "Auto4OMP requires one of two high resolution timers: a cycle counter or time-stamp counter.\n"
                "The cycle counter register is read with the builtin function __builtin_readcyclecounter(). "
                "To check if the system provides this timer, call __has_builtin(__builtin_readcyclecounter).\n"
                "The time-stamp counter is read with __rdtscp(). "
                "To check if the system provides it, run the following command: lscpu | grep rdtsc\n"
                "Cmake attempted to query for the two timers but did not find them. Auto4OMP may fail."
        )
    endif ()

    ## CUDA:
    ### CMake's FindCUDA is deprecated, but necessary (I think) as Auto4OMP looks for CUDA in the system.
    ### For now, use old policy. TODO: for future-proofing, is there a way to use the new policy?
    cmake_policy(PUSH) # Use PUSH and POP instead of the new block() to preserve compatibility with older CMake.
    if (POLICY CMP0146)
        cmake_policy(SET CMP0146 OLD)
    endif ()

    find_package(CUDA QUIET)
    if (${CUDA_FOUND})
        # If CUDA's installed, Auto4OMP will use it by default. The following fixes build errors.
        message(DEBUG "CUDA found.")

        # GPU compute capability:
        if (${CUDA_VERSION_MAJOR} GREATER 10)
            # CUDA 11+ dropped support for compute 3.5. Use 5.3 instead.
            # List of Nvidia GPUs and their compute capabilities: [7]
            set(
                    LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES 53
                    CACHE INTERNAL "${LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES_DOC}"
            )
        endif ()

        # GCC: CUDA requires gcc 12 or less. If installed, pass its path to NVCC.
        # According to LB4OMP's docs, it won't be used to produce binaries, only to pass NVCC's compiler checks.
        if (NOT AUTOPAS_NVCC_GNUC_PATH)
            find_program(GCC_PATH gcc)
            if (GCC_PATH)
                # Isolate gcc's version major from gcc's version output.
                set(GCC_VERSION_MAJOR OFF)
                execute_process(
                        COMMAND ${GCC_PATH} --version # Multi-line output.
                        COMMAND grep gcc # First line, including the version.
                        COMMAND cut -d " " -f4 # Full version, format ?.?.?.
                        COMMAND cut -d "." -f1 # Version major.
                        COMMAND tr -d '\n' # Truncates newlines.
                        OUTPUT_VARIABLE GCC_VERSION_MAJOR
                )
                if (${GCC_VERSION_MAJOR} LESS_EQUAL ${AUTOPAS_NVCC_MAX_GCC})
                    set(AUTOPAS_NVCC_GNUC_PATH "${GCC_PATH}" CACHE FILEPATH "${AUTOPAS_NVCC_GNUC_PATH_DOC}" FORCE)
                else ()
                    message(
                            STATUS
                            "CUDA enabled, but system's gcc is incompatible with NVCC. Looking for compatible gcc."
                    )
                endif ()
            endif ()
        endif ()
        if (NOT AUTOPAS_NVCC_GNUC_PATH)
            # If system's gcc is incompatible with NVCC, look for compatible gcc.
            # Look for gcc executable names based on Ubuntu's packages [8].
            # Names based on other package managers can be added here too.
            find_program(GCC_LESS NAMES gcc-12 gcc-11 gcc-10 gcc-9 gcc-8 gcc-7)
            if (GCC_LESS)
                set(AUTOPAS_NVCC_GNUC_PATH "${GCC_LESS}" CACHE FILEPATH "${AUTOPAS_NVCC_GNUC_PATH_DOC}" FORCE)
            else ()
                message(
                        WARNING
                        "CUDA enabled, but system's gcc is newer than ${AUTOPAS_NVCC_MAX_GCC}. NVCC may complain.\n"
                        "To fix, disable CUDA or install e.g. gcc-${AUTOPAS_NVCC_MAX_GCC} "
                        "and and pass its path (e.g. \"/usr/bin/gcc-${AUTOPAS_NVCC_MAX_GCC}\") "
                        "with the CMake variable LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER."
                )
            endif ()
        endif ()

        if (AUTOPAS_NVCC_GNUC_PATH)
            # If alternate gcc found or specified, pass it to Auto4OMP.
            set(
                    LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER "${AUTOPAS_NVCC_GNUC_PATH}"
                    CACHE INTERNAL "${AUTOPAS_NVCC_GNUC_PATH_DOC}"
            )
            message(DEBUG "Alternate GCC for NVCC at ${LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER}")
        else ()
            #If no gcc found, warn.
            message(
                    WARNING
                    "CUDA enabled, but no gcc found. NVCC may complain during Auto4OMP's build."
                    "To fix, install gcc ${AUTOPAS_NVCC_MAX_GCC} or less."
            )
        endif ()
    endif ()
    cmake_policy(POP)

    ## LLVM-lit:
    ### Look for FileCheck with the PATH environment variable.
    if (NOT AUTOPAS_LLVM_LIT_EXECUTABLE)
        message(DEBUG "Looking for LLVM-lit with the PATH environment variable.")
        find_program(AUTOPAS_LLVM_LIT_EXECUTABLE NAMES lit llvm-lit lit.py)
    endif ()

    ### Look for FileCheck with the find command.
    if (NOT AUTOPAS_LLVM_LIT_EXECUTABLE)
        message(DEBUG "Looking for LLVM-lit with the find command.")
        execute_process(
                COMMAND find /usr/lib -name lit.py # Lists paths containing lit.
                COMMAND sort -dfr # Puts paths through newer LLVM versions at the top.
                COMMAND grep -m1 "" # Keeps first path only, inspired from [9].
                COMMAND tr -d '\n' # Truncates newlines.
                OUTPUT_VARIABLE AUTOPAS_LLVM_LIT_EXECUTABLE
        )
    endif ()

    ### If a path to LLVM-lit executable was specified, pass it to Auto4OMP.
    if (AUTOPAS_LLVM_LIT_EXECUTABLE)
        message(DEBUG "LLVM-lit executable at ${AUTOPAS_LLVM_LIT_EXECUTABLE}")
        set(
                OPENMP_LLVM_LIT_EXECUTABLE "${AUTOPAS_LLVM_LIT_EXECUTABLE}"
                CACHE INTERNAL "${AUTOPAS_LLVM_LIT_EXECUTABLE_DOC}"
        )
    else ()
        message(
                STATUS
                "No path to LLVM-lit was found, Auto4OMP may warn.\n"
                "AutoPas attempted to look for lit with the PATH environment variable "
                "and the command \"find /usr/lib -name lit.py\", but did not find it.\n"
                "To fix, specify the path with -DAUTOPAS_LLVM_LIT_EXECUTABLE=\"path/to/lit.py\""
        )
    endif ()

    ## FileCheck:
    ### Look for FileCheck with the PATH environment variable.
    if (NOT AUTOPAS_FILECHECK_EXECUTABLE)
        message(DEBUG "Looking for FileCheck with the PATH environment variable.")
        find_program(AUTOPAS_FILECHECK_EXECUTABLE FileCheck)
    endif ()

    ### Look for FileCheck with the find command.
    if (NOT AUTOPAS_FILECHECK_EXECUTABLE)
        message(DEBUG "Looking for FileCheck with the find command.")
        execute_process(
                COMMAND find /usr/lib -name FileCheck # Lists paths containing FileCheck.
                COMMAND sort -dfr # Puts paths through newer LLVM versions at the top.
                COMMAND grep -m1 "" # Keeps first path only, inspired from [9].
                COMMAND tr -d '\n' # Truncates newlines.
                OUTPUT_VARIABLE AUTOPAS_FILECHECK_EXECUTABLE
        )
    endif ()

    ### If a path to FileCheck executable was found or specified, pass it to Auto4OMP.
    if (AUTOPAS_FILECHECK_EXECUTABLE)
        message(DEBUG "FileCheck executable at ${AUTOPAS_FILECHECK_EXECUTABLE}")
        set(
                OPENMP_FILECHECK_EXECUTABLE "${AUTOPAS_FILECHECK_EXECUTABLE}"
                CACHE INTERNAL "${AUTOPAS_FILECHECK_EXECUTABLE_DOC}"
        )
    else ()
        message(
                STATUS
                "No path to FileCheck was found, Auto4OMP may warn.\n"
                "AutoPas attempted to look for FileCheck with the PATH environment variable "
                "and the command \"find /usr/lib -name FileCheck\", but did not find it.\n"
                "To fix, specify the path with -DAUTOPAS_LLVM_LIT_EXECUTABLE=\"path/to/FileCheck\""
        )
    endif ()

    # Mark as standalone build to fix CMake's missing declaration errors.
    set(OPENMP_STANDALONE_BUILD ON)

    # Enable the FetchContent CMake module.
    include(FetchContent)

    # If OMP 4.5 is to be used, use the custom Auto4OMP zip to fix resulting build errors.
    if (AUTOPAS_AUTO4OMP_DEBUG)
        # Build the pre-packaged Auto4OMP that prints a message from kmp.h at runtime to confirm its use.
        FetchContent_Declare(
                auto4omp
                URL ${AUTOPAS_SOURCE_DIR}/libs/LB4OMP-debug.zip # [*]
                URL_HASH MD5=b1e85a5851894f99e290864e1a35c614 # Calculated with the md5sum command.
        )
    elseif (${AUTOPAS_OMP_VERSION} LESS 50 AND NOT AUTOPAS_AUTO4OMP_GIT_FORCE)
        # Build the pre-packaged Auto4OMP including the fix for build errors when specifying an old OMP version.
        FetchContent_Declare(
                auto4omp
                URL ${AUTOPAS_SOURCE_DIR}/libs/LB4OMP-custom.zip # [*]
                URL_HASH MD5=c62773279f8c73e72a6d1bcb36334fef # Calculated with the md5sum command.
        )
    elseif (AUTOPAS_AUTO4OMP_GIT OR AUTOPAS_AUTO4OMP_GIT_FORCE)
        # Download Auto4OMP from Git.
        FetchContent_Declare(
                auto4omp
                GIT_REPOSITORY https://github.com/unibas-dmi-hpc/LB4OMP # [*]
                GIT_TAG ${AUTOPAS_AUTO4OMP_GIT_TAG}
        )
        # If the custom zip must be used, warn.
        if (AUTOPAS_AUTO4OMP_GIT_FORCE AND ${AUTOPAS_OMP_VERSION} LESS 50)
            message(
                    WARNING
                    "OpenMP <5 is used, which will potentially result in build errors with Auto4OMP."
                    "Unless Auto4OMP releases a fix, it's not advised force a Git download."
                    "The bundled libs/LB4OMP-custom.zip contains a fix for this issue."
            )
        endif ()
    else ()
        # Build the pre-packaged Auto4OMP and make the CMake targets available. [***]
        ## Targets include omp, omptarget, ...
        FetchContent_Declare(
                auto4omp
                URL ${AUTOPAS_SOURCE_DIR}/libs/LB4OMP-master.zip # [*]
                URL_HASH MD5=b8090fa162eb2232870916271f36650f # Calculated with the md5sum command.
        )
    endif ()

    # Use the following old policies to suppress Auto4OMP's CMake warnings.
    ## cmake_policy() doesn't work, likely because Auto4OMP isn't in a subdirectory of this module's directory.
    set(CMAKE_POLICY_DEFAULT_CMP0146 OLD) # The FindCUDA module.
    set(CMAKE_POLICY_DEFAULT_CMP0148 OLD) # The FindPythonInterp and FindPythonLibs modules.

    # Integrate Auto4OMP into the project.
    FetchContent_MakeAvailable(auto4omp)

    # Reset Policies.
    set(CMAKE_POLICY_DEFAULT_CMP0146 NEW) # The FindCUDA module.
    set(CMAKE_POLICY_DEFAULT_CMP0148 NEW) # The FindPythonInterp and FindPythonLibs modules.

    # Mark Auto4OMP's variables advanced. [***]
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

    # Build and install Auto4OMP with make install.

    ## Add a target for Auto4OMP. It builds Auto4OMP by running make in its build directory.
    add_custom_target(
            auto4omp ALL
            WORKING_DIRECTORY ${auto4omp_BINARY_DIR}
            COMMAND ${CMAKE_MAKE_PROGRAM}
            COMMENT "Building Auto4OMP"
    )

    ## Set Auto4OMP's install directory if not specified.
    if (NOT AUTOPAS_AUTO4OMP_INSTALL_DIR)
        set(
                AUTOPAS_AUTO4OMP_INSTALL_DIR "${auto4omp_BINARY_DIR}/install"
                CACHE STRING ${AUTOPAS_AUTO4OMP_INSTALL_DIR_DOC} FORCE
        )
    endif ()
    set(AUTOPAS_AUTO4OMP_LIB_DIR "${AUTOPAS_AUTO4OMP_INSTALL_DIR}/lib")
    set(AUTOPAS_AUTO4OMP_INCLUDE_DIR "${AUTOPAS_AUTO4OMP_INSTALL_DIR}/include")

    ## Configure cmake install to use the local install directory.
    install(DIRECTORY "${auto4omp_BINARY_DIR}/runtime/src" DESTINATION "${AUTOPAS_AUTO4OMP_INSTALL_DIR}")

    ## Run cmake install post-build. Use a local directory instead of /usr, so no sudo privilege is required.
    add_custom_command(
            TARGET auto4omp POST_BUILD
            COMMAND ${CMAKE_COMMAND} --install "${auto4omp_BINARY_DIR}" --prefix "${AUTOPAS_AUTO4OMP_INSTALL_DIR}"
    )

    # Set booleans to indicate which Auto4OMP targets are defined, add custom targets for them. [11]

    ## libomp.so
    set(
            AUTOPAS_TARGET_omp ON
            CACHE BOOL "${AUTOPAS_TARGET_omp_DOC}" FORCE
    )
    add_library(autopas_auto4omp_omp SHARED IMPORTED)
    set_property(
            TARGET autopas_auto4omp_omp
            PROPERTY IMPORTED_LOCATION "${AUTOPAS_AUTO4OMP_LIB_DIR}/libomp.so"
    )

    ## libomptarget.so
    set(
            AUTOPAS_TARGET_omptarget ON
            CACHE BOOL "${AUTOPAS_TARGET_omptarget_DOC}" FORCE
    )
    add_library(autopas_auto4omp_omptarget SHARED IMPORTED)
    set_property(
            TARGET autopas_auto4omp_omptarget
            PROPERTY IMPORTED_LOCATION "${AUTOPAS_AUTO4OMP_LIB_DIR}/libomptarget.so"
    )

    ## libomptarget.rtl.cuda.so
    if (TARGET omptarget.rtl.cuda)
        set(
                AUTOPAS_TARGET_omptarget.rtl.cuda ON
                CACHE BOOL "${AUTOPAS_TARGET_omptarget.rtl.cuda_DOC}" FORCE
        )
        add_library(autopas_auto4omp_omptarget.rtl.cuda SHARED IMPORTED)
        set_property(
                TARGET autopas_auto4omp_omptarget.rtl.cuda
                PROPERTY IMPORTED_LOCATION "${AUTOPAS_AUTO4OMP_LIB_DIR}/libomptarget.rtl.cuda.so"
        )
    endif ()

    ## libomptarget.rtl.x86_64.so
    if (TARGET omptarget.rtl.x86_64)
        set(
                AUTOPAS_TARGET_omptarget.rtl.x86_64 ON
                CACHE BOOL "${AUTOPAS_TARGET_mptarget.rtl.x86_64_DOC}" FORCE
        )
        add_library(autopas_auto4omp_omptarget.rtl.x86_64 SHARED IMPORTED)
        set_property(
                TARGET autopas_auto4omp_omptarget.rtl.x86_64
                PROPERTY IMPORTED_LOCATION "${AUTOPAS_AUTO4OMP_LIB_DIR}/libomptarget.rtl.x86_64.so"
        )
    endif ()

    ## libomptarget-nvptx.a
    if (TARGET omptarget-nvptx)
        set(
                AUTOPAS_TARGET_omptarget-nvptx ON
                CACHE BOOL "${AUTOPAS_TARGET_omptarget-nvptx_DOC}" FORCE
        )
        add_library(autopas_auto4omp_omptarget-nvptx SHARED IMPORTED)
        set_property(
                TARGET autopas_auto4omp_omptarget-nvptx
                PROPERTY IMPORTED_LOCATION "${AUTOPAS_AUTO4OMP_LIB_DIR}/libomptarget-nvptx.a"
        )
    endif ()

    # Prioritize Auto4OMP's libomp by prepending its path to the environment variable LD_LIBRARY_PATH.
    ## src/autopas/CMakeLists.txt will prepend AUTOPAS_LIBOMP_PATH to LD_LIBRARY_PATH before AutoPas builds. [**]
    set(
            AUTOPAS_LIBOMP_DIR "${auto4omp_BINARY_DIR}/runtime/src"
            CACHE INTERNAL "Auto4OMP's libomp.so path."
    )
    set(
            AUTOPAS_LIBOMPTARGET_DIR "${auto4omp_BINARY_DIR}/libomptarget"
            CACHE INTERNAL "Auto4OMP's libomp.so path."
    )
    set(
            AUTOPAS_LIBOMP_PATH "${AUTOPAS_LIBOMP_DIR}/libomp.so"
            CACHE INTERNAL "Auto4OMP's libomp.so path."
    )
    set(
            AUTOPAS_LIBOMPTARGET_PATH "${AUTOPAS_LIBOMPTARGET_DIR}/libomptarget.so"
            CACHE INTERNAL "Auto4OMP's libomp.so path."
    )

endif ()

# Macro to prioritize Auto4OMP's installed libraries. Prepends the install dir to LD_LIBRARY_PATH for a given target.
function(auto4omp_prioritize target)
    # TODO: Try [*] instead.
    if (DEFINED ENV{LD_LIBRARY_PATH} AND NOT $ENV{LD_LIBRARY_PATH} STREQUAL "")
        # If the variable's not empty, prepend Auto4OMP's install directory.
        add_custom_command(
                TARGET ${target}
                PRE_LINK COMMAND export LD_LIBRARY_PATH="${AUTOPAS_AUTO4OMP_LIB_DIR}:$ENV{LD_LIBRARY_PATH}"
        )
    else ()
        # Else, set it to Auto4OMP's install directory.
        add_custom_command(
                TARGET ${target}
                PRE_LINK COMMAND export LD_LIBRARY_PATH="${AUTOPAS_AUTO4OMP_LIB_DIR}"
        )
    endif ()

    #[[
    # [*] Prioritize Auto4OMP's libomp.so [https://www.hpc.dtu.dk/?page_id=1180] TODO: options unused by Clang. Why?
    target_compile_options(
            ${target} BEFORE PRIVATE
            -I"${AUTOPAS_LIBOMP_DIR}" -I"${AUTOPAS_LIBOMPTARGET_DIR}"
            -Wl,-rpath="${AUTOPAS_LIBOMPTARGET_DIR}",-rpath="${AUTOPAS_LIBOMP_DIR}"
            -L"${AUTOPAS_LIBOMP_DIR}" -lomp -L"${AUTOPAS_LIBOMPTARGET_DIR}" -lomptarget
    )
    #]]
endfunction()

# Function to simulate find_package(OpenMP), finds Auto4OMP instead of standard OpenMP.
function(auto4omp_FindOpenMP)
    if (AUTOPAS_AUTO4OMP)
        # If Auto4OMP is used, don't use for system's OMP. Instead, manually point OMP's variables to Auto4OMP.

        # Set the booleans.
        set(
                OpenMP_FOUND TRUE
                CACHE INTERNAL "Variable indicating that OpenMP flags for all requested languages have been found."
        )
        set(
                OpenMP_C_FOUND TRUE
                CACHE INTERNAL "Variable indicating if OpenMP support for C was detected.")
        set(
                OpenMP_CXX_FOUND TRUE
                CACHE INTERNAL "Variable indicating if OpenMP support for C++ was detected."
        )

        # Set the flags.
        set(
                OpenMP_C_FLAGS "-fopenmp=libomp -pthread"
                CACHE INTERNAL "OpenMP compiler flags for C, separated by spaces."
        )
        set(
                OpenMP_CXX_FLAGS "-fopenmp=libomp -pthread"
                CACHE INTERNAL "OpenMP compiler flags for C++, separated by spaces."
        )

        # Get the needed library paths.

        ## libomp.so
        get_target_property(auto4omp_libomp_IMPORTED_LOCATION omp IMPORTED_LOCATION)

        ## libomptarget.so
        get_target_property(auto4omp_libomptarget_IMPORTED_LOCATION omptarget IMPORTED_LOCATION)

        ## libpthread.a
        find_package(Threads)
        if (NOT AUTOPAS_pthread_LIBRARY AND Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)
            set(
                    AUTOPAS_pthread_LIBRARY "${CMAKE_THREAD_LIBS_INIT}"
                    CACHE FILEPATH "${AUTOPAS_pthread_LIBRARY_DOC}" FORCE
            )
        endif ()
        if (NOT AUTOPAS_pthread_LIBRARY AND (NOT Threads_FOUND OR NOT CMAKE_USE_PTHREADS_INIT GREATER_EQUAL 0))
            execute_process(
                    COMMAND find /usr/lib -name libpthread.a
                    COMMAND grep -m1 "" # Keeps first path only, inspired from [9].
                    COMMAND tr -d '\n' # Truncates newlines.
                    OUTPUT_VARIABLE AUTOPAS_pthread_LIBRARY_FIND_CMD
            )
            if (NOT ${AUTOPAS_pthread_LIBRARY_FIND_CMD} STREQUAL "")
                set(
                        AUTOPAS_pthread_LIBRARY "${AUTOPAS_pthread_LIBRARY_FIND_CMD}"
                        CACHE FILEPATH "${AUTOPAS_pthread_LIBRARY_DOC}" FORCE
                )
            else ()
                message(
                        WARNING
                        "No pthread found, Auto4OMP may warn. AutoPas attempted to look for pthread.a with FindThreads "
                        "and the command \"find /usr/lib -name FileCheck\".\n"
                        "To fix, specify the path with -DAUTOPAS_pthread_LIBRARY=\"path/to/pthread.a\""
                )
            endif ()
        endif ()

        # Set the library paths.
        set(
                OpenMP_omp_LIBRARY "${auto4omp_libomp_IMPORTED_LOCATION}"
                CACHE INTERNAL "Location of libomp.so needed for OpenMP support."
        )
        set(
                OpenMP_omptarget_LIBRARY "${auto4omp_libomptarget_IMPORTED_LOCATION}"
                CACHE INTERNAL "Location of libomptarget.so needed for OpenMP support."
        )
        set(
                OpenMP_pthread_LIBRARY "${AUTOPAS_pthread_LIBRARY}"
                CACHE INTERNAL "Location of libpthread.a needed for OpenMP support."
        )

        ## Join the paths into a list and remove empty items.
        string(
                JOIN ";" AUTOPAS_OpenMP_LIBRARIES
                "${OpenMP_omp_LIBRARY}" "${OpenMP_omptarget_LIBRARY}" "${OpenMP_pthread_LIBRARY}"
        )
        list(REMOVE_ITEM AUTOPAS_OpenMP_LIBRARIES "")

        set(
                OpenMP_C_LIBRARIES "${AUTOPAS_OpenMP_LIBRARIES}"
                CACHE INTERNAL "A list of libraries needed to link with OpenMP code in C."
        )
        set(
                OpenMP_CXX_LIBRARIES "${AUTOPAS_OpenMP_LIBRARIES}"
                CACHE INTERNAL "A list of libraries needed to link with OpenMP code in C++."
        )

        # Figure out OMP's spec date based on the used version. The list is from CMake's FindOpenMP.cmake.
        if (${AUTOPAS_OMP_VERSION} EQUAL 52)
            set(OpenMP_C_SPEC_DATE_LOCAL "202111")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 51)
            set(OpenMP_C_SPEC_DATE_LOCAL "202011")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 50)
            set(OpenMP_C_SPEC_DATE_LOCAL "201811")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 45)
            set(OpenMP_C_SPEC_DATE_LOCAL "201511")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 40)
            set(OpenMP_C_SPEC_DATE_LOCAL "201307")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 31)
            set(OpenMP_C_SPEC_DATE_LOCAL "201107")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 30)
            set(OpenMP_C_SPEC_DATE_LOCAL "200805")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 25)
            set(OpenMP_C_SPEC_DATE_LOCAL "200505")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 20)
            set(OpenMP_C_SPEC_DATE_LOCAL "200203")
        elseif (${AUTOPAS_OMP_VERSION} EQUAL 10)
            set(OpenMP_C_SPEC_DATE_LOCAL "199810")
        endif ()
        set(
                OpenMP_C_SPEC_DATE ${OpenMP_C_SPEC_DATE_LOCAL}
                CACHE INTERNAL "Date of the OpenMP specification implemented by the C compiler."
        )
        set(
                OpenMP_CXX_SPEC_DATE ${OpenMP_C_SPEC_DATE_LOCAL}
                CACHE INTERNAL "Date of the OpenMP specification implemented by the C compiler."
        )

        # Get OMP's version major and minor.
        string(SUBSTRING ${AUTOPAS_OMP_VERSION} 0 1 AUTOPAS_OMP_VERSION_MAJOR)
        string(SUBSTRING ${AUTOPAS_OMP_VERSION} 1 1 AUTOPAS_OMP_VERSION_MINOR)

        # Set the versions.
        set(
                OpenMP_C_VERSION "${AUTOPAS_OMP_VERSION_MAJOR}.${AUTOPAS_OMP_VERSION_MINOR}"
                CACHE INTERNAL "OpenMP version implemented by the C compiler."
        )
        set(
                OpenMP_CXX_VERSION "${AUTOPAS_OMP_VERSION_MAJOR}.${AUTOPAS_OMP_VERSION_MINOR}"
                CACHE INTERNAL "OpenMP version implemented by the C++ compiler."
        )
        set(
                OpenMP_C_VERSION_MAJOR ${AUTOPAS_OMP_VERSION_MAJOR}
                CACHE INTERNAL "Major version of OpenMP implemented by the C compiler."
        )
        set(
                OpenMP_CXX_VERSION_MAJOR ${AUTOPAS_OMP_VERSION_MAJOR}
                CACHE INTERNAL "Major version of OpenMP implemented by the C++ compiler."
        )
        set(
                OpenMP_C_VERSION_MINOR ${AUTOPAS_OMP_VERSION_MINOR}
                CACHE INTERNAL "Major version of OpenMP implemented by the C compiler."
        )
        set(
                OpenMP_CXX_VERSION_MINOR ${AUTOPAS_OMP_VERSION_MINOR}
                CACHE INTERNAL "Major version of OpenMP implemented by the C++ compiler."
        )
    endif ()
endfunction()

#[=====================================================================================================================[
Sources:

    [1]   https://releases.llvm.org/4.0.1/tools/clang/LanguageExtensions.html
    [2]   https://stackoverflow.com/a/41380220
    [3]   https://stackoverflow.com/a/50306091
    [4]   https://clang.llvm.org/docs/LanguageExtensions.html#has-builtin
    [5]   https://gcc.gnu.org/onlinedocs/cpp/_005f_005fhas_005fbuiltin.html
    [6]   https://stackoverflow.com/a/35693504
    [7]   https://developer.nvidia.com/cuda-gpus
    [8]   https://packages.ubuntu.com/
    [9]   https://unix.stackexchange.com/a/294493
    [10]  https://stackoverflow.com/a/71817684
    [11]  https://cmake.org/cmake/help/latest/guide/importing-exporting/index.html

    [*]   Don't use LB4OMP v0.1 (the Auto4OMP release). Clone the newer master branch instead.
          v0.1 initializes atomics as follows: atomic<int> i = 0;
          Although this is handled since C++17 [10], CMake seems to set an older C++ standard somewhere.
          Building thus fails with modern Clang due to the deleted constructor atomic(const atomic&).
          Auto4OMP's master branch now correctly uses atomic<int> i(0), compatible with all C++ standards.
          A custom fix in kmp_settings.cpp was also added by AutoPas to fix errors with older OMP versions.
          A number of LB4OMP extension function definitions are skipped for older OMP versions (line 4875),
          but used for all OMP versions (line 5252). This results in undefined errors when specifying older versions.
          To fix, move the #endif from line 5250 to 5268, so that the functions are only used for OMP 5 and above.

    [**]  Make install requires access permissions to /usr/local. It installs:
          omp.h, omp-tools.h, ompt.h at /usr/local/include; libomp.so, libomptarget.so, libomptarget.rtl.cuda.so,
          libomptarget.rtl.x86_64.so, libomptarget-nvptx.a at /usr/local/lib.

          Instead, use the local files built with make:
          omp.h, omp-tools.h at auto4omp-build/runtime/src; auto4omp-build/runtime/src/libomp.so;
          and the libomptargets at auto4omp-build/libomptarget. For now, no need to include ompt.h. It is a legacy
          renamed copy of omp-tools.h, included for compatibility purposes. However, AutoPas does not use it.

          To make sure Auto4OMP's libraries are prioritized, run the terminal command "ldd <AutoPas-executable>".
          E.g., "ldd md-flexible". libomp.so in auto4omp-build should be linked.

    [***] Some code was inspired from autopas_spdlog.cmake, autopas_eigen.cmake and other AutoPas CMake files.
#]=====================================================================================================================]

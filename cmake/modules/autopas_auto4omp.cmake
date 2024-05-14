#[=====================================================================================================================[
File: autopas_auto4omp.cmake
Author: MehdiHachicha
Date: 15.04.2024

This CMake module loads Auto4OMP into the AutoPas project.
It attempts to automatically handle its required CMake arguments, and warns if one must be manually entered.
Ideally, "cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++" should work without issues.
If cuda's installed, a gcc compatible with NVCC should also be installed. E.g., gcc-12.
#]=====================================================================================================================]

# Version settings for convenience. If Auto4OMP supports new packages, update here.
set(
        AUTOPAS_NVCC_MAX_GCC 12
        CACHE STRING "Newest gcc version supported by NVCC."
)
set(
        AUTOPAS_OMP_VERSION 45
        CACHE STRING "Newest OpenMP version supported by AutoPas."
)

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
        "To find it, run the following command: find /usr -name lit.py. "
        "By default, Auto4OMP looks for FileCheck with the PATH environment variable."
)

set(
        AUTOPAS_NVCC_GNUC_PATH_DOC
        "Path to gcc ${AUTOPAS_NVCC_MAX_GCC} or less, required by NVCC, in case system's gcc is newer."
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

# AutoPas CMake variables for Auto4OMP:
option(AUTOPAS_AUTO4OMP ${AUTOPAS_AUTO4OMP_DOC} ON)
option(AUTOPAS_AUTO4OMP_GIT ${AUTOPAS_AUTO4OMP_GIT_DOC} OFF)
set(AUTOPAS_AUTO4OMP_INSTALL_DIR "" CACHE PATH ${AUTOPAS_AUTO4OMP_INSTALL_DIR_DOC})
set(AUTOPAS_NVCC_GNUC_PATH OFF CACHE FILEPATH ${AUTOPAS_NVCC_GNUC_PATH_DOC})
set(AUTOPAS_LLVM_LIT_EXECUTABLE OFF CACHE FILEPATH ${AUTOPAS_LLVM_LIT_EXECUTABLE_DOC})
set(AUTOPAS_FILECHECK_EXECUTABLE OFF CACHE FILEPATH ${AUTOPAS_FILECHECK_EXECUTABLE_DOC})

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
        option(LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER ${LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER_DOC} ON)
    else ()
        message(STATUS "No cycle counter found. Will look for time-stamp counter next.")
    endif ()

    ### Time-stamp counter:
    #### __rdtscp() reads from the time-stamp counter, supported by x86 and x86_64. Inspired from [6]
    execute_process(COMMAND lscpu COMMAND grep rdtsc OUTPUT_VARIABLE LSCPU_GREP_RDTSC_OUTPUT)
    if (NOT ${LSCPU_GREP_RDTSC_OUTPUT} STREQUAL "")
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
    # CMake's FindCUDA is deprecated, but necessary (I think) as Auto4OMP looks for CUDA in the system.
    # For now, use old policy. TODO: for future-proofing, is there a way to use the new policy?
    cmake_policy(PUSH) # Use PUSH and POP instead of the new block() to preserve compatibility with older CMakes.
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
                    CACHE INTERNAL ${LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES_DOC}
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
                    set(AUTOPAS_NVCC_GNUC_PATH ${GCC_PATH} CACHE FILEPATH ${AUTOPAS_NVCC_GNUC_PATH_DOC} FORCE)
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
                set(AUTOPAS_NVCC_GNUC_PATH ${GCC_LESS} CACHE FILEPATH ${AUTOPAS_NVCC_GNUC_PATH_DOC} FORCE)
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
                    LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER ${AUTOPAS_NVCC_GNUC_PATH}
                    CACHE INTERNAL ${AUTOPAS_NVCC_GNUC_PATH_DOC}
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
                COMMAND grep -m1 "" #Keeps first path only, inspired from [9].
                COMMAND tr -d '\n' # Truncates newlines.
                OUTPUT_VARIABLE AUTOPAS_LLVM_LIT_EXECUTABLE
        )
    endif ()

    ### If a path to LLVM-lit executable was specified, pass it to Auto4OMP.
    if (AUTOPAS_LLVM_LIT_EXECUTABLE)
        message(DEBUG "LLVM-lit executable at ${AUTOPAS_LLVM_LIT_EXECUTABLE}")
        set(
                OPENMP_LLVM_LIT_EXECUTABLE
                ${AUTOPAS_LLVM_LIT_EXECUTABLE}
                CACHE INTERNAL ${AUTOPAS_LLVM_LIT_EXECUTABLE_DOC}
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
                OPENMP_FILECHECK_EXECUTABLE ${AUTOPAS_FILECHECK_EXECUTABLE}
                CACHE INTERNAL ${AUTOPAS_FILECHECK_EXECUTABLE_DOC}
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

    # Mark as standalone build to fix CMake's missing declaration errors.
    set(OPENMP_STANDALONE_BUILD ON)

    # Enable the FetchContent CMake module.
    include(FetchContent)

    # If OMP 4.5 is to be used, use the custom Auto4OMP zip to fix resulting build errors.
    if (${AUTOPAS_OMP_VERSION} LESS 50)
        # Build the pre-packaged Auto4OMP including the fix for build errors when specifying an old OMP version.
        FetchContent_Declare(
                auto4omp
                URL ${AUTOPAS_SOURCE_DIR}/libs/LB4OMP-custom.zip # [*]
                URL_HASH MD5=98a1eede80be95a2e93a3ccd2139ecac # Calculated with the md5sum command.
        )
    elseif (AUTOPAS_AUTO4OMP_GIT)
        # Download Auto4OMP from GIT and make the CMake targets available. TODO: untested.
        FetchContent_Declare(
                auto4omp
                GIT_REPOSITORY https://github.com/unibas-dmi-hpc/LB4OMP # [*]
        )
    else ()
        # Build the pre-packaged Auto4OMP and make the CMake targets available. [***]
        FetchContent_Declare(
                auto4omp
                URL ${AUTOPAS_SOURCE_DIR}/libs/LB4OMP-master.zip # [*]
                URL_HASH MD5=b8090fa162eb2232870916271f36650f # Calculated with the md5sum command.
        )
    endif ()

    # Use the following old policies to suppress Auto4OMP's CMake warnings.
    # cmake_policy() doesn't work, likely because Auto4OMP isn't in a subdirectory of this module's directory.
    set(CMAKE_POLICY_DEFAULT_CMP0146 OLD) # The FindCUDA module.
    set(CMAKE_POLICY_DEFAULT_CMP0148 OLD) # The FindPythonInterp and FindPythonLibs modules.

    # Integrate Auto4OMP into the project.
    FetchContent_MakeAvailable(auto4omp)

    # Reset Policies.
    set(CMAKE_POLICY_DEFAULT_CMP0146 NEW) # The FindCUDA module.
    set(CMAKE_POLICY_DEFAULT_CMP0148 NEW) # The FindPythonInterp and FindPythonLibs modules.

    # Add a target for Auto4OMP. It builds Auto4OMP by running make in its build directory.
    add_custom_target(
            auto4omp ALL
            WORKING_DIRECTORY ${auto4omp_BINARY_DIR}
            COMMAND ${CMAKE_MAKE_PROGRAM}
            COMMENT "Building Auto4OMP"
    )

    # Specify Auto4OMP's install directory.
    if (NOT AUTOPAS_AUTO4OMP_INSTALL_DIR)
        set(
                AUTOPAS_AUTO4OMP_INSTALL_DIR "${auto4omp_BINARY_DIR}/install"
                CACHE STRING ${AUTOPAS_AUTO4OMP_INSTALL_DIR_DOC} FORCE
        )
    endif ()

    # Configure cmake install to use the local install directory.
    install(DIRECTORY "${auto4omp_BINARY_DIR}/runtime/src" DESTINATION "${AUTOPAS_AUTO4OMP_INSTALL_DIR}")

    # Run cmake install post-build. Use a local directory instead of /usr, so no sudo privilege is required.
    add_custom_command(
            TARGET auto4omp POST_BUILD
            COMMAND ${CMAKE_COMMAND} --install "${auto4omp_BINARY_DIR}" --prefix "${AUTOPAS_AUTO4OMP_INSTALL_DIR}"
    )

    # Linker flags to prioritize Auto4OMP's built libraries.
    set(AUTOPAS_AUTO4OMP_LINKER_FLAGS "-L\"${auto4omp_BINARY_DIR}/install/lib\" -lomp -lomptarget")

    # Add targets for Auto4OMP's built libraries. [**]
    add_library(auto4omp_libomp SHARED IMPORTED GLOBAL)
    set_target_properties(
            auto4omp_libomp
            PROPERTIES
            IMPORTED_LOCATION "${auto4omp_BINARY_DIR}/install/lib/libomp.so"
            IMPORTED_NO_SONAME TRUE
    )

    add_library(auto4omp_libomptarget SHARED IMPORTED GLOBAL)
    set_target_properties(
            auto4omp_libomptarget
            PROPERTIES
            IMPORTED_LOCATION "${auto4omp_BINARY_DIR}/install/lib/libomptarget.so"
            IMPORTED_NO_SONAME TRUE
    )

    if (EXISTS "${auto4omp_BINARY_DIR}/libomptarget/libomptarget.rtl.cuda.so")
        add_library(auto4omp_libomptarget_rtl_cuda SHARED IMPORTED GLOBAL)
        set_target_properties(
                auto4omp_libomptarget_rtl_cuda
                PROPERTIES
                IMPORTED_LOCATION "${auto4omp_BINARY_DIR}/libomptarget/libomptarget.rtl.cuda.so"
                IMPORTED_NO_SONAME TRUE
        )
        string(APPEND AUTOPAS_AUTO4OMP_LINKER_FLAGS " -lomptarget.rtl.cuda")
    endif ()

    if (EXISTS "${auto4omp_BINARY_DIR}/libomptarget/libomptarget-nvptx.a")
        add_library(auto4omp_libomptarget_nvptx STATIC IMPORTED GLOBAL)
        set_target_properties(
                auto4omp_libomptarget_nvptx
                PROPERTIES
                IMPORTED_LOCATION "${auto4omp_BINARY_DIR}/libomptarget/libomptarget-nvptx.a"
                IMPORTED_NO_SONAME TRUE
        )
        string(APPEND AUTOPAS_AUTO4OMP_LINKER_FLAGS " -lomptarget-nvptx")
    endif ()

    if (EXISTS "${auto4omp_BINARY_DIR}/libomptarget/libomptarget.rtl.x86_64.so")
        add_library(auto4omp_libomptarget_rtl_x86_64 SHARED IMPORTED GLOBAL)
        set_target_properties(
                auto4omp_libomptarget_rtl_x86_64
                PROPERTIES
                IMPORTED_LOCATION "${auto4omp_BINARY_DIR}/libomptarget/libomptarget.rtl.x86_64.so"
                IMPORTED_NO_SONAME TRUE
        )
        string(APPEND AUTOPAS_AUTO4OMP_LINKER_FLAGS " -lomptarget.rtl.x86_64")
    endif ()

    # Set the linker flags.
    set(CMAKE_C_FLAGS "${AUTOPAS_AUTO4OMP_LINKER_FLAGS} ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${AUTOPAS_AUTO4OMP_LINKER_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${AUTOPAS_AUTO4OMP_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${AUTOPAS_AUTO4OMP_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")

    # TODO: prioritize Auto4OMP's libomp and other libs. Did linking the built libraries achieve this? [**]

endif ()

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

    [***] Some code was inspired from autopas_spdlog.cmake, autopas_eigen.cmake and other AutoPas CMake files.
#]=====================================================================================================================]

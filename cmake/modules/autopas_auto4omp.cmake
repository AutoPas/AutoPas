# Declare the Auto4OMP description.
set(AUTOPAS_AUTO4OMP_DOC "Activates Auto4OMP to automatically select OpenMP scheduling algorithms from the LB4OMP portfolio during runtime.")

# Declare the Auto4OMP CMake option.
option(AUTOPAS_AUTO4OMP ${AUTOPAS_AUTO4OMP_DOC} ON)

# Declare the CUDA CMake option. OFF by default, as it complicates the build process and requires gcc 12 or less.
option(AUTOPAS_CUDA "Optional for CUDA capable Nvidia GPUs, used by Auto4OMP." OFF)

# Declare a CMake option to specify the gcc 12 or less path for CUDA's NVCC, in case the default gcc is newer.
option(AUTOPAS_NVCC_GNUC_PATH "Specifies a path to gcc 12 or less for NVCC (used by Auto4OMP's CUDA)." "")

# If Auto4OMP disabled, warn.
if (NOT AUTOPAS_AUTO4OMP)
    message(STATUS "Auto4OMP disabled.")

# If OpenMP is disabled, warn.
elseif (NOT AUTOPAS_OPENMP)
    message(STATUS "OpenMP must be enabled to use Auto4. Auto4OMP disabled.")

# If both OpenMP and Auto4OMP enabled, build Auto4OMP (which includes Auto4OMP).
else ()
    message(STATUS "Auto4OMP enabled.")

    # Set OpenMP's standalone build option to avoid CMake missing declaration errors from Auto4OMP.
    set(OPENMP_STANDALONE_BUILD ON)

    # Enable the FetchContent CMake module.
    include(FetchContent)

    ## Optional: download Auto4OMP from GIT.
    #FetchContent_Declare(
    #        auto4omp
    #        GIT_REPOSITORY https://github.com/unibas-dmi-hpc/LB4OMP
    #        GIT_TAG v0.1 # The Auto4OMP release.
    #)

    # Or: build the pre-packaged Auto4OMP and make the CMake targets available.
    FetchContent_Declare(
            auto4omp
            URL ${AUTOPAS_SOURCE_DIR}/libs/LB4OMP-0.1.zip
            URL_HASH MD5=7bfcf0f896ed99945515a33cd5734cf0 # Calculated with the md5sum command.
    )

    #TODO: Test Auto4OMP with gcc. Seems modern gcc produces errors. What's the newest supported gcc? Same with CUDA off? Use clang for now.

    # Auto4OMP's CMake options:

    # High resolution counters:
    # TODO: Find which compilers support which counters, do the checks in the if conditions. For now, use the read-cycle-counter since it's supported by modern clang.
    if (TRUE)
        option(LIBOMP_HAVE___BUILTIN_READCYCLECOUNTER "Indicates whether the compiler supports __builtin_readcyclecounter()." ON)
    endif()
    if (FALSE)
        option(LIBOMP_HAVE___RDTSC "Indicates whether the compiler supports __rdtscp()." ON)
    endif()

    # CUDA:
    if (AUTOPAS_CUDA)
        # Modern CUDA dropped support for compute_35, which is default in Auto4OMP. Use the newer compute_53 instead.
        option(LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES "Nvidia GPU's compute capability." 53)

        # CUDA requires gcc 12 or less. If installed, pass it to NVCC to pass its compiler checks. According to LB4OMP's docs, it won't be used to produce binaries.
        if (DEFINED ENV{__GNUC__}) # If gcc found:
            if ($ENV{__GNUC__} > 12) # If system's gcc is incompatible with NVCC:
                if (NOT AUTOPAS_NVCC_GNUC_PATH STREQUAL "") # If alternate gcc specified, pass it to Auto4OMP.
                    option(LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER "Path to gcc 12 or less, required by NVCC." AUTOPAS_NVCC_GNUC_PATH)
                else ()
                    message(STATUS "CUDA enabled, but system's gcc is newer than 12. NVCC will likely complain. To fix, disable CUDA or install e.g. gcc-12 and and pass its path (e.g. \"/usr/bin/gcc-12\") with the CMake option LIBOMPTARGET_NVPTX_ALTERNATE_HOST_COMPILER.")
                endif ()
            endif () # Do nothing if system's gcc is compatible.
        else ()
            message(STATUS "CUDA enabled, but no gcc found. NVCC will likely complain. To fix, install gcc 12 or less.")
        endif ()
    else ()
        # TODO: If CUDA off, make sure Auto4OMP doesn't look for it. If it finds it, it'll enable it by default, then produce errors because the above options weren't set.
    endif ()

    # Integrate Auto4OMP into the project.
    FetchContent_MakeAvailable(auto4omp)
    include_directories(SYSTEM "${auto4omp_SOURCE_DIR}/runtime/src")
    link_directories(SYSTEM "${auto4omp_SOURCE_DIR}/runtime/src")

endif ()
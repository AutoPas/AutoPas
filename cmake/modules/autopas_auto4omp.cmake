# Declare the Auto4OMP description.
set(AUTOPAS_AUTO4OMP_DOC "Activates Auto4OMP to automatically select OpenMP scheduling algorithms from the LB4OMP portfolio during runtime.")

# Declare the Auto4OMP CMake option.
option(AUTOPAS_AUTO4OMP ${AUTOPAS_OPENMP_DOC} ON)

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

    # Integrate into the project.
    FetchContent_MakeAvailable(auto4omp)
    include_directories(SYSTEM "${auto4omp_SOURCE_DIR}/runtime/src")
    link_directories(SYSTEM "${auto4omp_SOURCE_DIR}/runtime/src")

endif ()
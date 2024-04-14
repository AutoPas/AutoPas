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

# If both OpenMP and Auto4OMP enabled, build LB4OMP (which includes Auto4OMP).
else ()
    message(STATUS "Auto4OMP enabled.")

    # Enable the FetchContent CMake module.
    include(FetchContent)

    ## Optional: download from GIT. If used, comment out the declaration for the pre-packaged library.
    #FetchContent_Declare(
    #        LB4OMP
    #        GIT_REPOSITORY https://github.com/unibas-dmi-hpc/LB4OMP
    #        GIT_TAG v0.1 # The Auto4OMP release.
    #)

    # Build the pre-packaged LB4OMP and make the CMake targets available.
    FetchContent_Declare(
            LB4OMP
            URL ${AUTOPAS_SOURCE_DIR}/libs/LB4OMP-0.1.zip
            URL_HASH MD5=7bfcf0f896ed99945515a33cd5734cf0 # Calculated with the md5sum command.
    )

    # Fetch potential past population records.
    FetchContent_GetProperties(LB4OMP)

    # If not populated, populate and add the LB4OMP library.
    if (NOT LB4OMP_POPULATED)
        FetchContent_Populate(LB4OMP)
        add_library(LB4OMP STATIC IMPORTED GLOBAL)
    endif ()
endif ()
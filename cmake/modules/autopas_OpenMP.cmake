set(AUTOPAS_OPENMP_DOC "Activates OpenMP shared memory parallelization. (requires OpenMP 4.5)")
option(AUTOPAS_OPENMP ${AUTOPAS_OPENMP_DOC} ON)

# Include Auto4OMP's CMake module to access its variables. Must be included after AUTOPAS_OPENMP is set.
include(autopas_auto4omp)

if (AUTOPAS_OPENMP)
    message(STATUS "OpenMP enabled.")
    find_package(OpenMP REQUIRED)
    if (
        NOT AUTOPAS_AUTO4OMP
        # OpenMP version 4.5 was specified in 11.2015
        AND (OpenMP_CXX_SPEC_DATE LESS 201511)
        AND
            NOT
            (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
             AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6)
             )
    )
        # Disable AUTOPAS_OPENMP if no sufficient OpenMP version is available.
        set(AUTOPAS_OPENMP OFF CACHE BOOL ${AUTOPAS_OPENMP_DOC} FORCE )
        message(
            FATAL_ERROR
                "OpenMP version not supported (specification date: ${OpenMP_CXX_SPEC_DATE}). Required version: 4.5+"
        )
    endif ()
else ()
    message(STATUS "OpenMP disabled.")
endif ()

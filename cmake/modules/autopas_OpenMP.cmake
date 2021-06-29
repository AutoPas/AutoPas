option(AUTOPAS_OPENMP "Activates OpenMP shared memory parallelization." OFF)

if (AUTOPAS_OPENMP)
    message(STATUS "OpenMP enabled.")
    find_package(OpenMP REQUIRED)
    # OpenMP version 4.5 was specified in 11.2015
    if (
        (OpenMP_CXX_SPEC_DATE LESS 201511)
        AND
            NOT
            (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
             AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6)
             )
    )
        message(
            FATAL_ERROR
                "OpenMP version not supported (specification date: ${OpenMP_CXX_SPEC_DATE}). Required version: 4.5+"
        )
    endif ()
else ()
    message(STATUS "OpenMP disabled.")
endif ()

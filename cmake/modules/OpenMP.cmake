option(AUTOPAS_OPENMP "Activates OpenMP shared memory parallelization." OFF)

if (ARCHER)
    message(STATUS "archer detected, OpenMP enabled by default, so skipping OpenMP package search")
    set(AUTOPAS_OPENMP ON)
elseif (AUTOPAS_OPENMP)
    message(STATUS "OpenMP enabled.")
    find_package(OpenMP REQUIRED)

    # OpenMP version 4.5 was specified in 11.2015
    if (OpenMP_CXX_SPEC_DATE LESS 201511)
        message(
            FATAL_ERROR
                "OpenMP version not supported (specification date: ${OpenMP_CXX_SPEC_DATE}). Required version: 4.5+"
        )
    endif ()
else ()
    message(STATUS "OpenMP disabled.")
endif ()

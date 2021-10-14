option(AUTOPAS_INTERNODE_TUNING "Activates distributed tuning in MPI parallel simulations" OFF)

# messages first
if (AUTOPAS_INTERNODE_TUNING)
    message(STATUS "Internode tuning enabled. Adding MPI.")
else ()
    message(STATUS "Internode tuning disabled.")
endif ()

if (MD_FLEXIBLE_USE_MPI)
    message(STATUS "Adding MPI for md-flexible.")
endif ()

# actual action
if (AUTOPAS_INTERNODE_TUNING OR MD_FLEXIBLE_USE_MPI)
    find_package(MPI REQUIRED)
    # set flag that will be used to include the macro with the same name
    set(AUTOPAS_INCLUDE_MPI true)
    message(STATUS "MPI compiler found: ${MPI_CXX_COMPILER}")
endif ()

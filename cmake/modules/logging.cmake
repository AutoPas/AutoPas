# option for more verbose log messages.
option(AUTOPAS_VERBOSE_LOGGING "Print Filename and line in every log line." OFF)

if (AUTOPAS_VERBOSE_LOGGING)
    target_compile_definitions(autopas PUBLIC AUTOPAS_VERBOSE_LOG)
    message(STATUS "Verbose log messages enabled.")
endif ()

# option for colored log messages.
option(
        AUTOPAS_COLORED_LOGGING "Print colored logging messages (only applies to cout and cerr)." OFF
)
if (AUTOPAS_COLORED_LOGGING)
    target_compile_definitions(autopas PUBLIC AUTOPAS_COLORED_CONSOLE_LOGGING)
    message(STATUS "Colored log messages enabled.")
endif ()

# option for GaussianClusterLogger
option(AUTOPAS_Log_GaussianCluster "Generate a csv file about the gaussian cluster model that can be used for plotting." OFF)
if (AUTOPAS_Log_GaussianCluster)
    target_compile_definitions(autopas PUBLIC AUTOPAS_Log_GaussianCluster)
    message(STATUS "GaussianClusterLogger enabled.")
endif ()

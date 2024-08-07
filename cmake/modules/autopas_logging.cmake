# option for colored log messages.
option(AUTOPAS_COLORED_LOGGING "Print colored logging messages (only applies to cout and cerr)." OFF)
if (AUTOPAS_COLORED_LOGGING)
    target_compile_definitions(autopas PUBLIC AUTOPAS_COLORED_CONSOLE_LOGGING)
    message(STATUS "Colored log messages enabled.")
endif ()

# option for GaussianClusterLogger
option(AUTOPAS_LOG_GAUSSIANCLUSTER "Generate a csv file about the gaussian cluster model that can be used for plotting." OFF)
if (AUTOPAS_LOG_GAUSSIANCLUSTER OR AUTOPAS_LOG_ALL)
    target_compile_definitions(autopas PUBLIC AUTOPAS_LOG_GAUSSIANCLUSTER)
    message(STATUS "GaussianClusterLogger enabled.")
endif ()

# option for PredictionsLogger
option(AUTOPAS_LOG_PREDICTIONS "Generate a csv file about the predictive tuning strategy that can be used for plotting." OFF)
if (AUTOPAS_LOG_PREDICTIONS OR AUTOPAS_LOG_ALL)
    target_compile_definitions(autopas PUBLIC AUTOPAS_LOG_PREDICTIONS)
    message(STATUS "PredictionsLogger enabled.")
endif ()

# option for IterationLogger
option(AUTOPAS_LOG_ITERATIONS "Generate a csv tracking the performance of iteratePairwise calls." OFF)
if (AUTOPAS_LOG_ITERATIONS OR AUTOPAS_LOG_ALL)
    target_compile_definitions(autopas PUBLIC AUTOPAS_LOG_ITERATIONS)
    message(STATUS "IterationLogger enabled.")
endif ()

# option for FLOPLogger
option(AUTOPAS_LOG_FLOPS "Generate a csv tracking the FLOP count and hit rate." OFF)
if (AUTOPAS_LOG_FLOPS OR AUTOPAS_LOG_ALL)
    target_compile_definitions(autopas PUBLIC AUTOPAS_LOG_FLOPS)
    message(STATUS "FLOPLogger enabled.")
endif ()

# option for TuningResultLogger
option(AUTOPAS_LOG_TUNINGRESULTS "Generate a csv tracking the decisions of the auto tuner." OFF)
if (AUTOPAS_LOG_TUNINGRESULTS OR AUTOPAS_LOG_ALL)
    target_compile_definitions(autopas PUBLIC AUTOPAS_LOG_TUNINGRESULTS)
    message(STATUS "TuningResultLogger enabled.")
endif ()

# option for TuningDataLogger
option(AUTOPAS_LOG_TUNINGDATA "Generate a csv tracking data collected by the auto tuner." OFF)
if (AUTOPAS_LOG_TUNINGDATA OR AUTOPAS_LOG_ALL)
    target_compile_definitions(autopas PUBLIC AUTOPAS_LOG_TUNINGDATA)
    message(STATUS "TuningDataLogger enabled.")
endif ()

# option for LiveInfoLogger
option(AUTOPAS_LOG_LIVEINFO "Generate a csv tracking the live information of the autopas object." OFF)
if (AUTOPAS_LOG_LIVEINFO OR AUTOPAS_LOG_ALL)
    target_compile_definitions(autopas PUBLIC AUTOPAS_LOG_LIVEINFO)
    message(STATUS "LiveInfoLogger enabled.")
endif ()
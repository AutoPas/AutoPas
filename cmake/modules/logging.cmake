# option for more verbose log messages.
option(AUTOPAS_VERBOSE_LOGGING "Print Filename and line in every log line." OFF)
# option for colored log messages.
option(
    AUTOPAS_COLORED_LOGGING "Print colored logging messages (only applies to cout and cerr)." OFF
)

if (AUTOPAS_VERBOSE_LOGGING)
    target_compile_definitions(autopas PUBLIC AUTOPAS_VERBOSE_LOG)
    message(STATUS "Verbose log messages enabled.")
endif ()

if (AUTOPAS_COLORED_LOGGING)
    target_compile_definitions(autopas PUBLIC AUTOPAS_COLORED_CONSOLE_LOGGING)
    message(STATUS "Colored log messages enabled.")
endif ()

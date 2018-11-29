# option for more verbose log messages.
option(VERBOSE_LOGGING "Print Filename and line in every log line." OFF)

if (NOT VERBOSE_LOGGING)
    return()
endif ()

target_compile_definitions(autopas
        PUBLIC
        AUTOPAS_VERBOSE_LOG
        )

message(STATUS "Verbose log messages enabled.")

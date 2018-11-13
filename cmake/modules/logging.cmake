# option for more verbose log messages.
option(VERBOSE_LOGGING "Print Filename and line in every log line." OFF)
target_compile_definitions(autopas
        INTERFACE
        $<$<BOOL:${VERBOSE_LOGGING}>:AUTOPAS_VERBOSE_LOG>
)

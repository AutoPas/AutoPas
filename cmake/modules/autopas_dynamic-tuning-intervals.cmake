set(AUTOPAS_DYNAMIC_TUNING_INTERVALS_DOC "Changes from static retuning intervals to dynamic retune interval estimation.")
option(AUTOPAS_DYNAMIC_TUNING_INTERVALS ${AUTOPAS_DYNAMIC_TUNING_INTERVALS_DOC} OFF)

if (AUTOPAS_DYNAMIC_TUNING_INTERVALS)
    message(STATUS "Dynamic tuning interval estimation enabled.")
    target_compile_definitions(autopas PUBLIC AUTOPAS_DYNAMIC_TUNING_INTERVALS_ENABLED)
else ()
    message(STATUS "OpenMP disabled.")
endif ()
 
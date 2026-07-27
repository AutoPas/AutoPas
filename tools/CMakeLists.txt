add_subdirectory(generators)

option(AUTOPAS_BUILD_TOOLS "Build stand-alone tools (ruleChecker, tuningLogToSQL)." OFF)

if (NOT AUTOPAS_BUILD_TOOLS)
    return()
endif ()

add_subdirectory(ruleChecker)
add_subdirectory(tuningLogToSQL)

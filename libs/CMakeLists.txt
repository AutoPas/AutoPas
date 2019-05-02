#spdlog
find_package(spdlog CONFIG QUIET)
if (spdlog_FOUND)
    message(STATUS "spdlog found - using system version")
else()
    message(STATUS "spdlog not found - using bundled version")
    TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC spdlog-1.3.1/include)
endif ()

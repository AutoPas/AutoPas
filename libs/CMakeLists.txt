
#spdlog
find_package(spdlog CONFIG)
if (spdlog_FOUND)
    message(STATUS "spdlog found - using system version")
else()
    message(STATUS "spdlog not found - using bundled version")
    TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC spdlog-0.16.3/include)
endif ()

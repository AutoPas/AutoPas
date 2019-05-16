#spdlog
find_package(spdlog CONFIG QUIET)
if (spdlog_FOUND)
    message(STATUS "spdlog found - using system version")
else()
    message(STATUS "spdlog not found - using bundled version")
    TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC spdlog-0.16.3/include)
endif ()

#eigen3
TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC eigen-3.3.7)

find_package(Threads REQUIRED)

find_package(spdlog)

if (spdlog_FOUND)
    message(STATUS "spdlog found - using system version")
else()
    message(STATUS "spdlog not found - using bundled version")
    ADD_SUBDIRECTORY(spdlog-0.16.3)
endif ()
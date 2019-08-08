#spdlog
message(STATUS "spdlog - using bundled version")
TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC spdlog-1.3.1/include)

#eigen3
include(autopas_eigen)

#spdlog
message(STATUS "spdlog - using bundled version")
TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC spdlog-1.3.1/include)

#active-harmony
TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC harmony/include)

#eigen3
include(autopas_eigen)

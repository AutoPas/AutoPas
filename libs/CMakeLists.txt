#spdlog
message(STATUS "spdlog - using bundled version")
TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC spdlog-1.3.1/include)

#eigen3
TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC eigen-3.3.7)

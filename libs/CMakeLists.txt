#spdlog
message(STATUS "spdlog - using bundled version")
message(STATUS "frugally-deep - using bundled version with eigen, json and functionalPlus")
TARGET_INCLUDE_DIRECTORIES(autopas SYSTEM PUBLIC spdlog-1.3.1/include FunctionalPlus-master/include
        eigen-eigen-323c052e1731/include json-single/include frugally-deep-master/include)

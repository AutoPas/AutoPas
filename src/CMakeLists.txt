file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_LIBRARY(autopas STATIC
        ${MY_SRC}
        )

target_link_libraries(autopas rt) # required for Time.h
target_link_libraries(autopas spdlog) # required for Logger

TARGET_INCLUDE_DIRECTORIES(autopas INTERFACE .)

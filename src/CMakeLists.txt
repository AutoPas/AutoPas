file(GLOB_RECURSE MY_SRC
        "*.cpp"
        )

ADD_LIBRARY(autopas STATIC
        ${MY_SRC})
TARGET_INCLUDE_DIRECTORIES(autopas INTERFACE .)

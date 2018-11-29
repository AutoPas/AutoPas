file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(cluster-test
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(cluster-test
        autopas)

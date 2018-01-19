file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(md-main
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(md-main
        autopas)

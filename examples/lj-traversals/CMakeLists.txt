file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(lj-traversals
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(lj-traversals
        autopas)

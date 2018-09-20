file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(sph-traversals
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-traversals
        autopas)
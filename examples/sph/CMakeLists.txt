file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(sph-main
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-main
        autopas)

file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(sph-main-verlet
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-main-verlet
        autopas)

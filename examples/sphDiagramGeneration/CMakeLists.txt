file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(sph-diagram-generation
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(sph-diagram-generation
        autopas)

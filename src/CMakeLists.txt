SET(MY_SRC
        utils/Timer.cpp
        particles/MoleculeLJ.cpp
        )

ADD_LIBRARY(autopas STATIC
        ${MY_SRC})
TARGET_INCLUDE_DIRECTORIES(autopas INTERFACE .)

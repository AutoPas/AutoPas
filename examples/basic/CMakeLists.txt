SET(MY_SRC
        main.cpp
        )

ADD_EXECUTABLE(basic
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(basic
        autopas)

SET(MY_SRC
        md-main.cpp
        )

ADD_EXECUTABLE(md-main
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(md-main
        autopas)

file(GLOB_RECURSE MY_SRC
        "*.cpp"
        )


ADD_EXECUTABLE(basic
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(basic
        autopas)

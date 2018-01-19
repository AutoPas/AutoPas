file(GLOB_RECURSE MY_SRC
        "*.cpp"
        )


ADD_EXECUTABLE(runTests
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(runTests
        autopas
        libgtest
        )

ADD_TEST(NAME testAutopas
        COMMAND runTests)
SET(MY_SRC
        allTests.cpp
        myMoleculeTest.cpp
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
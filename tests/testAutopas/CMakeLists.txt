file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(runTests
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(runTests
        autopas
        libgtest
        libgmock
        )

include(FindGTest)


GTEST_ADD_TESTS(runTests
        ""
        ${MY_SRC}
        )

#ADD_TEST(NAME testAutopas
#        COMMAND runTests)

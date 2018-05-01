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


#these tests produce nice output for ctest and make test
GTEST_ADD_TESTS(runTests
        ""
        ${MY_SRC}
        )

#this test is needed for a proper generation of xml files of the test
ADD_TEST(NAME testAutopas
        COMMAND runTests)

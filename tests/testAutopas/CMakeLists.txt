FILE(GLOB_RECURSE MY_SRC
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

if(CMAKE_VERSION VERSION_LESS "3.10")
    include(FindGTest)
    # uses regexp expressions (not very robust)
    GTEST_ADD_TESTS(runTests
        ""
        ${MY_SRC}
        )
else()
    # this cmake module was only introduced in 3.10
    include(GoogleTest)
    # more robust, queries the compiled executable
    GTEST_DISCOVER_TESTS(runTests)
endif()

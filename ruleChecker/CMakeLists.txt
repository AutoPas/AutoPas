file(
        GLOB_RECURSE
        ruleChecker_SRC
        "*.cpp"
        "*.h"
)

add_executable(ruleChecker ${ruleChecker_SRC})

target_include_directories(ruleChecker PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ruleChecker PUBLIC autopas autopasTools)

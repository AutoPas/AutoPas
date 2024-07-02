if(AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING)

    file(
            GLOB_RECURSE
            ruleChecker_SRC
            "*.cpp"
            "*.h"
    )

    add_executable(ruleChecker ${ruleChecker_SRC})

    target_include_directories(ruleChecker PUBLIC .)

    target_link_libraries(ruleChecker PUBLIC autopas autopasTools)

endif()

set(AUTOPAS_DISABLE_RULES_BASED_TUNING
        OFF
        CACHE
        BOOL "Disables rule-based tuning forcing no compilation of ANTLR or uuid. May speed up compilation time and fix compatibility problems"
        )

if(NOT AUTOPAS_DISABLE_RULES_BASED_TUNING)

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

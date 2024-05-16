set(AUTOPAS_ENABLE_RULES_BASED_TUNING
        ON
        CACHE
        BOOL "Enables rules-based tuning, which, if using the bundled version, will compile ANTLR and, if uuid is not installed already, will compile uuid which can cause errors."
        )

if(AUTOPAS_ENABLE_RULES_BASED_TUNING)

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

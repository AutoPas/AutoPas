file(
    GLOB_RECURSE
    MY_SRC
    "*.cpp"
    "*.h"
)

file(
    GLOB_RECURSE
    CU_SRC
    "*.cu"
    "*.cuh"
)
set(AUTOPAS_ENABLE_RULES_BASED_TUNING
        ON
        CACHE
        BOOL "Enables rules-based tuning, which, if using the bundled version, will compile ANTLR and, if uuid is not installed already, will compile uuid which can cause errors."
        )

# If rules-based tuning is disabled, we must remove all files which link to antlr. As this includes generated files,
# where we can't easily change the file, we disable the files here.
if (NOT AUTOPAS_ENABLE_RULES_BASED_TUNING)
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/parser_generated/autopas_generated_rule_syntax/RuleLanguageBaseVisitor.h")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/parser_generated/autopas_generated_rule_syntax/RuleLanguageBaseVisitor.cpp")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/parser_generated/autopas_generated_rule_syntax/RuleLanguageVisitor.h")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/parser_generated/autopas_generated_rule_syntax/RuleLanguageVisitor.cpp")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/parser_generated/autopas_generated_rule_syntax/RuleLanguageLexer.h")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/parser_generated/autopas_generated_rule_syntax/RuleLanguageLexer.cpp")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/parser_generated/autopas_generated_rule_syntax/RuleLanguageParser.h")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/parser_generated/autopas_generated_rule_syntax/RuleLanguageParser.cpp")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedProgramParser.h")
    list(REMOVE_ITEM MY_SRC "${AUTOPAS_SOURCE_DIR}/src/autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedProgramParser.cpp")

    message(STATUS "Here")

    message(STATUS ${MY_SRC})
endif ()


add_library(autopas STATIC ${MY_SRC})

target_link_libraries(
    autopas
    PUBLIC
        rt # required for Time.h
        ${CMAKE_THREAD_LIBS_INIT} # required for Logger and ExceptionHandler
        $<$<BOOL:${AUTOPAS_OPENMP}>:OpenMP::OpenMP_CXX>
        spdlog::spdlog
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:MPI::MPI_CXX>
        $<$<OR:$<BOOL:${AUTOPAS_INTERNODE_TUNING}>,$<BOOL:${MD_FLEXIBLE_USE_MPI}>>:MPI::MPI_CXX>
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:MPI::MPI_CXX>
        Eigen3
    PRIVATE
        # harmony are only needed privately when building AutoPas.
        harmony
)

# If rules-based tuning is not disabled, also link antlr4cpp
if (AUTOPAS_ENABLE_RULES_BASED_TUNING)
    target_link_libraries(autopas PUBLIC antlr4cpp)
endif()

# OpenMP stuff needs to be here because autopas_OpenMP.cmake needs to run before this file to create to OpenMP
# target. this can be resolved by upgrading to CMake 3.13 and enforcing CMP0079.
target_compile_definitions(
    autopas
    PUBLIC
        $<$<BOOL:${AUTOPAS_OPENMP}>:AUTOPAS_USE_OPENMP>
        $<$<NOT:$<BOOL:${AUTOPAS_OPENMP}>>:EIGEN_DONT_PARALLELIZE>
        $<$<BOOL:${AUTOPAS_INCLUDE_MPI}>:AUTOPAS_INCLUDE_MPI>
        $<$<BOOL:${AUTOPAS_INTERNODE_TUNING}>:AUTOPAS_INTERNODE_TUNING>
        $<$<BOOL:${AUTOPAS_ENABLE_RULES_BASED_TUNING}>:AUTOPAS_ENABLE_RULES_BASED_TUNING>
        _USE_MATH_DEFINES
)


target_include_directories(autopas PUBLIC ${AUTOPAS_SOURCE_DIR}/src/)
target_include_directories(autopas PUBLIC ${AUTOPAS_BINARY_DIR}/src/)

option(AUTOPAS_ENABLE_CLANG_TIDY "Add clang-tidy automatically to builds")

if (AUTOPAS_ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy" DOC "Path to clang-tidy executable")
    if (CLANG_TIDY_EXE)
        message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
        set(CLANG_TIDY_CHECKS "-*,modernize-*,-clang-analyzer-alpha.*")
        message(STATUS "cmake source dir: ${AUTOPAS_SOURCE_DIR}")
        set(
            CMAKE_CXX_CLANG_TIDY
            "${CLANG_TIDY_EXE};-checks=${CLANG_TIDY_CHECKS};-header-filter='${AUTOPAS_SOURCE_DIR}/*'"
            CACHE STRING "" FORCE
        )
    else ()
        message(AUTHOR_WARNING "clang-tidy not found!")
        set(CMAKE_CXX_CLANG_TIDY "" CACHE STRING "" FORCE) # delete it
    endif ()
else ()
    message(
        STATUS
            "clang tidy disabled, to enable specify -DAUTOPAS_ENABLE_CLANG_TIDY=\"ON\" when calling cmake"
    )
endif ()

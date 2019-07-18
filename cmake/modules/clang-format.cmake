# additional target to perform clang-format run, requires clang-format

# get all project files
set(
    INCLUDE_DIRS
    "${PROJECT_SOURCE_DIR}/src/" "${PROJECT_SOURCE_DIR}/examples/" "${PROJECT_SOURCE_DIR}/tests/"
)
file(
    GLOB_RECURSE
    ALL_SOURCE_FILES
    "*.cpp"
    "*.h"
    "*.cuh"
    "*.cu"
)
foreach (TMP_PATH ${ALL_SOURCE_FILES})
    set(_found FALSE)
    foreach (_incdir ${INCLUDE_DIRS})
        string(FIND ${TMP_PATH} ${_incdir} INCLUDE_DIR_FOUND)
        if (${INCLUDE_DIR_FOUND} EQUAL 0)
            set(_found TRUE)
        endif ()
    endforeach (_incdir)
    if (NOT ${_found})
        list(REMOVE_ITEM ALL_SOURCE_FILES ${TMP_PATH})
    endif ()
endforeach (TMP_PATH)

find_program(CLANG_FORMAT NAMES clang-format-8)

if (CLANG_FORMAT)
    message(STATUS "clang format found, added clangformat target")
    set(dummyfiles)
    foreach (_file ${ALL_SOURCE_FILES})
        string(
            REPLACE
                "."
                "_"
                file_cf
                ${_file}
        )
        string(
            REPLACE
                ".."
                "."
                file_cf
                ${file_cf}
        )
        set(file_cf ".dummy/cf/${file_cf}_cf")
        add_custom_command(
            OUTPUT ${file_cf}
            COMMAND
                ${CLANG_FORMAT}
                -style=file
                -i
                ${_file}
            DEPENDS ${_file}
        )
        list(APPEND dummyfiles ${file_cf})
    endforeach ()
    add_custom_command(OUTPUT .dummy/cf/clang_dummy COMMAND true DEPENDS ${dummyfiles})
    add_custom_target(clangformat DEPENDS .dummy/cf/clang_dummy)
else ()
    message(
        STATUS
            "clang-format-8 not found, not adding clang format target. Other Versions not supported!"
    )
endif ()

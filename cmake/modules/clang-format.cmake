# additional target to perform clang-format run, requires clang-format

# get all project files
file(GLOB_RECURSE ALL_SOURCE_FILES *.cpp *.h)


set (EXCLUDE_DIR "/libs/")
file (GLOB_RECURSE ALL_SOURCE_FILES "*.cpp" "*.h")
foreach (TMP_PATH ${ALL_SOURCE_FILES})
    string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM ALL_SOURCE_FILES ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)


add_custom_target(
        clangformat
        COMMAND /usr/bin/clang-format
        -style=Google
        -i
        ${ALL_SOURCE_FILES}
)
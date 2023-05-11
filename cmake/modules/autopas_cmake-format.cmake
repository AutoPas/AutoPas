# additional target to perform cmake-format run, requires cmake-format

# get all CMake related files
set(
    INCLUDE_DIRS
    "${PROJECT_SOURCE_DIR}/cmake/"
    "${PROJECT_SOURCE_DIR}/src/"
    "${PROJECT_SOURCE_DIR}/examples/"
    "${PROJECT_SOURCE_DIR}/tests/"
    "${PROJECT_SOURCE_DIR}/tools/"
    "${PROJECT_SOURCE_DIR}/applicationLibrary/"
)

# set ALL_CMake_FILES to cmake files in PROJECT_SOURCE_DIR, as we cannot do a recurse there.
set(ALL_CMake_FILES "${PROJECT_SOURCE_DIR}/CMakeLists.txt" "${PROJECT_SOURCE_DIR}/version.cmake")

foreach (TMP_PATH ${INCLUDE_DIRS})
    # search for files for each path in INCLUDE_DIRS and append them to ALL_CMake_FILES
    file(
        GLOB_RECURSE
        ALL_CMake_FILES_TMP
        "${TMP_PATH}/CMakeLists.txt"
        "${TMP_PATH}/*.cmake"
    )

    list(APPEND ALL_CMake_FILES ${ALL_CMake_FILES_TMP})
endforeach (TMP_PATH)

find_program(CMAKE_FORMAT NAMES cmake-format)

if (CMAKE_FORMAT)
    message(STATUS "CMake format found, added cmakeformat target")
    set(dummyfiles)
    foreach (_file ${ALL_CMake_FILES})
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
                ${CMAKE_FORMAT}
                --in-place
                --config-file
                ${PROJECT_SOURCE_DIR}/cmake/cmake-format.py
                ${_file}
            DEPENDS ${_file}
        )
        list(APPEND dummyfiles ${file_cf})
    endforeach ()
    add_custom_command(OUTPUT .dummy/cf/cmake_dummy COMMAND true DEPENDS ${dummyfiles})
    add_custom_target(cmakeformat DEPENDS .dummy/cf/cmake_dummy)
else ()
    message(STATUS "CMake format not found, not adding cmakeformat target")
endif ()

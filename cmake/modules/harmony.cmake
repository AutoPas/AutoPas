message(STATUS "harmony - using bundled version")

# Enable FetchContent CMake module
include(FetchContent)

# Extract and build harmony
FetchContent_Declare(
    harmony
    URL ${AUTOPAS_SOURCE_DIR}/libs/harmony.zip
    URL_HASH MD5=a8768c2886bdc2e44e3b6b7d4f94729c
    PREFIX
    ${CMAKE_CURRENT_BINARY_DIR}/harmony
    # since we only unpack a header lib src == include
    SOURCE_DIR
    ${CMAKE_CURRENT_BINARY_DIR}/harmony/include
)

# Check if population has already been performed by project including AutoPas as dependency
FetchContent_GetProperties(harmony)
if (NOT harmony_POPULATED)
    # Fetch the content using previously declared details
    FetchContent_Populate(harmony)

    # Check if library has already been built
    if (EXISTS "${harmony_SOURCE_DIR}/lib/libharmony.a")
        message(STATUS "harmony - using previous build")
    else ()
        # harmony is a make file project run make all and make install manually
        find_program(MAKE_EXE NAMES gmake nmake make)
        message(STATUS "harmony - building library...")
        execute_process(
            COMMAND ${MAKE_EXE} CFLAGS=-w
            WORKING_DIRECTORY ${harmony_SOURCE_DIR}
            OUTPUT_FILE ${harmony_BINARY_DIR}/build_output.log
            ERROR_FILE ${harmony_BINARY_DIR}/build_output.log
            RESULT_VARIABLE result
        )

        if (result)
            message(
                FATAL_ERROR
                    "Failed harmony build, see build log at:\n"
                    "    ${harmony_BINARY_DIR}/build_output.log"
            )
        endif ()

        execute_process(
            COMMAND ${MAKE_EXE} install DESTDIR=${harmony_PREFIX}
            WORKING_DIRECTORY ${harmony_SOURCE_DIR}
            OUTPUT_FILE ${harmony_BINARY_DIR}/install_output.log
            ERROR_FILE ${harmony_BINARY_DIR}/install_output.log
            RESULT_VARIABLE result
        )

        if (result)
            message(
                FATAL_ERROR
                    "Failed harmony install, see install log at:\n"
                    "    ${harmony_BINARY_DIR}/install_output.log"
            )
        endif ()

        message(STATUS "harmony - build complete")
    endif ()

    # Bring the populated content into the build
    add_library(
        harmony
        STATIC
        IMPORTED
        GLOBAL
    )

    target_link_libraries(harmony INTERFACE ${harmony_SOURCE_DIR}/lib/libharmony.a)

    set_target_properties(
        harmony
        PROPERTIES "IMPORTED_LOCATION" "${harmony_SOURCE_DIR}/lib/libharmony.a"
    )

    target_include_directories(harmony SYSTEM INTERFACE "${harmony_SOURCE_DIR}/include")

    # Set macro needed to set environment variable for ActiveHarmony
    target_compile_definitions(harmony INTERFACE HARMONY_HOME="HARMONY_HOME=${harmony_SOURCE_DIR}")

endif ()

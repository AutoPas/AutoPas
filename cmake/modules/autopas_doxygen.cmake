# make doc_doxygen optional if someone does not have / like doxygen
set(AUTOPAS_BUILD_TARGET_DOC_DOC "Create \"make doc_doxygen\" target (requires Doxygen)")
cmake_dependent_option(AUTOPAS_BUILD_TARGET_DOC
        ${AUTOPAS_BUILD_TARGET_DOC_DOC}
        ON
        AUTOPAS_STANDALONE_BUILD
        OFF
)

# do nothing if nothing should be done
if (NOT AUTOPAS_BUILD_TARGET_DOC)
    return()
endif ()

# Other versions will work but might throw warnings or incomplete documentation.
set(DOXYGEN_RECOMMENDED_VERSION 1.9.2)

# check if Doxygen is installed
find_package(Doxygen COMPONENTS dot OPTIONAL_COMPONENTS mscgen dia)
if (DOXYGEN_FOUND)
    if (${DOXYGEN_VERSION} VERSION_LESS ${DOXYGEN_RECOMMENDED_VERSION})
        # e.g. because of anonymous template arguments with default values in Math.h
        message(WARNING "Doxygen - Versions before ${DOXYGEN_RECOMMENDED_VERSION} might produce incomplete documentation.")
    endif ()
    # set input and output files
    set(DOXY_CONF_DIR docs)
    set(DOXYGEN_IN ${DOXY_CONF_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    set(DOXY_MAIN_PAGE README.md)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    add_custom_target(
        doc_doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )

    if(AUTOPAS_BUILD_EXAMPLES)
        set(DOXYGEN_MD-FLEXIBLE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile_md-flexible)
        # request to configure the file
        configure_file(examples/md-flexible/Doxyfile.in ${DOXYGEN_MD-FLEXIBLE_OUT} @ONLY)
        add_custom_target(
                doc_doxygen_md-flexible
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_MD-FLEXIBLE_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating documentation for md-flexible with Doxygen"
                VERBATIM
        )
    endif()

    message(STATUS "Doxygen - Configured")
else ()
    message(
        WARNING
            "Doxygen needs to be installed to generate the doxygen documentation, you might also have to install dot (graphviz)"
    )
    set(AUTOPAS_BUILD_TARGET_DOC OFF CACHE BOOL ${AUTOPAS_BUILD_TARGET_DOC_DOC} FORCE )
endif ()

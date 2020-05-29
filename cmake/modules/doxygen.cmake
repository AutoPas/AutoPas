# make doc_doxygen optional if someone does not have / like doxygen
option(AUTOPAS_BUILD_TARGET_DOC "Create \"make doc_doxygen\" target (requires Doxygen)" ON)

# do nothing if nothing should be done
if (NOT AUTOPAS_BUILD_TARGET_DOC)
    return()
endif ()

# check if Doxygen is installed
find_package(Doxygen COMPONENTS dot OPTIONAL_COMPONENTS mscgen dia)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXY_CONF_DIR docs)
    set(DOXYGEN_IN ${DOXY_CONF_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    set(DOXY_MAIN_PAGE README.md)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    # note the option ALL which allows to build the docs together with the application
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
        # note the option ALL which allows to build the docs together with the application
        add_custom_target(
                doc_doxygen_md-flexible
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_MD-FLEXIBLE_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating documentation for md-flexible with Doxygen"
                VERBATIM
        )
    endif()

    message(STATUS "Doxygen configured")
else ()
    message(
        WARNING
            "Doxygen needs to be installed to generate the doxygen documentation, you might also have to install dot (graphviz)"
    )
    set(AUTOPAS_BUILD_TARGET_DOC OFF CACHE BOOL "" FORCE)
endif ()

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation (requires Doxygen)" ON)

# check if Doxygen is installed
find_package(Doxygen
        REQUIRED dot
        OPTIONAL_COMPONENTS mscgen dia)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXY_CONF_DIR docs)
    set(DOXYGEN_IN ${DOXY_CONF_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    set(DOXY_MAIN_PAGE ${DOXY_CONF_DIR}/mainpage.dox)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message(STATUS "Doxygen configured")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target(doc_doxygen
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
else (DOXYGEN_FOUND)
    message(STATUS "Doxygen needs to be installed to generate the doxygen documentation, you might also have to install doc (graphviz)")
endif (DOXYGEN_FOUND)
message("Generating ALL Documentation:")
find_package(Doxygen REQUIRED)
find_package(Python3 REQUIRED)
message(DEBUG "Python3 version: ${PYTHON3_VERSION}")
find_package(Sphinx REQUIRED)
#find breathe
# todo(s.schulz): The following does not seem to work. It returns non zero even
# though the package is installed and the correct Python interpreter is found.
# The python interpreter is in a non standard path and included using
# environment variables and the module was installed in the users directory
# (pip3 install --user breathe).
#execute_process(
#	COMMAND Python3::Interpreter -c "import breathe"
#	RESULT_VARIABLE BREATHE_NOT_FOUND
#	OUTPUT_QUIET
#	)
#if(NOT ${BREATHE_NOT_FOUND} EQUAL 0)
#	message( FATAL_ERROR "Python \"breathe\" package not found")
#endif()
# todo(s.schulz): Also find recommonmark
#do we need to make sure sphinx_rtd_theme is installed?

# Find all public headers
get_target_property(ALL_PUBLIC_HEADER_DIR ALL INTERFACE_INCLUDE_DIRECTORIES)
file(GLOB_RECURSE ALL_PUBLIC_HEADERS ${ALL_PUBLIC_HEADER_DIR})

# setup configuring of doxygen
# Caveat(s.schulz): Make sure there are not spaces within the paths.
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/example ${PROJECT_SOURCE_DIR}/tests/feature ${PROJECT_SOURCE_DIR}/tests/unit ${PROJECT_SOURCE_DIR}/README.md)
string(REPLACE ";" "\" \"" DOXYGEN_INPUT_DIR "${DOXYGEN_INPUT_DIR}")
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
set(DOXYGEN_INTERNAL_DOCS NO)
set(DOXYGEN_EXAMPLE_PATH ${PROJECT_SOURCE_DIR}/example)
# configure doxygen config file
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# create doxygen output directory
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

# todo(s.schulz): Add dependency on other files (Fortran module etc.)
add_custom_command( 
    OUTPUT ${DOXYGEN_INDEX_FILE}
    DEPENDS ${ALL_PUBLIC_HEADERS} ${ALL_EXE_NAME}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
    COMMENT "Generating ALL Docs (Doxygen)")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

# setup Sphinx
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

set(SPHINX_IN ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in)
set(SPHINX_OUT ${CMAKE_CURRENT_BINARY_DIR}/conf.py)
set(SPHINX_BREATHE_PROJECT_PATH ${DOXYGEN_OUTPUT_DIR}/xml)
configure_file(${SPHINX_IN} ${SPHINX_OUT} @ONLY)


# todo(s.schulz): Add dependencies to all .rst files
add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                  COMMAND
                  ${SPHINX_EXECUTABLE} -b html
                  -c "${CMAKE_CURRENT_BINARY_DIR}"
                  ${SPHINX_SOURCE} ${SPHINX_BUILD}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  DEPENDS 
                  ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
                  ${DOXYGEN_INDEX_FILE}
                  MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py.in
                  COMMENT "Generating ALL Docs (Sphinx)"
                  )

add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

include(GNUInstallDirs)
install(DIRECTORY ${SPHINX_BUILD}
        DESTINATION ${CMAKE_INSTALL_DOCDIR})
install(DIRECTORY ${DOXYGEN_OUTPUT_DIR}/html
        DESTINATION ${CMAKE_INSTALL_DOCDIR})

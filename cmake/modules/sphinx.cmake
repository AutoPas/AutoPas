find_package(Sphinx)

if (NOT SPHINX_FOUND)
    message(STATUS "Sphinx not found, not adding target doc")
    return()
endif()

set(SPHINX_SOURCE ${AUTOPAS_SOURCE_DIR}/docs/source/)
set(SPHINX_BUILD ${AUTOPAS_BINARY_DIR}/docs/sphinx)

add_custom_target(Sphinx
        COMMAND
        ${SPHINX_EXECUTABLE} -b html
        -Dbreathe_projects.AutoPas=${CMAKE_CURRENT_BINARY_DIR}/doc_doxygen/xml/
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${AUTOPAS_BINARY_DIR}
        COMMENT "Generating documentation with Sphinx")
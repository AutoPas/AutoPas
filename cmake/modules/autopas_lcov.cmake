if (NOT AUTOPAS_ENABLE_COVERAGE)
  return()
endif()

include(ExternalProject)

ExternalProject_Add(
    lcov
    PREFIX            ${CMAKE_CURRENT_BINARY_DIR}/lcov
    URL               ${PROJECT_SOURCE_DIR}/libs/lcov-2.1.zip
    URL_HASH          MD5=c276a95c6a3bf9af1d14a7350196a87d
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   make PREFIX=${CMAKE_CURRENT_BINARY_DIR}/lcov/install install
    BUILD_IN_SOURCE   True
)

# Update PATH to include the installed lcov
set(LCOV_EXECUTABLE ${CMAKE_BINARY_DIR}/lcov/install/bin/lcov)
set(GENHTML_EXECUTABLE ${CMAKE_BINARY_DIR}/lcov/install/bin/genhtml)
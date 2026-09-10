project(python)

find_package(
  Python
  COMPONENTS Interpreter Development
  REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(pypmt pypmt.cpp)
set_target_properties(pypmt PROPERTIES OUTPUT_NAME "pmt")
target_include_directories(pypmt PRIVATE ${PMT_INCLUDE_DIR})
target_link_libraries(pypmt PRIVATE pmt)

set(PYTHON_SITE_PACKAGES
    "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/pmt"
)

install(TARGETS pypmt DESTINATION "${PYTHON_SITE_PACKAGES}")
install(FILES __init__.py DESTINATION "${PYTHON_SITE_PACKAGES}")
install(FILES demo.py DESTINATION ${CMAKE_INSTALL_BINDIR})

# Setup site-packages in the build directory so that the CI can use it.
set(PYTHON_SITE_PACKAGES_BUILD "${CMAKE_CURRENT_BINARY_DIR}/site-packages/pmt")
add_custom_command(
  TARGET pypmt
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pypmt>
          ${PYTHON_SITE_PACKAGES_BUILD})
file(COPY __init__.py DESTINATION ${PYTHON_SITE_PACKAGES_BUILD})

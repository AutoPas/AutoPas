macro(add_sensor)
  ## Add sensor is a macro that allows to build and link a new sensor
  ## to PMT
  ## It accepts the following arguments:
  ## SENSOR_NAME : the name of the sensor
  ## HEADER : the main header file of the sensor
  ## SRC_FILES : the list of the files needed to create the sensor library
  ## INCLUDE_DIRECTORIES : the list of directories to include in the built
  ## LINK_LIBRARIES : the list of libraries to link to the sensor

  set(oneValueArgs SENSOR_NAME HEADER)
  set(multiValueArgs SRC_FILES LINK_LIBRARIES INCLUDE_DIRECTORIES)
  cmake_parse_arguments(ADD_SENSOR "" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  # Get current directory
  get_filename_component(SENSOR_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  # Set helper variables
  set(LIBRARY_NAME pmt-${ADD_SENSOR_SENSOR_NAME})

  # Build the necessary sensor files
  add_library(${LIBRARY_NAME} OBJECT ${ADD_SENSOR_SRC_FILES})
  target_link_libraries(${LIBRARY_NAME} PUBLIC ${ADD_SENSOR_LINK_LIBRARIES})
  if(ADD_SENSOR_INCLUDE_DIRECTORIES)
    target_include_directories(${LIBRARY_NAME}
                               PRIVATE ${ADD_SENSOR_INCLUDE_DIRECTORIES})
  endif()

  # Link to the global scope
  install(FILES ${ADD_SENSOR_HEADER} DESTINATION include/pmt)
  target_link_libraries(pmt PRIVATE pmt-${ADD_SENSOR_SENSOR_NAME})
  list(APPEND PMT_HEADER_FILES "${SENSOR_DIR}/${ADD_SENSOR_HEADER}")

  set(PMT_HEADER_FILES
      ${PMT_HEADER_FILES}
      PARENT_SCOPE)

endmacro()

project(pmt-cray)

include(AddSensor)
set(LINK_LIBRARIES Threads::Threads)

add_sensor(
  SENSOR_NAME
  "cray"
  HEADER
  Cray.h
  SRC_FILES
  Cray.cpp
  FilenamesHelper.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES})

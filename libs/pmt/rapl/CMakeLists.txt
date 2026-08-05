project(pmt-rapl)

set(LINK_LIBRARIES Threads::Threads)

add_sensor(
  SENSOR_NAME
  "rapl"
  HEADER
  Rapl.h
  SRC_FILES
  Rapl.cpp
  RaplImpl.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES})

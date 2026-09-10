project(pmt-xilinx)

include(AddSensor)

set(LINK_LIBRARIES Threads::Threads)

add_sensor(
  SENSOR_NAME
  "xilinx"
  HEADER
  Xilinx.h
  SRC_FILES
  Xilinx.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES})

project(pmt-amdsmi)

find_package(amd_smi REQUIRED)

include(AddSensor)

set(LINK_LIBRARIES Threads::Threads amd_smi)
set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR})

add_sensor(
  SENSOR_NAME
  "amdsmi"
  HEADER
  AMDSMI.h
  SRC_FILES
  AMDSMI.cpp
  AMDSMIImpl.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES}
  INCLUDE_DIRECTORIES
  ${INCLUDE_DIRECTORIES})

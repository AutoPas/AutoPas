project(pmt-tegra)

include(AddSensor)

set(LINK_LIBRARIES Threads::Threads)
set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR})

add_sensor(
  SENSOR_NAME
  "tegra"
  HEADER
  Tegra.h
  SRC_FILES
  Tegra.cpp
  TegraImpl.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES}
  INCLUDE_DIRECTORIES
  ${INCLUDE_DIRECTORIES})

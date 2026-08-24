project(pmt-likwid)

find_package(Likwid REQUIRED)

include(AddSensor)
set(LINK_LIBRARIES Threads::Threads ${LIKWID_LIBRARY})
set(INCLUDE_DIRECTORIES ${LIKWID_INCLUDE_DIR})

add_sensor(
  SENSOR_NAME
  "likwid"
  HEADER
  Likwid.h
  SRC_FILES
  Likwid.cpp
  LikwidImpl.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES}
  INCLUDE_DIRECTORIES
  ${INCLUDE_DIRECTORIES})

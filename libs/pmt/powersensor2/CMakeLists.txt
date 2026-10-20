project(pmt-powersensor2)

# Enable OpenMP
find_package(OpenMP REQUIRED)

FetchContent_Declare(
  powersensor2
  GIT_REPOSITORY https://git.astron.nl/RD/PowerSensor2
  GIT_TAG v3.0)

FetchContent_Populate(powersensor2)

set(POWERSENSOR2_SOURCE "${powersensor2_SOURCE_DIR}/PowerSensor.cc")

include(AddSensor)

set(LINK_LIBRARIES OpenMP::OpenMP_C)
set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR} ${powersensor2_SOURCE_DIR})

add_sensor(
  SENSOR_NAME
  "powersensor2"
  HEADER
  PowerSensor2.h
  SRC_FILES
  PowerSensor2.cpp
  ${POWERSENSOR2_SOURCE}
  LINK_LIBRARIES
  ${LINK_LIBRARIES}
  INCLUDE_DIRECTORIES
  ${INCLUDE_DIRECTORIES})

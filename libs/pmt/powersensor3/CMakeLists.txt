project(pmt-powersensor3)

# Enable OpenMP
find_package(OpenMP REQUIRED)

FetchContent_Declare(
  powersensor3
  GIT_REPOSITORY https://github.com/nlesc-recruit/PowerSensor3
  GIT_TAG v1.5.0)

FetchContent_Populate(powersensor3)

set(POWERSENSOR3_ROOT_SRC_DIR "${powersensor3_SOURCE_DIR}/host/src")
set(POWERSENSOR3_ROOT_INC_DIR "${powersensor3_SOURCE_DIR}/host/include")
set(POWERSENSOR3_SOURCE "${POWERSENSOR3_ROOT_SRC_DIR}/PowerSensor.cc"
                        "${POWERSENSOR3_ROOT_SRC_DIR}/sensors.cc")

include(AddSensor)

set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR} ${POWERSENSOR3_ROOT_INC_DIR})

add_sensor(
  SENSOR_NAME
  "powersensor3"
  HEADER
  PowerSensor3.h
  SRC_FILES
  PowerSensor3.cpp
  PowerSensor3Impl.cpp
  ${POWERSENSOR3_SOURCE}
  INCLUDE_DIRECTORIES
  ${INCLUDE_DIRECTORIES})

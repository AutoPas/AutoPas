project(pmt-rocm)

find_package(rocm_smi REQUIRED)
add_library(ROCM::smi INTERFACE IMPORTED)
target_include_directories(ROCM::smi INTERFACE ${ROCM_SMI_INCLUDE_DIRS})
target_link_libraries(ROCM::smi INTERFACE ${ROCM_SMI_LIBRARIES})

include(AddSensor)

set(LINK_LIBRARIES Threads::Threads ROCM::smi)
set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR})

add_sensor(
  SENSOR_NAME
  "rocm"
  HEADER
  ROCM.h
  SRC_FILES
  ROCM.cpp
  ROCMImpl.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES}
  INCLUDE_DIRECTORIES
  ${INCLUDE_DIRECTORIES})

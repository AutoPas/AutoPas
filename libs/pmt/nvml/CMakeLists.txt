project(pmt-nvml)

find_package(CUDAToolkit REQUIRED)

option(PMT_NVML_LEGACY_MODE "Build NVML in legacy mode" OFF)

if(CUDAToolkit_VERSION VERSION_LESS 12.0)
  set(PMT_NVML_LEGACY_MODE ON)
endif()

include(AddSensor)
set(LINK_LIBRARIES cudawrappers::nvml Threads::Threads)

add_sensor(
  SENSOR_NAME
  "nvml"
  HEADER
  NVML.h
  SRC_FILES
  NVML.cpp
  NVMLImpl.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES})

if(${PMT_NVML_LEGACY_MODE})
  target_compile_definitions(pmt-nvml PRIVATE PMT_NVML_LEGACY_MODE)
endif()

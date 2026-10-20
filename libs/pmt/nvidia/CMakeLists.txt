project(pmt-nvidia)

include(AddSensor)

find_package(CUDAToolkit REQUIRED)
set(LINK_LIBRARIES cudawrappers::nvml Threads::Threads)

add_sensor(
  SENSOR_NAME
  "nvidia"
  HEADER
  NVIDIA.h
  SRC_FILES
  NVIDIA.cpp
  LINK_LIBRARIES
  ${LINK_LIBRARIES})

if(PMT_BUILD_NVML)
  target_compile_definitions(${PROJECT_NAME} PRIVATE "PMT_BUILD_NVML")
endif()
if(PMT_BUILD_TEGRA)
  target_compile_definitions(${PROJECT_NAME} PRIVATE "PMT_BUILD_TEGRA")
endif()

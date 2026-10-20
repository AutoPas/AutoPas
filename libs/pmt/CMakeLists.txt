cmake_minimum_required(VERSION 3.17.5)

project(pmt)

# Set cmake module path
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# Set build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Default build type." FORCE)
endif()

# Set CXX Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set debug print info
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# Enable PThreads
find_package(Threads REQUIRED)

# Enable compile commands (needed for iwyu)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options to disable implementations
# that require external dependencies
option(PMT_BUILD_CRAY "" OFF)
option(PMT_BUILD_LIKWID "" OFF)
option(PMT_BUILD_NVML "" OFF)
option(PMT_BUILD_POWERSENSOR2 "" OFF)
option(PMT_BUILD_POWERSENSOR3 "" OFF)
option(PMT_BUILD_PYTHON "" OFF)
option(PMT_BUILD_RAPL "" OFF)
option(PMT_BUILD_ROCM "" OFF)
option(PMT_BUILD_AMDSMI "" OFF)
option(PMT_BUILD_TEGRA "" OFF)
option(PMT_BUILD_XILINX "" OFF)
option(PMT_BUILD_NVIDIA "" OFF)
option(PMT_BUILD_BINARY "" OFF)

include(FetchContent)

if(${PMT_BUILD_NVML} OR ${PMT_BUILD_TEGRA})
  set(PMT_BUILD_NVIDIA ON)
  find_package(CUDAToolkit REQUIRED)
  FetchContent_Declare(
    cudawrappers
    GIT_REPOSITORY https://github.com/nlesc-recruit/cudawrappers
    GIT_TAG main)
  FetchContent_MakeAvailable(cudawrappers)
endif()

# Config file
configure_file(${PROJECT_SOURCE_DIR}/cmake/config/pmt.h.in
               ${PROJECT_BINARY_DIR}/pmt.h)

# Installation
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
      ${CMAKE_BINARY_DIR}
      CACHE PATH "pmt install prefix" FORCE)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

include(GNUInstallDirs)

# Set toplevel source directory as include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Use this variable wherever <pmt.h> is included
set(PMT_INCLUDE_DIR "${PROJECT_BINARY_DIR}")

# PMT interface and common functionality
add_subdirectory(common)

# PMT library
add_library(pmt SHARED $<TARGET_OBJECTS:pmt-common>)

# The PMT instances below will append their public header file to this list, so
# that the header file structure can be replicated in the binary directory
# later.
set(PMT_HEADER_FILES "")

# Dummy
add_subdirectory(dummy)

# PMT components
set(PMT_COMPONENTS
    cray
    nvml
    tegra
    python
    rapl
    rocm
    amdsmi
    xilinx
    likwid
    nvidia
    powersensor2
    powersensor3)

foreach(component ${PMT_COMPONENTS})
  string(TOUPPER ${component} COMPONENT_UPPERCASE)
  if(${PMT_BUILD_${COMPONENT_UPPERCASE}})
    add_subdirectory(${component})
  endif()
endforeach()

# PMT executable
if(${PMT_BUILD_BINARY})
  add_subdirectory(bin)
endif()

# Public header file
set_target_properties(
  ${PROJECT_NAME} PROPERTIES PUBLIC_HEADER ${CMAKE_CURRENT_BINARY_DIR}/pmt.h)

# Replicate header file structure in binary directory
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/pmt")
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/pmt/common")
file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/common/PMT.h"
     "${PROJECT_BINARY_DIR}/pmt/common/PMT.h" SYMBOLIC)
foreach(PMT_HEADER_FILE ${PMT_HEADER_FILES})
  set(LINK_SRC "${CMAKE_CURRENT_SOURCE_DIR}/${PMT_HEADER_FILE}")
  get_filename_component(LINK_DST ${PMT_HEADER_FILE} NAME)
  file(CREATE_LINK ${LINK_SRC} "${PROJECT_BINARY_DIR}/pmt/${LINK_DST}" SYMBOLIC)
endforeach()

if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set_target_properties(${PROJECT_NAME} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                   ${PMT_INCLUDE_DIR})
else()
  # Install library
  target_include_directories(
    pmt PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-config
    COMPONENT ${PROJECT_NAME}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/)

  # Install cmake targets
  install(EXPORT ${PROJECT_NAME}-config
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
endif()

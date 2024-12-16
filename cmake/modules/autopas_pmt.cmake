option(pmt_ForceBundled "Do not look for an installed version, always used bundled." ON)

# RAPL sensor is disabled in PMT (PMT_BUILD_RAPL is OFF by default),
# however in AutoPas, it is set to ON by default whenever energy measurement is enabled.
set(PMT_BUILD_RAPL ON CACHE BOOL "RAPL is by default enabled when PMT is enabled" FORCE)

if (NOT ${pmt_ForceBundled})
    set(expectedVersion ${expectedVersion} QUIET)
    if (pmt_FOUND)
        message(STATUS "pmt - using installed version ${pmt_VERSION}")
        set_target_properties(pmt::pmt PROPERTIES "IMPORTED_GLOBAL" "TRUE")
        return()
    else()
        message(STATUS "pmt - no system version compatible to version ${expectedVersion} found")
        message(STATUS "pmt - if you want to use your version point the cmake variable pmt_DIR to the directory containing pmtConfig.cmake in order to find package")
    endif()
endif()

message(STATUS "pmt - using bundled version (commit e84518f7) and patch")

include(FetchContent)

FetchContent_Declare(
    pmt
    URL
        # pmt-master.zip contains commit e84518f7 (master as on Dec 16, 2024)
        ${AUTOPAS_SOURCE_DIR}/libs/pmt-master.zip
        URL_HASH MD5=79a6d8132db70bcd7151f70fdb17c50f
    
    # Applying patch with following changes:  
    # removal of 100 ms minimum time interval for energy samples, and
    # removing asynchronous energy measurement
    PATCH_COMMAND
        ${CMAKE_COMMAND} -E echo "Applying patch  ${CMAKE_BUILD_DIR} " &&
        ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR}/_deps/pmt-build git apply ${CMAKE_SOURCE_DIR}/libs/patches/patch-file-pmt-for-autopas.patch
)

FetchContent_MakeAvailable(pmt)

# sensors available in pmt that are not currently required in AutoPas
mark_as_advanced(
    PMT_BUILD_CRAY
    PMT_BUILD_NVML
    PMT_BUILD_POWERSENSOR2
    PMT_BUILD_POWERSENSOR3
    PMT_BUILD_PYTHON
    PMT_BUILD_ROCM
    PMT_BUILD_TEGRA
    PMT_BUILD_XILINX
    PMT_BUILD_NVIDIA
    PMT_BUILD_BINARY
)

if (IS_DIRECTORY "${pmt_SOURCE_DIR}")
    set_property(DIRECTORY ${pmt_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
endif()

target_compile_options(pmt PUBLIC -DCMAKE_INSTALL_PREFIX="./build")

target_include_directories(pmt SYSTEM PUBLIC "${pmt_SOURCE_DIR}")

option(pmt_ForceBundled "Do not look for an installed version, always used bundled." ON)
set(PMT_BUILD_RAPL ON CACHE BOOL "Enabling RAPL always when energy measurement is requested" FORCE)

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

message(STATUS "pmt - using bundled version 1.3.1 (commit e93666bd)")

include(FetchContent)

FetchContent_Declare(
    pmt
    URL
        # pmt_stable_1.0.zip contains pmt_1.3.1 version with changes including
        # removal of 100 ms minimum time interval for energy samples, and
        # removing asynchronous energy measurement
        ${AUTOPAS_SOURCE_DIR}/libs/pmt_stable_1.0.zip
        URL_HASH MD5=9a3e1d5c73bb672ff93ce9ffb3f277b7
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

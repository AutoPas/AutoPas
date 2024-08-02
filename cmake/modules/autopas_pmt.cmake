option(pmt_ForceBundled "Do not look for an installed version, always used bundled." ON)

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
        ${AUTOPAS_SOURCE_DIR}/libs/pmt-1.3.1.zip
        URL_HASH MD5=fce9602dd5c51677ae7309e7913e6939
)

FetchContent_MakeAvailable(pmt)

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

get_target_property(propval pmt INTERFACE_INCLUDE_DIRECTORIES)
message(" Directories: ${propval}")

target_include_directories(pmt SYSTEM PUBLIC "${pmt_SOURCE_DIR}")
#target_include_directories(pmt SYSTEM PUBLIC "${propval}")

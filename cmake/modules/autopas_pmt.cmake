option(AUTOPAS_ENABLE_ENERGY_MEASUREMENTS
        "Enables energy measurements and tuning using PMT library and enables the RAPL interface, by default, allowing for energy measurement on Intel and AMD hardware."
        ON
)

if (AUTOPAS_ENABLE_ENERGY_MEASUREMENTS)

    option(pmt_ForceBundled "Do not look for an installed version, always used bundled." ON)

    # RAPL sensor is disabled in PMT (PMT_BUILD_RAPL is OFF by default),
    # however in AutoPas, it is set to ON by default whenever energy measurement is enabled.
    set(PMT_BUILD_RAPL ON CACHE BOOL "RAPL is by default enabled when PMT is enabled" FORCE)

    # LIKWID can be enabled by setting the cmake option PMT_BUILD_LIKWID to ON. This is particularly useful on CoolMuc-4 where RAPL does not work.
    set(PMT_BUILD_LIKWID OFF CACHE BOOL "LIKWID is OFF by default when PMT is enabled, but user can change this." )

    if (NOT ${pmt_ForceBundled})
        if (pmt_FOUND)
            message(STATUS "pmt - using installed version ${pmt_VERSION}")
            set_target_properties(pmt::pmt PROPERTIES "IMPORTED_GLOBAL" "TRUE")
            return()
        else()
            message(STATUS "pmt not found!")
            message(STATUS "pmt - if you want to use your version point the cmake variable pmt_DIR to the directory containing pmtConfig.cmake in order to find package")
        endif()
    endif()

    message(STATUS "pmt - using bundled version (commit 7a56fa3a) and patch")

    include(FetchContent)

    FetchContent_Declare(
            pmt
            URL
            # pmt-master.zip contains commit 7a56fa3a (master as on Dec 19, 2024) with patch applied
            # The patch applies the following changes:
            # removal of 100 ms minimum time interval for energy samples, and
            # removing asynchronous energy measurement
            # better error handling
            # the patch can found under AutoPas/libs/patches/patch-file-pmt-for-autopas.patch
            ${AUTOPAS_SOURCE_DIR}/libs/pmt-master.zip
            URL_HASH MD5=3c60096bf151e11cde6efc6e5ede1195
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

    target_compile_options(pmt PUBLIC -w -DCMAKE_INSTALL_PREFIX="./build")

    target_include_directories(pmt SYSTEM PUBLIC "${pmt_SOURCE_DIR}")

endif ()

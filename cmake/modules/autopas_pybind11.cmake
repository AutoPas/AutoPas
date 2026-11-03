# Gets pybind11 by (in order of priority): an installed version via find_package, which is also what a pip build
# provides via pyproject.toml, or a download from upstream. It is not bundled in libs/, as it is only needed for the
# Python interface, which is off by default.
set(expectedVersion 3.1.0)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 ${expectedVersion} QUIET)

if (pybind11_FOUND)
    message(STATUS "pybind11 - using installed version ${pybind11_VERSION}: ${pybind11_DIR}")
else ()
    message(STATUS "pybind11 - no installed version >= ${expectedVersion} found; downloading v${expectedVersion}")
    include(FetchContent)
    FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG v${expectedVersion}
            GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(pybind11)
endif ()

# The Python interface autopas._core. Its pure Python part lives in python/autopas.
#
# This compiles md-flexible's sources into the module. Its configuration, YAML parser, particle generators, integrator
# and thermostat are what the Python interface builds on. md-flexible itself lives below examples/, which pulls in all
# of the example applications and requires all application libraries, so its target can not be reused here.

include(autopas_yaml-cpp)
include(autopas_pybind11)

if (NOT TARGET molecularDynamicsLibrary)
    message(FATAL_ERROR "AUTOPAS_BUILD_PYTHON_INTERFACE=ON requires AUTOPAS_BUILD_MD_LIBRARY=ON.")
endif ()
string(TOLOWER "${MD_FLEXIBLE_MODE}" MD_FLEXIBLE_MODE_lower)
if (NOT MD_FLEXIBLE_MODE_lower MATCHES "singlesite")
    message(FATAL_ERROR "The Python interface only supports MD_FLEXIBLE_MODE=SINGLESITE.")
endif ()
if (AUTOPAS_INCLUDE_MPI OR MD_FLEXIBLE_USE_MPI)
    message(FATAL_ERROR "The Python interface does not support MPI.")
endif ()

# Every functor costs a large instantiation of all traversals. The defaults leave out those which are only useful on
# specific hardware; Highway already picks the best instruction set.
option(AUTOPAS_PYTHON_FUNCTOR_AUTOVEC "Expose the auto vectorized Lennard-Jones functor to Python." ON)
option(AUTOPAS_PYTHON_FUNCTOR_HWY "Expose the Highway Lennard-Jones functor to Python." ON)
option(AUTOPAS_PYTHON_FUNCTOR_AVX "Expose the AVX Lennard-Jones functor to Python." OFF)
option(AUTOPAS_PYTHON_FUNCTOR_SVE "Expose the SVE Lennard-Jones functor to Python." OFF)
option(AUTOPAS_PYTHON_FUNCTOR_ATM "Expose the Axilrod-Teller-Muto functor to Python." ON)
option(
    AUTOPAS_PYTHON_CALC_GLOBALS
    "Native functors calculate potential energy and virial. Slightly slower, but the history of a simulation in Python is empty without."
    ON
)

set(MD_FLEXIBLE_DIR ${AUTOPAS_SOURCE_DIR}/examples/md-flexible)
file(GLOB_RECURSE MD_FLEXIBLE_SRC CONFIGURE_DEPENDS "${MD_FLEXIBLE_DIR}/src/*.cpp")
# The main loop is replaced by PySimulation, and with it everything which only that main loop uses.
list(FILTER MD_FLEXIBLE_SRC EXCLUDE REGEX "/src/(main|Simulation|GlobalVariableLogger)\\.cpp$")

file(GLOB_RECURSE PYTHON_INTERFACE_SRC CONFIGURE_DEPENDS "*.cpp" "*.h")

# NO_EXTRAS turns off link time optimization, which pybind11 adds by default. For a module made of all traversal
# instantiations this takes longer than compiling all of it.
pybind11_add_module(_core NO_EXTRAS ${PYTHON_INTERFACE_SRC} ${MD_FLEXIBLE_SRC})
if (NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|RelWithDebInfo|None)$")
    pybind11_strip(_core)
endif ()

target_include_directories(_core PRIVATE ${MD_FLEXIBLE_DIR} ${AUTOPAS_SOURCE_DIR})
# yaml-cpp comes first, because the include directories of autopas may contain another version of it, see src/autopas.
target_link_libraries(_core PRIVATE yaml-cpp::yaml-cpp autopas molecularDynamicsLibrary)
target_compile_definitions(
    _core
    PRIVATE
        $<$<BOOL:${AUTOPAS_PYTHON_FUNCTOR_AUTOVEC}>:MD_FLEXIBLE_FUNCTOR_AUTOVEC>
        $<$<BOOL:${AUTOPAS_PYTHON_FUNCTOR_HWY}>:MD_FLEXIBLE_FUNCTOR_HWY>
        $<$<BOOL:${AUTOPAS_PYTHON_FUNCTOR_AVX}>:MD_FLEXIBLE_FUNCTOR_AVX>
        $<$<BOOL:${AUTOPAS_PYTHON_FUNCTOR_SVE}>:MD_FLEXIBLE_FUNCTOR_SVE>
        $<$<BOOL:${AUTOPAS_PYTHON_FUNCTOR_ATM}>:MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC>
        $<$<BOOL:${AUTOPAS_PYTHON_CALC_GLOBALS}>:MD_FLEXIBLE_CALC_GLOBALS>
        AUTOPAS_PYTHON_VERSION="${AutoPas_VERSION_MAJOR}.${AutoPas_VERSION_MINOR}.${AutoPas_VERSION_PATCH}"
)

# Keep the path to the OpenMP runtime in the installed module, it is typically not in a default location on macOS.
set_target_properties(_core PROPERTIES INSTALL_RPATH_USE_LINK_PATH ON)

# scikit-build-core builds the wheel from what CMake installs (see pyproject.toml). The component keeps the install
# rules of the bundled dependencies out of the wheel. A plain CMake build keeps the module in the build tree.
if (SKBUILD)
    install(TARGETS _core LIBRARY DESTINATION autopas COMPONENT AutoPasPython)
endif ()

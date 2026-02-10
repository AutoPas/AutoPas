# Dependencies needed by the main library target
include(autopas_spdlog)
include(autopas_harmony)
include(autopas_eigen)
include(autopas_pmt)
include(autopas_pybind11)

add_subdirectory(fake-dlclose)

# the other libraries in this folder are needed from elsewhere (examples, test, ...)

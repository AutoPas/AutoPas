# Dependencies needed by the main library target
include(autopas_spdlog)
include(autopas_harmony)
include(autopas_eigen)

add_subdirectory(fake-dlclose)
include(autopas_nsimd)
include(autopas_xsimd)

# the other libraries in this folder are needed from elsewhere (examples, test, ...)

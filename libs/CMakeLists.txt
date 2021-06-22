# Dependencies needed by the main library target
include(autopas_spdlog)
include(harmony)
include(autopas_eigen)

add_subdirectory(fake-dlclose)

# the other libraries in this folder are needed from elsewhere (examples, test, ...)
# Enable FetchContent CMake module
include(FetchContent)

# Build spdlog and make the cmake targets available
FetchContent_Declare(
        kokkos
        URL
        # kokkos master:
        # https://github.com/kokkos/kokkos
        # kokkos commit c28a8b03288b185f846ddfb1b7c08213e12e2634:
        ${AUTOPAS_SOURCE_DIR}/libs/kokkos-master.zip
        URL_HASH MD5=87f59d0536f527158d5ad53bc133cd3f
)

# Disable stuff we don't need
#option(SPDLOG_BUILD_EXAMPLE "" OFF)
#option(SPDLOG_BUILD_TESTS "" OFF)
#option(SPDLOG_INSTALL "" OFF)

# hide options from ccmake
#mark_as_advanced(
#)

FetchContent_MakeAvailable(kokkos)

#if (IS_DIRECTORY "${kokkos_SOURCE_DIR}")
#    set_property(DIRECTORY ${kokkos_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
#endif ()

# Disable warnings
#target_compile_options(kokkos PRIVATE -w)

#get_target_property(propval kokkos INTERFACE_INCLUDE_DIRECTORIES)
#target_include_directories(kokkos SYSTEM PUBLIC "${propval}")

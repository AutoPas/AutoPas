# Gets googletest by (in order of priority): reusing targets a parent project provides, an installed
# version via find_package, or the bundled 1.17.0 copy in libs/googletest.
# Set googletest_ForceBundled=ON to force bundled.
option(googletest_ForceBundled "Ignore any provided/installed googletest and always use the bundled copy." OFF)
mark_as_advanced(googletest_ForceBundled)

if (NOT googletest_ForceBundled AND NOT AUTOPAS_FORCE_ALL_BUNDLED)
    # Path 1: reuse a googletest a parent project already defined.
    # AutoPas's tests link gmock, which pulls in gtest, so gmock is what has to be there.
    if (TARGET gmock OR TARGET GTest::gmock)
        message(STATUS "AutoPas: Reusing googletest provided by parent project")
        set(googletestResolved TRUE)
    elseif (TARGET gtest OR TARGET GTest::gtest)
        # A parent project provides gtest but not gmock. We cannot use an installed or bundled googletest to add gmock
        # without that pulling in a conflicting gtest, so we fail.
        message(
            FATAL_ERROR
                "A parent project provides gtest but no gmock, which AutoPas's tests link. Build its googletest with BUILD_GMOCK=ON, or set AUTOPAS_BUILD_TESTS=OFF."
        )
    else ()
        # Path 2: installed version
        set(expectedVersion 1.17.0)
        find_package(GTest ${expectedVersion} QUIET)
        if (GTest_FOUND)
            message(STATUS "gtest - using installed version ${GTest_VERSION}")
            autopas_promote_global(GTest::gtest)
            autopas_promote_global(GTest::gmock)
            set(googletestResolved TRUE)
        else ()
            message(STATUS "gtest - no installed version >= ${expectedVersion} found; using bundled copy")
        endif ()
    endif ()
endif ()

if (NOT googletestResolved)
    message(STATUS "gtest - using bundled version 1.17.0")
    find_package(Threads REQUIRED)

    option(INSTALL_GTEST "" OFF)

    # hide options from ccmake
    mark_as_advanced(
            BUILD_GMOCK
            INSTALL_GTEST
            GTEST_HAS_ABSL  # Only relevant when also using the Abseil library, which we don't.
    )

    # Prevent overriding the parent project's compiler/linker settings on Windows.
    # => Compiles gtest with correct mt(d)/md(d)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # EXCLUDE_FROM_ALL, such that gtest is not built by default (target `all`).
    add_subdirectory(${AUTOPAS_SOURCE_DIR}/libs/googletest ${CMAKE_CURRENT_BINARY_DIR}/googletest EXCLUDE_FROM_ALL)
endif ()

# googletest names its targets `gtest` and `gmock`, whereas an installed version only provides them
# namespaced. AutoPas's tests link the plain names, so make both names resolve to the same libraries.
autopas_alias_dependency(gtest GTest::gtest)
autopas_alias_dependency(gmock GTest::gmock)

# The installed package may have been built without gmock, as may the bundled copy if BUILD_GMOCK was
# turned off. Either way, AutoPas's tests cannot be built.
if (NOT TARGET gmock)
    message(
        FATAL_ERROR
            "The googletest used by AutoPas provides no gmock target, which AutoPas's tests link. Build it with BUILD_GMOCK=ON, or set AUTOPAS_BUILD_TESTS=OFF."
    )
endif ()
set(AUTOPAS_ENABLE_TREELITE_BASED_TUNING
    OFF
    CACHE
    BOOL "Enables tuning strategies that use Treelite. If enabled, Treelite will be built and linked."
)

if (AUTOPAS_ENABLE_TREELITE_BASED_TUNING)
  set(expectedVersion 4.6.1)
  set(bundledCommit b3930088a2b9d9385009ae95b96fd7fb55841d62)

  option(treelite_ForceBundled "Forcibly use the bundled version of Treelite (v${expectedVersion})" ON)

  if (NOT treelite_ForceBundled)
    find_package(treelite ${expectedVersion} CONFIG QUIET)
    if (treelite_FOUND)
      message(STATUS "Treelite - using installed system version")
      # Promote target to global visibility
      set_target_properties(treelite::treelite PROPERTIES "IMPORTED_GLOBAL" "TRUE")
      return()
    else()
      message(STATUS "Treelite - no compatible system version ${expectedVersion} found")
    endif()
  endif()

  # ForceBundled=ON or System version not found -> install bundled patched version
  message(STATUS "Treelite - using bundled version ${expectedVersion} (commit ${bundledCommit}, AutoPas patch)")

  # Enable FetchContent CMake module
  include(FetchContent)

  # Build treelite and make the cmake targets available
  FetchContent_Declare(
    treelite
    # treelite-4.6.1-patched.zip contains commit b3930088a2b9d9385009ae95b96fd7fb55841d62
    # with patch from libs/patches/patch-file-treelite-for-autopas.patch applied
    # The patch adds the prefix TREELITE_* to Treelite CMake options to prevent accidental name
    # collision with other options, and gets local archived Treelite dependencies
    URL ${AUTOPAS_SOURCE_DIR}/libs/treelite-${expectedVersion}-patched.zip
    URL_HASH MD5=60945421f670da753257f1bbf3eb1cd5
    EXCLUDE_FROM_ALL
  )

  # Tell patched Treelite where AutoPas stores offline dependency archives
  # (rapidjson-ab1842a2.zip, json-3.11.3.tar.xz, mdspan-0.6.0.zip).
  set(TREELITE_LOCAL_DEPS_DIR "${AUTOPAS_SOURCE_DIR}/libs" CACHE PATH "" FORCE)

  # Configure Treelite build options before MakeAvailable
  set(TREELITE_BUILD_CPP_TEST OFF CACHE BOOL "" FORCE)
  set(TREELITE_BUILD_DOXYGEN OFF CACHE BOOL "" FORCE)
  set(TREELITE_TEST_COVERAGE OFF CACHE BOOL "" FORCE)
  set(TREELITE_ENABLE_ALL_WARNINGS OFF CACHE BOOL "" FORCE)
  set(TREELITE_USE_SANITIZER OFF CACHE BOOL "" FORCE)
  set(TREELITE_HIDE_CXX_SYMBOLS OFF CACHE BOOL "" FORCE)
  set(TREELITE_DETECT_CONDA_ENV OFF CACHE BOOL "" FORCE)

  # Enabling OpenMP for Treelite may result in oversubscription if Treelite's GTIL predictor is configured to use more than 1 threads AND the predictor is used inside a parallelised region.
  if (AUTOPAS_OPENMP)
    set(TREELITE_USE_OPENMP ON CACHE BOOL "" FORCE)
  else ()
    set(TREELITE_USE_OPENMP OFF CACHE BOOL "" FORCE)
  endif()
  
  FetchContent_MakeAvailable(treelite)

  # Ensure the in-tree Treelite target exposes headers
  if (TARGET treelite)
  target_include_directories(treelite PUBLIC
    $<BUILD_INTERFACE:${treelite_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${treelite_BINARY_DIR}/include>
  )
  endif()

  # Create an alias for Treelite target
  if (TARGET treelite AND NOT TARGET treelite::treelite)
    add_library(treelite::treelite ALIAS treelite)
  endif()

endif()

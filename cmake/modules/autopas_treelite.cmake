# set(AUTOPAS_ENABLE_TREELITE_BASED_TUNING
#     OFF
#     CACHE
#     BOOL "Enables tuning strategies that use Treelite. If enabled, Treelite will be built (bundled) and linked."
# )
#
# if (AUTOPAS_ENABLE_TREELITE_BASED_TUNING)
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

  # System version not found -> install bundled version
  message(STATUS "Treelite - using bundled version ${expectedVersion} (commit ${bundledCommit})")

  # Enable FetchContent CMake module
  include(FetchContent)

  # Build treelite and make the cmake targets available
  FetchContent_Declare(
    treelite
    URL ${AUTOPAS_SOURCE_DIR}/libs/treelite-${expectedVersion}.zip
    URL_HASH MD5=b414be3a875a35de981723cbd01fce09
    EXCLUDE_FROM_ALL
  )

  # Configure Treelite build options before MakeAvailable
  set(BUILD_CPP_TEST OFF CACHE BOOL "" FORCE)
  set(BUILD_DOXYGEN OFF CACHE BOOL "" FORCE)
  set(TEST_COVERAGE OFF CACHE BOOL "" FORCE)
  set(ENABLE_ALL_WARNINGS OFF CACHE BOOL "" FORCE)
  set(USE_SANITIZER OFF CACHE BOOL "" FORCE)
  set(HIDE_CXX_SYMBOLS OFF CACHE BOOL "" FORCE)

  # Enabling this may result in oversubscription when AutoPas already uses OpenMP
  set(USE_OPENMP ON CACHE BOOL "" FORCE)
  
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

# endif()

file(
        GLOB_RECURSE
        MDFlex_SRC
        "src/*.cpp"
        "src/*.h"
)

add_executable(md-flexible ${MDFlex_SRC})

target_include_directories(md-flexible PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(autopas_yaml-cpp)
include(autopas_all)

option(MD_FLEXIBLE_USE_MPI "Enable MPI parallelism for md-flexible." OFF)

target_link_libraries(md-flexible PUBLIC autopas autopasTools molecularDynamicsLibrary yaml-cpp  ${ALL_LIB} $<$<BOOL:${MD_FLEXIBLE_USE_MPI}>:MPI::MPI_CXX>)

option(MD_FLEXIBLE_FUNCTOR_AUTOVEC "Compile AutoVec Functor." OFF)
option(MD_FLEXIBLE_FUNCTOR_AVX "If instruction set is available, compile AVX Functor." OFF)
option(MD_FLEXIBLE_FUNCTOR_SVE "If instruction set is available, compile SVE Functor." OFF)
option(MD_FLEXIBLE_FUNCTOR_HWY "Compile Highway Functor, which uses the best available SIMD instruction set supported by google highway." ON)
option(MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC "Compile auto vectorized Axilrod-Teller-Muto Functor." OFF)
option(MD_FLEXIBLE_FUNCTOR_ATM_HWY "Compile highway Axilrod-Teller-Muto Functor." OFF)
option(MD_FLEXIBLE_CALC_GLOBALS "Compiles md-flexible Functors with calculation of globals, including calculation of shift offset for LJ functors." OFF)

target_compile_definitions(
        md-flexible
        PUBLIC
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AUTOVEC}>:MD_FLEXIBLE_FUNCTOR_AUTOVEC>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AVX}>:MD_FLEXIBLE_FUNCTOR_AVX>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_SVE}>:MD_FLEXIBLE_FUNCTOR_SVE>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC}>:MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_ATM_HWY}>:MD_FLEXIBLE_FUNCTOR_ATM_HWY>
        $<$<BOOL:${MD_FLEXIBLE_CALC_GLOBALS}>:MD_FLEXIBLE_CALC_GLOBALS>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_HWY}>:MD_FLEXIBLE_FUNCTOR_HWY>
)

# --- display warnings if single-/multi-site mode is used with a functor type that has not been implemented ---
string(TOLOWER "${MD_FLEXIBLE_MODE}" MD_FLEXIBLE_MODE_lower)
if (MD_FLEXIBLE_MODE_lower MATCHES "multisite")
    if (MD_FLEXIBLE_FUNCTOR_AVX)
        message(WARNING "AVX Lennard-Jones functor has not been implemented for Multi-Site Molecules")
    endif()
    if (MD_FLEXIBLE_FUNCTOR_SVE)
        message(WARNING "SVE Lennard-Jones functor has not been implemented for Multi-Site Molecules")
    endif()
    if (MD_FLEXIBLE_FUNCTOR_HWY)
        message(WARNING "Highway functor has not been implemented for Multi-Site Molecules")
    endif()
    if (MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC)
        message(WARNING "Axilrod-Teller-Muto functor has not been implemented for Multi-Site Molecules")
    endif()
    if (MD_FLEXIBLE_FUNCTOR_ATM_HWY)
        message(WARNING "Axilrod-Teller-Muto functor has not been implemented for Multi-Site Molecules")
    endif()
endif()

# --- copy script files to build dir ---
file(
        GLOB_RECURSE SCRIPTS
        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
        "*.sh"
        "*.gp"
        "*.py"
        "*.yaml"
        "*.rule"
        "*.frule"
)

foreach (script ${SCRIPTS})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${script} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach (script)

# -----------------test-----------------

add_subdirectory(tests)

# add check for current target
# cmake-format: off

# randomly generated imbalanced scenario tested with LC only
add_test(
        NAME md-flexible.test-static
        COMMAND
        md-flexible
        --container linked
        --cutoff 1.
        --distribution-mean 5.0
        --distribution-stddeviation 2.0
        --data-layout soa
        --functor "Lennard-Jones"
        --iterations 10
        --particle-generator gauss
        --particles-total 10
        --traversal c08,sliced
        --verlet-rebuild-frequency 5
        --verlet-skin-radius 0
        --boundary-type none
        --deltaT 0.
        --no-end-config
        --no-progress-bar
        CONFIGURATIONS checkExamples
)


# The AVX example should only be tested if AVX is available.
set(HAS_AVX true)
if (NOT AUTOPAS_USE_AUTOVEC)
    # If no vectorization is used, don't enable the example.
    message(STATUS "No vectorization used, not adding md-flexible-avx.test-unaligned to ctest.")
    set(HAS_AVX false)
elseif (AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "SSE")
    # If vectorization is set to SSE: don't enable the example
    message(STATUS "Only SSE specified, not adding md-flexible-avx.test-unaligned to ctest.")
    set(HAS_AVX false)
elseif (AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "SVE")
    message(STATUS "SVE vectorization level, not adding md-flexible-avx.test-unaligned to ctest.")
    set(HAS_AVX false)
elseif (AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "NATIVE")
    # If Vectorization is set to native: we try to figure out the vectorization level.
    message(STATUS "Native vectorization level, trying to detect AVX capability.")
    try_compile(
            HAS_AVX
            ${CMAKE_BINARY_DIR}
            ${AUTOPAS_SOURCE_DIR}/cmake/tests/has_avx_test.cpp
            # We are abusing COMPILE_DEFINITIONS here. There is no other sane way of passing "-march=native" ...
            COMPILE_DEFINITIONS "-march=native"
            OUTPUT_VARIABLE HAS_AVX_ERROR
    )
else ()
    message(STATUS "Adding AVX example to test, as proper vectorization is manually specified.")
endif ()

if (HAS_AVX)
    set(TEST_FUNCTOR "Lennard-Jones-AVX")
    message(STATUS "AVX detected. Adding AVX test.")
    #dangerous, as lc_c04_combined_SoA uses unaligned values.
    add_test(
            NAME md-flexible-avx.test-unaligned
            COMMAND
            md-flexible
            --no-end-config
            --functor "${TEST_FUNCTOR}"
            --vectorization-pattern 1xVectorLength
            --deltaT 0
            --particle-generator uniform
            --log-level debug
            --traversal lc_c04_combined_SoA
            --particles-total 71
            --no-progress-bar
            CONFIGURATIONS checkExamples
    )
    set_tests_properties(
            md-flexible-avx.test-unaligned
            PROPERTIES
            ENVIRONMENT  "${LSAN_OPTIONS_STR}"
    )
else ()
  message(STATUS "AVX could not be detected. Not testing the AVX Functor.")
  set(TEST_FUNCTOR "Lennard-Jones")
endif ()

#stable, periodic particle grid tested with all configurations
add_test(
        NAME md-flexible.test-sim
        COMMAND
        md-flexible
        --container all
        --traversal all
        --cutoff 1.5
        --functor "${TEST_FUNCTOR}"
        --vectorization-pattern 1xVectorLength
        --tuning-phases 1
        --particle-generator grid
        --particles-per-dimension 10
        --particle-spacing 1.1225
        --verlet-rebuild-frequency 4
        --verlet-skin-radius 0.2
        --boundary-type periodic
        --deltaT 0.005
        --no-end-config
        --no-progress-bar
        CONFIGURATIONS checkExamples
)

add_test(
        NAME md-flexMeasurePerf
        COMMAND measurePerf.sh md-flexible -silent
        CONFIGURATIONS checkExamples
)
# cmake-format: on

set_tests_properties(
        md-flexible.test-static
        md-flexible.test-sim
        md-flexMeasurePerf
        # This adds lsanignoreList.txt to ignore MPI internal leaks
        PROPERTIES
            ENVIRONMENT  "${LSAN_OPTIONS_STR}"
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)

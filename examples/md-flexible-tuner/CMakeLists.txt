file(
        GLOB_RECURSE
        MDFlexTuner_SRC
        "../md-flexible/src/*.cpp"
        "../md-flexible/src/*.h"
        "src/*.cpp"
        "src/*.h"
)

list(FILTER MDFlexTuner_SRC EXCLUDE REGEX "[.][.]/md-flexible/src/main[.]cpp$")

add_executable(md-flexible-tuner ${MDFlexTuner_SRC})

target_include_directories(md-flexible-tuner PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../md-flexible")
target_include_directories(md-flexible-tuner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(autopas_yaml-cpp)
include(autopas_all)

option(MD_FLEXIBLE_USE_MPI "Enable MPI parallelism for md-flexible-tuner." OFF)

target_link_libraries(md-flexible-tuner PUBLIC autopas autopasTools molecularDynamicsLibrary yaml-cpp  ${ALL_LIB} $<$<BOOL:${MD_FLEXIBLE_USE_MPI}>:MPI::MPI_CXX>)

option(MD_FLEXIBLE_FUNCTOR_AUTOVEC "Compile AutoVec Functor without calculation of global values." OFF)
option(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS "Compile AutoVec Functor with calculation of global values." OFF)
option(MD_FLEXIBLE_FUNCTOR_AVX "If instruction set is available, compile AVX Functor without calculation of global values." ON)
option(MD_FLEXIBLE_FUNCTOR_SVE "If instruction set is available, compile SVE Functor without calculation of global values." ON)

target_compile_definitions(
        md-flexible-tuner
        PUBLIC
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AUTOVEC}>:MD_FLEXIBLE_FUNCTOR_AUTOVEC>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS}>:MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AVX}>:MD_FLEXIBLE_FUNCTOR_AVX>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_SVE}>:MD_FLEXIBLE_FUNCTOR_SVE>
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
endif()

# --- copy script files to build dir ---
file(
        GLOB_RECURSE SCRIPTS
        "../md-flexible/input/*.yaml"
)

message(STATUS "SCRIPTS=${SCRIPTS} CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")

foreach (script ${SCRIPTS})
    configure_file(${script} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach (script)

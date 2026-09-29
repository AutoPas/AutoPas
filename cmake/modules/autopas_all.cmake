# cmake module for adding ALL (A Loadbalancing Library)
#
# When enabled, it is taken from (in order of priority): a target a parent project provides, an installed version via 
# find_package, or a download of v0.9.4 from upstream. As an optional dependency for example code, where the native
# load balancing alternative (InvertedPressure) performs typically similar, we do not maintain a bundled version.

option(MD_FLEXIBLE_ENABLE_ALLLBL "Enable load balancing via ALL for MD-Flex" OFF)

if (NOT MD_FLEXIBLE_ENABLE_ALLLBL)
    message(STATUS "ALL load balancing library support disabled")
    set(ALL_LIB "")
    return()
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMD_FLEXIBLE_ENABLE_ALLLBL")
set(ALL_LIB "ALL")

# Path 1: reuse an ALL target a parent project already defined.
if (TARGET ALL OR TARGET ALL::ALL)
    message(STATUS "AutoPas: Reusing ALL provided by parent project")
    autopas_alias_dependency(ALL ALL::ALL)
    return()
endif ()

# Path 2: installed version
set(expectedVersion 0.9.4)
find_package(ALL ${expectedVersion} QUIET)
if (ALL_FOUND)
    message(STATUS "ALL - using installed system version ${ALL_VERSION}")
    autopas_promote_global(ALL)
    autopas_promote_global(ALL::ALL)
    autopas_alias_dependency(ALL ALL::ALL)
    return()
endif ()

# Path 3 + fallback: download from upstream. 
message(STATUS "ALL - no system version >= ${expectedVersion} found; downloading v${expectedVersion} from upstream")

include(FetchContent)
FetchContent_Declare(
        allfetch
        GIT_REPOSITORY https://gitlab.jsc.fz-juelich.de/SLMS/loadbalancing.git
        # tag v0.9.4, which is commit f71d9ef98bdc9561e9f4c29068e41f16098a5a2f
        GIT_TAG v0.9.4
        GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(allfetch)

autopas_alias_dependency(ALL ALL::ALL)

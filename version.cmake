# NOTE: copy-paste from cmake's configuration.

# AutoPas version number components.
set(AutoPas_VERSION_MAJOR 2)
set(AutoPas_VERSION_MINOR 0)
set(AutoPas_VERSION_PATCH 0)
#set(AutoPas_VERSION_RC 0)
set(AutoPas_VERSION_IS_DIRTY 0)

# Start with the full version number used in tags.  It has no dev info.
set(AutoPas_VERSION "${AutoPas_VERSION_MAJOR}.${AutoPas_VERSION_MINOR}.${AutoPas_VERSION_PATCH}")
if(DEFINED AutoPas_VERSION_RC)
    set(AutoPas_VERSION "${AutoPas_VERSION}-rc${AutoPas_VERSION_RC}")
endif()

# Releases define a small patch level.
if("${AutoPas_VERSION_PATCH}" VERSION_LESS 20000000)
    set(AutoPas_VERSION_IS_RELEASE 1)
else()
    set(AutoPas_VERSION_IS_RELEASE 0)
endif()

if(NOT AutoPas_VERSION_NO_GIT)
    # If this source was exported by 'git archive', use its commit info.
    set(git_info [==[$Format:%h %s$]==])

    # Otherwise, try to identify the current development source version.
    if(NOT git_info MATCHES "^([0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f]?[0-9a-f]?)[0-9a-f]* "
            AND EXISTS ${AutoPas_SOURCE_DIR}/.git)
        find_package(Git QUIET)
        if(GIT_FOUND)
            macro(_git)
                execute_process(
                        COMMAND ${GIT_EXECUTABLE} ${ARGN}
                        WORKING_DIRECTORY ${AutoPas_SOURCE_DIR}
                        RESULT_VARIABLE _git_res
                        OUTPUT_VARIABLE _git_out OUTPUT_STRIP_TRAILING_WHITESPACE
                        ERROR_VARIABLE _git_err ERROR_STRIP_TRAILING_WHITESPACE
                )
            endmacro()
        endif()
        if(COMMAND _git)
            # Get the commit checked out in this work tree.
            _git(log -n 1 HEAD "--pretty=format:%h %s" --)
            set(git_info "${_git_out}")
        endif()
    endif()

    # Extract commit information if available.
    if(git_info MATCHES "^([0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f]?[0-9a-f]?)[0-9a-f]* (.*)$")
        # Have commit information.
        set(git_hash "${CMAKE_MATCH_1}")
        set(git_subject "${CMAKE_MATCH_2}")

        # If this is not the exact commit of a release, add dev info.
        if(NOT "${git_subject}" MATCHES "^[Cc][Mm]ake ${AutoPas_VERSION}$")
            set(AutoPas_VERSION "${AutoPas_VERSION}-${git_hash}")
        endif()

        # If this is a work tree, check whether it is dirty.
        if(COMMAND _git)
            _git(update-index -q --refresh)
            _git(diff-index --name-only HEAD --)
            if(_git_out)
                set(AutoPas_VERSION_IS_DIRTY 1)
            endif()
        endif()
    else()
        # No commit information.
        if(NOT AutoPas_VERSION_IS_RELEASE)
            # Generic development version.
            set(AutoPas_VERSION "${AutoPas_VERSION}-git")
        endif()
    endif()
endif()

# Extract the version suffix component.
if(AutoPas_VERSION MATCHES "-(.*)$")
    set(AutoPas_VERSION_SUFFIX "${CMAKE_MATCH_1}")
else()
    set(AutoPas_VERSION_SUFFIX "")
endif()
if(AutoPas_VERSION_IS_DIRTY)
    set(AutoPas_VERSION ${AutoPas_VERSION}-dirty)
endif()

message(STATUS "AutoPas_Version: ${AutoPas_VERSION}")

configure_file(
        "${AUTOPAS_SOURCE_DIR}/src/autopas/Version.h.in"
        "${AUTOPAS_BINARY_DIR}/src/autopas/Version.h"
)

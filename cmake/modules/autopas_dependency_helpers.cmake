# Helpers for handling third-party dependency modules

# Umbrella switch for the per-dependency <dep>_ForceBundled options
option(AUTOPAS_FORCE_ALL_BUNDLED "Ignore all provided/installed dependencies and always use the bundled copies." OFF)


# A dependency can appear under a plain target name and/or a namespaced `ns::name`, depending on whether
# it comes from a parent project, an installed package (find_package) or the bundled copy. These helpers
# let each module reuse whichever target is available and expose both names, so consumers can link either
# and only ever a single copy of the dependency ends up in the build.


# Promote an imported target to global visibility so it can be linked from sibling directory scopes
# (e.g. a dependency resolved in libs/ but linked from src/). No-op for non-imported, already-global, or
# non-existent targets. Only call this for a target created in the current scope (i.e. by find_package
# here): promoting a target created in an outer scope, such as one a parent project provided, errors in CMake.
function(autopas_promote_global target)
    if (NOT TARGET ${target})
        return()
    endif ()
    # an ALIAS cannot take properties; the underlying target is promoted under its own name instead
    get_target_property(aliased ${target} ALIASED_TARGET)
    if (aliased)
        return()
    endif ()
    get_target_property(imported ${target} IMPORTED)
    get_target_property(global ${target} IMPORTED_GLOBAL)
    if (imported AND NOT global)
        set_target_properties(${target} PROPERTIES IMPORTED_GLOBAL TRUE)
    endif ()
endfunction()

# Make both `plain` and `namespaced` (a `ns::name`) resolve to the same library by creating whichever is
# missing from the one that exists. A namespaced name can only be a real ALIAS; the plain name instead
# becomes a forwarding INTERFACE target, which (unlike an ALIAS) also works when the existing target is a
# non-global imported target from an outer scope (e.g. provided by a parent project). No-op if both, or neither, already exist.
function(autopas_alias_dependency plain namespaced)
    if (TARGET ${plain} AND NOT TARGET ${namespaced})
        # Avoid aliasing an alias
        get_target_property(aliased ${plain} ALIASED_TARGET)
        if (aliased)
            add_library(${namespaced} ALIAS ${aliased})
        else ()
            add_library(${namespaced} ALIAS ${plain})
        endif ()
    elseif (TARGET ${namespaced} AND NOT TARGET ${plain})
        add_library(${plain} INTERFACE)
        target_link_libraries(${plain} INTERFACE ${namespaced})
    endif ()
endfunction()

# Warn if the copy of `dependency` a parent project provides is older than `minVersion`, print message if this cannot
# be determined. AutoPas reuses such a target as it is, whereas an installed version has the minimum enforced by 
# find_package. Targets carry no standard
# version, so it is taken from `<package>_VERSION`, which a parent's find_package leaves visible in AutoPas' scope,
# or else, for a parent that vendors the dependency, from the VERSION property upstream sets on its library target,
# looked up on the target names passed after the first three arguments. If neither is available, the version cannot
# be checked, which is only reported as a status message.
function(autopas_warn_if_parent_version_too_old dependency package minVersion)
    # Option 1: The target provides a `<package>_VERSION`
    set(version "${${package}_VERSION}")
    # Option 2: Get it from the target property
    if (NOT version)
        foreach (target IN LISTS ARGN)
            if (TARGET ${target})
                get_target_property(type ${target} TYPE)
                # Reading VERSION from an INTERFACE library is a configure error before CMake 3.19
                if (NOT type STREQUAL "INTERFACE_LIBRARY")
                    get_target_property(version ${target} VERSION)
                    if (version)
                        break()
                    endif ()
                endif ()
            endif ()
        endforeach ()
    endif ()

    if (NOT version)
        message(STATUS "${dependency} - cannot determine the version of the parent project's copy, so it is not "
                       "checked against AutoPas' minimum of ${minVersion}")
    elseif (version VERSION_LESS minVersion)
        message(
            WARNING
                "${dependency} - the parent project provides version ${version}, but AutoPas requires at least "
                "${minVersion}. AutoPas reuses it regardless, so it may fail to compile or misbehave; if so, upgrade "
                "the parent project's ${dependency}."
        )
    else ()
        message(STATUS "${dependency} - parent project provides version ${version}")
    endif ()
endfunction()

# To avoid violating the one defintion rule, error if forcing bundled where a parent project provides their own version
# of the dependency.
function(autopas_error_if_forced_bundled_collides dependency forceOption)
    foreach (target IN LISTS ARGN)
        if (TARGET ${target})
            message(
                FATAL_ERROR
                    "${dependency} - the bundled copy is forced (${forceOption} or AUTOPAS_FORCE_ALL_BUNDLED), but a "
                    "parent project already provides ${dependency} as the target '${target}', and AutoPas cannot "
                    "build its own copy alongside it. Unset the option to reuse the parent's ${dependency}; if its "
                    "version is incompatible with AutoPas, align the versions instead."
            )
        endif ()
    endforeach ()
endfunction()
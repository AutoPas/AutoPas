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
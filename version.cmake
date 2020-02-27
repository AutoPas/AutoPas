set(AutoPas_VERSION_MAJOR 0 CACHE STRING "AutoPas major version number." FORCE)
set(AutoPas_VERSION_MINOR 1 CACHE STRING "AutoPas minor version number." FORCE)
set(AutoPas_VERSION_PATCH 0 CACHE STRING "AutoPas patch version number." FORCE)
set(
    AutoPas_VERSION
    ${AutoPas_VERSION_MAJOR}.${AutoPas_VERSION_MINOR}.${AutoPas_VERSION_PATCH}
    CACHE STRING "AutoPas version" FORCE
)
mark_as_advanced(
    AutoPas_VERSION
    AutoPas_VERSION_MAJOR
    AutoPas_VERSION_MINOR
    AutoPas_VERSION_PATCH
)

# - Try to find likwid
# This module tries to find the likwid library on your system
#
# Once done this will define
#  LIKWID_FOUND       - system has likwid
#  LIKWID_INCLUDE_DIR - the likwid include directory
#  LIKWID_LIBRARY     - link these to use likwid

find_package(PackageHandleStandardArgs)

find_library(LIKWID_LIBRARY likwid ENV LD_LIBRARY_PATH)

get_filename_component(LIKWID_LIB_DIR ${LIKWID_LIBRARY} PATH)
get_filename_component(_LIKWID_INC_DIR ${LIKWID_LIB_DIR}/../include ABSOLUTE)

# search for headers relative to the library
find_path(LIKWID_INCLUDE_DIR likwid.h PATHS ${_LIKWID_INC_DIR})

find_package_handle_standard_args(Likwid DEFAULT_MSG LIKWID_LIBRARY
                                  LIKWID_INCLUDE_DIR)

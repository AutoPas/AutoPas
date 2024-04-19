message(STATUS "libxerces - using system version")

# Load the XercesC package

# TODO: try to use a bundled version of xerces-c instead of installing it globally
find_package(XercesC 3.2.0 REQUIRED)

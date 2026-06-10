project(voro)
cmake_minimum_required(VERSION 2.8)
include(GNUInstallDirs)

include_directories (include)
add_subdirectory (src)

install (DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.hh")

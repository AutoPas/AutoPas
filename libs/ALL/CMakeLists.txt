cmake_minimum_required (VERSION 3.14)

option(CM_ALL_VTK_OUTPUT "VTK output routine" OFF)
option(CM_ALL_VORONOI "Voronoi-based loadbalancing scheme (uses bundled Voro++)" OFF)
option(CM_ALL_EXTERNAL_VORONOI "Voronoi-based loadbalancing scheme (requires external Voro++)" OFF)
option(CM_ALL_FORTRAN "VTK output routine" OFF)
option(CM_ALL_DEBUG "Enable debug information" OFF)
option(CM_ALL_USE_F08 "Use some Fortran 2008 features (mpi_f08)" OFF)
option(CM_ALL_FORTRAN_ERROR_ABORT "Abort on any error when using the
Fortran interface" OFF)
option(CM_ALL_TESTS "Enables test suite" OFF)
option(CM_ALL_TESTS_INTEGRATION "Enables integration tests" OFF)
option(CM_ALL_AUTO_DOC "Enables creation of auto-documentation")

# leading and trailing whitespace should be silently ignored
# thanks to an old FindMPI script
# cmake_policy(SET CMP0004 OLD)

enable_language(CXX)
enable_language(C)

if (CM_ALL_USE_F08)
    if (NOT CM_ALL_FORTRAN)
        message(FATAL_ERROR "Flag for Fortran 08 MPI interface set, while Fortran interface is not enabled!")
    endif ()
endif ()

if (CM_ALL_FORTRAN)
    enable_language(Fortran)
# todo(s.schulz): This should be extended for other compilers, if they
#                 support standard checking. The Intel Fortran compiler
#                 does not.
    if(CM_ALL_USE_F08)
        add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-std=f2008>)
    else()
        add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-std=f2003>)
    endif()
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-fbacktrace>)
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:Intel>>:-traceback>)
endif (CM_ALL_FORTRAN)

# set standards compliance flags
# set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (CM_ALL_FORTRAN)
    add_compile_options("$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-Wall;-Wextra>")
endif (CM_ALL_FORTRAN)
add_compile_options("$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:GNU>>:-Wall;-Wextra>")
add_compile_options("$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:-Wall;-Wextra>")

project(ALL
    VERSION "0.9.3"
    DESCRIPTION "A Loadbalacing Library"
    HOMEPAGE_URL "http://slms.pages.jsc.fz-juelich.de/websites/all-website/")
set(ALL_SOVERSION 0)

include(GNUInstallDirs)

# add custom find-scripts to module path
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if (CM_ALL_FORTRAN)
    find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
else()
    find_package(MPI REQUIRED COMPONENTS C CXX)
endif (CM_ALL_FORTRAN)

# Prepare the Config.cmake.in
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Config.cmake.in
"@PACKAGE_INIT@\n
\n
include(CMakeFindDependencyMacro)\n
find_dependency(MPI)\n"
    )

if(CM_ALL_VTK_OUTPUT)
    message(STATUS "Using VTK output")
    find_package(VTK)
    if(NOT VTK_FOUND)
        message(FATAL_ERROR "VTK not found, help CMake to find it by setting VTK_LIBRARY and VTK_INCLUDE_DIR")
    else()
        if (VTK_MAJOR_VERSION GREATER_EQUAL 9)
            find_package(VTK REQUIRED COMPONENTS
                CommonCore
                CommonDataModel
                FiltersProgrammable
                IOParallelXML
                IOXML
                ParallelMPI
            )
        else()
            find_package(VTK REQUIRED COMPONENTS
                vtkCommonCore
                vtkCommonDataModel
                vtkFiltersProgrammable
                vtkIOParallelXML
                vtkIOXML
                vtkParallelMPI
            )
        endif()
    endif()
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/Config.cmake.in "find_dependency(VTK REQUIRED)\n")
endif(CM_ALL_VTK_OUTPUT)

if(CM_ALL_EXTERNAL_VORONOI)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/Config.cmake.in "find_dependency(VORO REQUIRED)\n")
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindVORO.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ALL)
endif(CM_ALL_EXTERNAL_VORONOI)

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/Config.cmake.in
"include(\"\${CMAKE_CURRENT_LIST_DIR}/ALLTargets.cmake\")\n
check_required_components(ALL)\n"
    )

if (CM_ALL_DEBUG)
    message("Using ALL debug information")
    add_compile_definitions("ALL_DEBUG_ENABLED")
endif(CM_ALL_DEBUG)

if(CM_ALL_FORTRAN_ERROR_ABORT)
    add_compile_definitions("ALL_FORTRAN_ERROR_ABORT")
endif()

if(CM_ALL_USE_F08 AND MPI_Fortran_HAVE_F08_MODULE)
    add_compile_definitions("ALL_USE_F08")
elseif(CM_ALL_USE_F08 AND NOT MPI_Fortran_HAVE_F08_MODULE)
    message(FATAL_ERROR "Fortran 2008 MPI chosen, but MPI installation does not support the Fortran 2008 modules")
endif(CM_ALL_USE_F08 AND MPI_Fortran_HAVE_F08_MODULE)

if(CM_ALL_VORONOI AND CM_ALL_EXTERNAL_VORONOI)
    message(FATAL_ERROR "You cannot use CM_ALL_VORONOI and CM_ALL_EXTERNAL_VORONOI together, pick one")
elseif(CM_ALL_VORONOI)
    message(STATUS "compiling voro++ version in contrib/voro++")
    add_subdirectory(contrib/voro++)
elseif(CM_ALL_EXTERNAL_VORONOI)
    find_package(VORO REQUIRED)
endif()
if(CM_ALL_TESTS_INTEGRATION)
    if(NOT CM_ALL_TESTS)
        message(WARNING "Integration tests enabled, but testing itself not. No tests will be generated!")
    endif()
endif()

set(STB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/contrib/stb)

add_subdirectory(src)
add_subdirectory(example)

if (CM_ALL_TESTS)
    enable_testing()
    include(FindPackageHandleStandardArgs)
    add_subdirectory(tests)
endif(CM_ALL_TESTS)

if(CM_ALL_AUTO_DOC)
    add_subdirectory("docs")
endif(CM_ALL_AUTO_DOC)


# Export targets and create package
install(EXPORT ALLTargets
    NAMESPACE ALL::
    FILE ALLTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ALL
    )

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_BINARY_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/ALLConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ALL
    )
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ALLConfigVersion.cmake"
    VERSION "${ALL_VERSION}"
    COMPATIBILITY SameMajorVersion)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ALLConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ALLConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ALL)

# Create Makefile.libALL.in for inclusion in external Makefile projects
# todo(s.schulz): The include directories and libraries should be auto
# detected... but that is not easily done right now. Ideally, this should also
# provide the very long VTK linker flags.
set(ALL_MAKEFILE_INC ${CMAKE_CURRENT_BINARY_DIR}/Makefile.libALL.in)
get_target_property(_liball_defines ALL INTERFACE_COMPILE_DEFINITIONS)
list(TRANSFORM _liball_defines PREPEND " -D")
list(JOIN _liball_defines " " _liball_defines_string)
file(WRITE ${ALL_MAKEFILE_INC} "LIBALL_DEFINES := ${_liball_defines_string}\n")
file(APPEND ${ALL_MAKEFILE_INC} "LIBALL_INCLUDE := -I\$(ALL_DIR)/${CMAKE_INSTALL_INCLUDEDIR}\n")
file(APPEND ${ALL_MAKEFILE_INC} "LIBALL_LIB :=\n")
if(CM_ALL_FORTRAN)
    file(APPEND ${ALL_MAKEFILE_INC} "LIBALL_LIB += \$(ALL_DIR)/${CMAKE_INSTALL_LIBDIR}/libALL_fortran.a \$(ALL_DIR)/${CMAKE_INSTALL_LIBDIR}/libALL.a -lstdc++\n")
    file(APPEND ${ALL_MAKEFILE_INC} "LIBALL_INCLUDE += -I\$(ALL_DIR)/${CMAKE_INSTALL_LIBDIR}\n")
endif()
if(CM_ALL_VORONOI)
    file(APPEND ${ALL_MAKEFILE_INC} "LIBALL_LIB += \$(ALL_DIR)/${CMAKE_INSTALL_LIBDIR}/libvoro.a\n")
elseif(CM_ALL_EXTERNAL_VORONOI)
    file(APPEND ${ALL_MAKEFILE_INC} "LIBALL_LIB += ${VORO_LIBRARIES}\n")
endif()
install(FILES ${ALL_MAKEFILE_INC} DESTINATION ${CMAKE_INSTALL_LIBDIR}/make)

# vim: sw=4 ts=4 et

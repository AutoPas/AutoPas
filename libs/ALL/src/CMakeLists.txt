set(ALL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(ALL_HEADER_FILES
    ${ALL_INCLUDE_DIR}/ALL_CustomExceptions.hpp
    ${ALL_INCLUDE_DIR}/ALL_Defines.h
    ${ALL_INCLUDE_DIR}/ALL_ForceBased.hpp
    ${ALL_INCLUDE_DIR}/ALL_Functions.hpp
    ${ALL_INCLUDE_DIR}/ALL_Histogram.hpp
    ${ALL_INCLUDE_DIR}/ALL.hpp
    ${ALL_INCLUDE_DIR}/ALL_LB.hpp
    ${ALL_INCLUDE_DIR}/ALL_Point.hpp
    ${ALL_INCLUDE_DIR}/ALL_Staggered.hpp
    ${ALL_INCLUDE_DIR}/ALL_Tensor.hpp
    ${ALL_INCLUDE_DIR}/ALL_Voronoi.hpp)

add_library (ALL ALL.cpp ALL_fortran.cpp ${ALL_HEADER_FILES})
set_target_properties(ALL PROPERTIES SOVERSION ${ALL_SOVERSION})

if(CM_ALL_VTK_OUTPUT)
    target_include_directories(ALL PUBLIC ${VTK_INCLUDE_DIRS})
    target_link_libraries(ALL PUBLIC ${VTK_LIBRARIES})
    target_compile_definitions(ALL PUBLIC ALL_VTK_OUTPUT)
endif(CM_ALL_VTK_OUTPUT)
if(CM_ALL_VORONOI)
    target_link_libraries(ALL PUBLIC voro)
    target_compile_definitions(ALL PUBLIC ALL_VORONOI_ACTIVE)
elseif(CM_ALL_EXTERNAL_VORONOI)
    target_link_libraries(ALL PUBLIC VORO::voro++)
    target_compile_definitions(ALL PUBLIC ALL_VORONOI_ACTIVE)
endif(CM_ALL_VORONOI)

target_include_directories(ALL PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(ALL PUBLIC ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
target_include_directories (ALL
    PUBLIC
    "$<BUILD_INTERFACE:${ALL_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

# For inclusion in Fortran code, explicit linking against the C++ standard
# library is necessary.
target_link_libraries(ALL INTERFACE stdc++)

install(TARGETS ALL
    EXPORT ALLTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${ALL_HEADER_FILES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if (CM_ALL_FORTRAN)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../modules)
    add_library (ALL_fortran ALL_module.F90)
    set_property(TARGET ALL_fortran PROPERTY LINKER_LANGUAGE Fortran)
    set_target_properties(ALL_fortran PROPERTIES SOVERSION ${ALL_SOVERSION})

    target_include_directories(ALL_fortran PUBLIC ${MPI_Fortran_INCLUDE_PATH})
    target_include_directories(ALL_fortran
        INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>")
    target_link_libraries(ALL_fortran PUBLIC ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
    target_link_libraries(ALL_fortran PUBLIC ALL)

    install(TARGETS ALL_fortran
        EXPORT ALLTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/ DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif(CM_ALL_FORTRAN)

add_executable(Amalgamate Amalgamate.c)
target_include_directories(Amalgamate PRIVATE ${STB_INCLUDE_DIR})
set(AMALGAMATED_FILENAME ALL_Amalgamated.hpp)
add_custom_command(
    TARGET ALL
    POST_BUILD
    COMMAND $<TARGET_FILE:Amalgamate> ${ALL_INCLUDE_DIR}/ALL.hpp ${ALL_INCLUDE_DIR} > ${AMALGAMATED_FILENAME}
    )
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${AMALGAMATED_FILENAME}
    TYPE INCLUDE
    )
# vim: sw=4 ts=4 et

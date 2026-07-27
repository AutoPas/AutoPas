# using Boost test environment for tests
find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)

# each test listed in alphabetical order
foreach( PROG
         all_class_creation
         border_shift_1d
         point_access
         point_arithmetic
         point_cross
         point_dot
         point_tetrahedron
       )

    file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
    add_executable(unit_${PROG} ${${PROG}_SOURCES})
    target_link_libraries(unit_${PROG} ALL Boost::unit_test_framework)
    target_compile_definitions(unit_${PROG} PRIVATE BOOST_TEST_DYN_LINK)
    add_test(unit_${PROG} unit_${PROG})
    set_tests_properties(unit_${PROG} PROPERTIES LABELS "ALL;unit")
endforeach(PROG) 

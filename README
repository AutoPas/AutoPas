AUTOPAS Library

==========BUILD===========
build instructions:

mkdir build
cd build
cmake ..
make

==========TESTS===========
* to run tests:
	make test
  or using the ctest environment:
	ctest
to get verbose output:
	ctest --verbose
to run specific tests:
  use the --gtest_filter variable:
	./tests/testAutopas/runTests --gtest_filter=ArrayMathTest.testAdd*
  use the GTEST_FILTER environment variable:
	GTEST_FILTER="ArrayMathTest.testAdd*" ctest --verbose



=======DOCUMENTATION======
* to build the documentation:
    make doc_doxygen
* requirements:
    doxygen

==========CODING==========
* We use google code style.
* code style can be build with
	make clangformat
* requirements:
	clang-format

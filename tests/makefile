INCLUDEFLAGS = -I../src
CXXFLAGS =	-std=c++11 -O3 -g -Wall $(INCLUDEFLAGS) #-fmessage-length=0

LIBS = -lgtest -lgtest_main -lpthread

default:
	g++ allTests.cpp myMoleculeTest.cpp  $(INCLUDEFLAGS) $(LIBS) -o runTests
	./runTests

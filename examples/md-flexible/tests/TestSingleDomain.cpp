/**
 * @file TestSingleDomain.cpp
 * @author J. Körner
 * @date 19/05/21
 */
#include "TestSingleDomain.h"

#include "src/domainDecomposition/SingleDomain.h"

TestSingleDomain::TestSingleDomain() : AutoPasTestBase() { }

TEST_F(TestSingleDomain, testGetLocalDomain){
	std::vector<double> globalBoxMin = { 1.0, 1.0, 1.0 };
	std::vector<double> globalBoxMax = { 10.0, 10.0, 10.0};

	SingleDomain domainDecomposition(0,"", 3, globalBoxMin, globalBoxMax);

	EXPECT_EQ(globalBoxMin, domainDecomposition.getLocalBoxMin());
	EXPECT_EQ(globalBoxMax, domainDecomposition.getLocalBoxMax());
}

/**
 * @file TestNoDecomposition.cpp
 * @author J. Körner
 * @date 19/05/21 */
#include "TestNoDecomposition.h"

#include "src/domainDecomposition/NoDecomposition.h"

TEST_F(TestNoDecomposition, testGetLocalDomain) {
  std::vector<double> globalBoxMin = {1.0, 1.0, 1.0};
  std::vector<double> globalBoxMax = {10.0, 10.0, 10.0};

  char** argv;

  NoDecomposition domainDecomposition(0, argv, 3, globalBoxMin, globalBoxMax);

  EXPECT_EQ(globalBoxMin, domainDecomposition.getLocalBoxMin());
  EXPECT_EQ(globalBoxMax, domainDecomposition.getLocalBoxMax());
}

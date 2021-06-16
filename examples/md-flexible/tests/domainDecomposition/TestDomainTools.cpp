/**
 * @file TestDomainTools.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "TestDomainTools.h"

#include "src/domainDecomposition/DomainTools.h"

TEST_F(TestDomainTools, testIsInsideDomain) {
  std::vector<double> globalBoxMin = {1.0, 1.0, 1.0};
  std::vector<double> globalBoxMax = {10.0, 10.0, 10.0};

  std::vector<double> inside = {5.0, 5.0, 5.0};
  std::vector<double> onBoundary = {1.0, 1.0, 10.0};
  std::vector<double> outside = {0.0, 0.0, 0.0};

  EXPECT_EQ(true, DomainTools::isInsideDomain(inside, globalBoxMin, globalBoxMax));
  EXPECT_EQ(true, DomainTools::isInsideDomain(onBoundary, globalBoxMin, globalBoxMax));
  EXPECT_EQ(false, DomainTools::isInsideDomain(outside, globalBoxMin, globalBoxMax));
}

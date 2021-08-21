/**
 * @file TestDomainTools.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "TestDomainTools.h"

#include "src/domainDecomposition/DomainTools.h"

TEST_F(TestDomainTools, testIsInsideDomain) {
  std::array<double, 3> globalBoxMin = {1.0, 1.0, 1.0};
  std::array<double, 3> globalBoxMax = {10.0, 10.0, 10.0};

  std::array<double, 3> inside = {5.0, 5.0, 5.0};
  std::array<double, 3> onLowerBoundary = {1.0, 2.0, 5.0};
  std::array<double, 3> onUpperBoundary = {2.0, 10.0, 5.0};
  std::array<double, 3> outside = {0.0, 0.0, 0.0};

  EXPECT_EQ(true, DomainTools::isInsideDomain(inside, globalBoxMin, globalBoxMax));
  EXPECT_EQ(true, DomainTools::isInsideDomain(onLowerBoundary, globalBoxMin, globalBoxMax));
  EXPECT_EQ(false, DomainTools::isInsideDomain(onUpperBoundary, globalBoxMin, globalBoxMax));
  EXPECT_EQ(false, DomainTools::isInsideDomain(outside, globalBoxMin, globalBoxMax));
}

TEST_F(TestDomainTools, testBalanceAdjacentDomains) {
  const double workA = 6;
  const double workB = 4;
  const double minA = 2;
  const double maxB = 24;

  EXPECT_EQ(10.8, DomainTools::balanceAdjacentDomains(workA, workB, minA, maxB, 0));
}

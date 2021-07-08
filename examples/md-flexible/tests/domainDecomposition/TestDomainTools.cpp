/**
 * @file TestDomainTools.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "TestDomainTools.h"

#include "autopas/utils/ArrayUtils.h"
#include "src/domainDecomposition/DomainTools.h"

TEST_F(TestDomainTools, testIsInsideDomain) {
  const std::array<double, 3> globalBoxMin = {1.0, 1.0, 1.0};
  const std::array<double, 3> globalBoxMax = {10.0, 10.0, 10.0};

  const std::array<double, 3> inside = {5.0, 5.0, 5.0};
  const std::array<double, 3> onLowerBoundary = {1.0, 2.0, 5.0};
  const std::array<double, 3> onUpperBoundary = {2.0, 10.0, 5.0};
  const std::array<double, 3> outside = {0.0, 0.0, 0.0};

  EXPECT_EQ(true, DomainTools::isInsideDomain(inside, globalBoxMin, globalBoxMax));
  EXPECT_EQ(true, DomainTools::isInsideDomain(onLowerBoundary, globalBoxMin, globalBoxMax));
  EXPECT_EQ(false, DomainTools::isInsideDomain(onUpperBoundary, globalBoxMin, globalBoxMax));
  EXPECT_EQ(false, DomainTools::isInsideDomain(outside, globalBoxMin, globalBoxMax));
}

TEST_F(TestDomainTools, testGenerateDecomposition) {
  const unsigned int subdomainCount = 45;

  std::array<int, 3> decomposition;
  DomainTools::generateDecomposition(subdomainCount, decomposition);

  EXPECT_EQ(3, decomposition[0]);
  EXPECT_EQ(3, decomposition[1]);
  EXPECT_EQ(5, decomposition[2]);
}

TEST_F(TestDomainTools, testConvertIdToIndex) {
  const std::array<int, 3> domainId = {2, 3, 4};
  const std::array<int, 3> decomposition = {3, 4, 5};

  EXPECT_EQ(59, DomainTools::convertIdToIndex(domainId, decomposition));
}

TEST_F(TestDomainTools, testConvertIndexToId) {
  const std::array<int, 3> decomposition = {3, 4, 5};

  std::array<int, 3> domainId = DomainTools::convertIndexToId(59, decomposition);
  EXPECT_EQ(2, domainId[0]);
  EXPECT_EQ(3, domainId[1]);
  EXPECT_EQ(4, domainId[2]);

  domainId = DomainTools::convertIndexToId(34, decomposition);
  EXPECT_EQ(1, domainId[0]);
  EXPECT_EQ(2, domainId[1]);
  EXPECT_EQ(4, domainId[2]);
}

TEST_F(TestDomainTools, testGetAccumulatedTail) {
  const std::array<int, 3> decomposition = {4, 3, 2};
  const std::array<int, 3> indices = {0, 1, 2};

  EXPECT_EQ(6, DomainTools::getAccumulatedTail(indices[0], decomposition));
  EXPECT_EQ(2, DomainTools::getAccumulatedTail(indices[1], decomposition));
  EXPECT_EQ(1, DomainTools::getAccumulatedTail(indices[2], decomposition));
}

TEST_F(TestDomainTools, testGetExtentOfSubdomain) {
  const std::array<int, 3> decomposition = {3, 4, 5};

  std::array<int, 6> extent = DomainTools::getExtentOfSubdomain(59, decomposition);
  EXPECT_EQ(2, extent[0]);
  EXPECT_EQ(3, extent[1]);
  EXPECT_EQ(3, extent[2]);
  EXPECT_EQ(4, extent[3]);
  EXPECT_EQ(4, extent[4]);
  EXPECT_EQ(5, extent[5]);

  extent = DomainTools::getExtentOfSubdomain(34, decomposition);
  EXPECT_EQ(1, extent[0]);
  EXPECT_EQ(2, extent[1]);
  EXPECT_EQ(2, extent[2]);
  EXPECT_EQ(3, extent[3]);
  EXPECT_EQ(4, extent[4]);
  EXPECT_EQ(5, extent[5]);
}

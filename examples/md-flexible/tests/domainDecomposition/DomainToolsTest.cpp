/**
 * @file DomainToolsTest.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "DomainToolsTest.h"

#include "src/domainDecomposition/DomainTools.h"

TEST_F(DomainToolsTest, testIsInsideDomain) {
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

TEST_F(DomainToolsTest, testGenerateDecomposition) {
  const size_t subdomainCount = 27;
  std::array<int, 3> decomposition;

  std::array<bool, 3> subdivideDomains = {true, true, true};
  DomainTools::generateDecomposition(subdomainCount, subdivideDomains);
  EXPECT_EQ(decomposition[0], 3);
  EXPECT_EQ(decomposition[1], 3);
  EXPECT_EQ(decomposition[2], 3);

  subdivideDomains = {false, true, true};
  DomainTools::generateDecomposition(subdomainCount, subdivideDomains);
  EXPECT_EQ(decomposition[0], 1);
  EXPECT_EQ(decomposition[1], 9);
  EXPECT_EQ(decomposition[2], 3);

  subdivideDomains = {true, false, true};
  DomainTools::generateDecomposition(subdomainCount, subdivideDomains);
  EXPECT_EQ(decomposition[0], 9);
  EXPECT_EQ(decomposition[1], 1);
  EXPECT_EQ(decomposition[2], 3);

  subdivideDomains = {true, true, false};
  DomainTools::generateDecomposition(subdomainCount, subdivideDomains);
  EXPECT_EQ(decomposition[0], 9);
  EXPECT_EQ(decomposition[1], 3);
  EXPECT_EQ(decomposition[2], 1);

  subdivideDomains = {true, false, false};
  DomainTools::generateDecomposition(subdomainCount, subdivideDomains);
  EXPECT_EQ(decomposition[0], 27);
  EXPECT_EQ(decomposition[1], 1);
  EXPECT_EQ(decomposition[2], 1);

  subdivideDomains = {false, true, false};
  DomainTools::generateDecomposition(subdomainCount, subdivideDomains);
  EXPECT_EQ(decomposition[0], 1);
  EXPECT_EQ(decomposition[1], 27);
  EXPECT_EQ(decomposition[2], 1);

  subdivideDomains = {false, false, true};
  DomainTools::generateDecomposition(subdomainCount, subdivideDomains);
  EXPECT_EQ(decomposition[0], 1);
  EXPECT_EQ(decomposition[1], 1);
  EXPECT_EQ(decomposition[2], 27);
}

/**
 * @file TestRegularGrid.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "TestRegularGrid.h"

#include "src/domainDecomposition/RegularGrid.h"

TEST_F(TestRegularGrid, testGetGlobalDomain) {
  std::vector<double> globalBoxMin = {1.0, 1.0, 1.0};
  std::vector<double> globalBoxMax = {10.0, 10.0, 10.0};

  char** argv;
  
  RegularGrid domainDecomposition(0, argv, 3, globalBoxMin, globalBoxMax);

  EXPECT_EQ(globalBoxMin, domainDecomposition.getGlobalBoxMin());
  EXPECT_EQ(globalBoxMax, domainDecomposition.getGlobalBoxMax());
}

TEST_F(TestRegularGrid, testGetLocalDomain) {
  std::vector<double> globalBoxMin = {1.0, 1.0, 1.0};
  std::vector<double> globalBoxMax = {10.0, 10.0, 10.0};

  char** argv;

  RegularGrid domainDecomposition(0, argv, 3, globalBoxMin, globalBoxMax);

  std::vector<double> globalBoxExtend = autopas::utils::ArrayMath::sub(globalBoxMax, globalBoxMin);
  std::vector<int> decomposition = domainDecomposition.getDecomposition();

  // @todo: finish this test
  //std::vector expectedlocalBoxExtend = 

  double boxLengthX = globalBoxExtend[0] / decomposition[0];
  double boxLengthY = globalBoxExtend[1] / decomposition[1];
  double boxLengthZ = globalBoxExtend[2] / decomposition[2];
  
  EXPECT_EQ(globalBoxMin, domainDecomposition.getLocalBoxMin());
  EXPECT_EQ(globalBoxMax, domainDecomposition.getLocalBoxMax());
}

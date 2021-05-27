/**
 * @file TestRegularGrid.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "TestRegularGrid.h"

#include "mpi.h"
#include "src/domainDecomposition/DomainTools.h"
#include "src/domainDecomposition/RegularGrid.h"

extern int extArgc;
extern char** extArgv;

namespace {
  std::vector<double> sub(std::vector<double> a, std::vector<double> b) {
    std::vector<double> difference = a;
    if (a.size() == b.size()){
      for (int i = 0; i < a.size(); ++i){
        difference[i] = a[i] - b[i];
      }
    }
    return difference;
  }

  std::vector<double> div(std::vector<double> a, std::vector<int> b) {
    std::vector<double> result = a;
    if (a.size() == b.size()){
      for (int i = 0; i < a.size(); ++i){
        result[i] = a[i] / (double)b[i];
      }
    }
    return result;
  }
}

TEST_F(TestRegularGrid, testGetLocalDomain) {
  std::vector<double> globalBoxMin = {1.0, 1.0, 1.0};
  std::vector<double> globalBoxMax = {10.0, 10.0, 10.0};

  char** argv;

  std::cout << "Arguments: " << extArgc << ", " << extArgv[1] << std::endl;
  RegularGrid domainDecomposition(extArgc, extArgv, 3, globalBoxMin, globalBoxMax);

  std::vector<double> globalBoxExtend = sub(globalBoxMax, globalBoxMin);

  int numberOfProcesses;
  MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

  std::vector<int> decomposition;
  DomainTools::generateDecomposition(numberOfProcesses, 3, decomposition);

  std::vector<double> expectedLocalBoxExtend = div(globalBoxExtend, decomposition);

  std::vector<double> resultingLocalBoxExtend =
    sub(domainDecomposition.getLocalBoxMax(), domainDecomposition.getLocalBoxMin());

  EXPECT_NEAR(expectedLocalBoxExtend[0], resultingLocalBoxExtend[0], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[1], resultingLocalBoxExtend[1], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[2], resultingLocalBoxExtend[2], 1e-10);
}

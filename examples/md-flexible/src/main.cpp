/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "Simulation.h"
#include "autopas/utils/WrapMPI.h"

/**
 * The main function for md-flexible.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  autopas::AutoPas_MPI_Init(&argc, &argv);

  MDFlexConfig configuration(argc, argv);

  const std::vector<double> boxMin(configuration.boxMin.value.begin(), configuration.boxMin.value.end());
  const std::vector<double> boxMax(configuration.boxMax.value.begin(), configuration.boxMax.value.end());

  RegularGridDecomposition domainDecomposition(3, boxMin, boxMax, configuration.cutoff.value,
                                               configuration.verletSkinRadius.value);

  const std::vector<double> localBoxMin = domainDecomposition.getLocalBoxMin();
  const std::vector<double> localBoxMax = domainDecomposition.getLocalBoxMax();

  for (int i = 0; i < localBoxMin.size(); ++i) {
    configuration.boxMin.value[i] = localBoxMin[i];
    configuration.boxMax.value[i] = localBoxMax[i];
  }

  Simulation simulation(configuration, domainDecomposition, argc, argv);
  simulation.run();

  if (domainDecomposition.getDomainIndex() == 0) {
    std::cout << std::endl << "Using " << autopas::autopas_get_max_threads() << " Threads" << std::endl;
  }

  autopas::AutoPas_MPI_Finalize();
  return EXIT_SUCCESS;
}

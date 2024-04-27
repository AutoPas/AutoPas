#include <fstream>
#include <iostream>
#include <chrono>

#include "SimulationExtForTuning.h"
#include "autopas/utils/WrapMPI.h"

using namespace std::chrono;

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<ParticleType>;

/**
 * The main function for md-flexible-tuner.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  autopas::AutoPas_MPI_Init(&argc, &argv);
  {
    MDFlexConfig configuration(argc, argv);

    auto domainDecomposition = std::make_shared<RegularGridDecomposition>(configuration);

    if (not configuration.checkpointfile.value.empty()) {
      configuration.flushParticles();
      configuration.loadParticlesFromCheckpoint(domainDecomposition->getDomainIndex(),
                                                domainDecomposition->getSubdomainCount());
    }

    // print start configuration and parallelization info
    if (domainDecomposition->getDomainIndex() == 0) {
      std::cout << configuration.to_string() << std::endl;
      std::cout << std::endl << "Using " << autopas::autopas_get_max_threads() << " Threads" << std::endl;
#if defined(AUTOPAS_INCLUDE_MPI)
      std::cout << "MPI is running with " << domainDecomposition->getNumberOfSubdomains() << " ranks." << std::endl;
#else
      std::cout << "MPI is disabled." << std::endl;
#endif
    }

  SimulationExtForTuning simulation(configuration, domainDecomposition);
  std::ofstream out("values.txt",std::ios_base::out | std::ios_base::trunc);
  if (out.is_open()) {
    //out << "Nr. Particles | Difference" << '\n';
    out << "Nr. Particles | SortThres 0 | SortThres Max | Difference" << '\n';
    int optimalSortingThreshold = 0;
    int numParticles = 0;
    int runs = 3;
    for (int run = 1; run <= runs; run++) {
      std::cout << "-------------------" << std::endl;
      int priorDifference = 0;
      int currentDifference = 0;

    //while (priorDifference >= 0) {
    for (numParticles = 2; numParticles <= 20 /*configuration.getParticles().size();*/; numParticles+=1) {
      //std::cout << "Number of Particles: " << numParticles << std::endl;
      out << "       " << numParticles;
      int repeat = 100000;

    //SortingThreshold = 0
      milliseconds msStart = duration_cast< milliseconds >(
          system_clock::now().time_since_epoch()
      );
      configuration.sortingThreshold.value = 0;
      for (int ii = 0; ii < repeat; ii++) {
          simulation.processCellPair(numParticles, /*sortingThreshold =*/ 0);
      }
      milliseconds msEnd = duration_cast< milliseconds >(
      system_clock::now().time_since_epoch()
      );
      milliseconds msDifferenceThreshold0 = msEnd - msStart;
      // std::cout  << "Number of Particles: " << numParticles << " SortingThreshold: 0. The Simulation took: " << msDifferenceThreshold0.count() << " milliseconds." << std::endl;
      out << "           " << msDifferenceThreshold0.count();

    //SortingThreshold = Max
      msStart = duration_cast< milliseconds >(
          system_clock::now().time_since_epoch()
      );
      configuration.sortingThreshold.value = configuration.getParticles().size();
      for (int ii = 0; ii < repeat; ii++) {
        simulation.processCellPair(numParticles, /*sortingThreshold =*/ INT_MAX);
      }
      msEnd = duration_cast< milliseconds >(
      system_clock::now().time_since_epoch()
      );
      milliseconds msDifferenceThresholdMax = msEnd - msStart;
      //std::cout  << "Number of Particles: " << numParticles << " SortingThreshold: Max. The Simulation took: " << msDifferenceThresholdMax.count() << " milliseconds." << std::endl;
      out << "            " << msDifferenceThresholdMax.count();

      currentDifference = msDifferenceThreshold0.count() - msDifferenceThresholdMax.count();
      std::cout   << "Number of Particles: " << numParticles << " diff =" << currentDifference << std::endl;
      //out << "       " << numParticles << "           " << currentDifference << '\n';
      out << "              " << currentDifference << '\n';

      if (currentDifference < 0) {
        if (abs(currentDifference) < priorDifference) {
          //std::cout << std::endl << "SortingThreshold: " << numParticles << std::endl << std::endl;
          //out << '\n' << "SortingThreshold: " << numParticles << '\n';
          optimalSortingThreshold += numParticles;
        }
        else {
          //std::cout << std::endl << "SortingThreshold: " << numParticles-1 << std::endl << std::endl;
          //out << '\n' << "SortingThreshold: " << numParticles-1 << '\n';
          optimalSortingThreshold += numParticles-1;
        }
        break;
      }
      priorDifference = currentDifference;
    }
    out << '\n';
    }
    if (optimalSortingThreshold == 0) {
      optimalSortingThreshold = numParticles * runs;
    }
    optimalSortingThreshold = (optimalSortingThreshold + runs/2) / runs;
    std::cout << std::endl << "Your optimal SortingThreshold is " << optimalSortingThreshold << std::endl << std::endl;
    out << '\n' << "Optimal SortingThreshold: " << optimalSortingThreshold << '\n';

  }
  out.close();
  }
  autopas::AutoPas_MPI_Finalize();
  return EXIT_SUCCESS;
}

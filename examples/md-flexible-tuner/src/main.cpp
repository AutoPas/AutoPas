#include <fstream>
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

    //Calculate optimal SortingThreshold for ProcessCellPair:

    int optimalSortingThreshold = 0;
    int numParticles = 20;           //Or int numParticles = configuration.getParticles().size();
    int runs = 3;                    //Run the calculation several times
    int repeat = 100000;             //Repeat calling ProcessCell/ProcessCellPair per run

    for (int run = 1; run <= runs; run++) {
      std::cout << "-------------------" << std::endl;
      int priorDifference = 0;
      int currentDifference = 0;

      for (numParticles = 2; numParticles <= 20 /*configuration.getParticles().size();*/; numParticles+=1) {

        //SortingThreshold = 0
        milliseconds msStart = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
       configuration.sortingThreshold.value = 0;
        for (int ii = 0; ii < repeat; ii++) {
            simulation.processCellPair(numParticles, /*sortingThreshold =*/ 0);
        }
        milliseconds msEnd = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
        milliseconds msDifferenceThreshold0 = msEnd - msStart;

        //SortingThreshold = Max
        msStart = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
        configuration.sortingThreshold.value = configuration.getParticles().size();
        for (int ii = 0; ii < repeat; ii++) {
          simulation.processCellPair(numParticles, /*sortingThreshold =*/ INT_MAX);
        }
        msEnd = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
         milliseconds msDifferenceThresholdMax = msEnd - msStart;

        currentDifference = msDifferenceThreshold0.count() - msDifferenceThresholdMax.count();
        std::cout   << "Number of Particles: " << numParticles << " diff =" << currentDifference << std::endl;
 
        //When the optimal SortingThreshold lies between two numbers, calculate which one it is closer to
        if (currentDifference < 0) {
          if (abs(currentDifference) < priorDifference) {
            optimalSortingThreshold += numParticles;
          }
          else {
            optimalSortingThreshold += numParticles-1;
          }
          break;
        }
        priorDifference = currentDifference;
      }
      if (numParticles == 20) {
        optimalSortingThreshold += numParticles;
      }
    }
    optimalSortingThreshold = (optimalSortingThreshold + runs/2) / runs;
    std::cout << std::endl << "Optimal SortingThreshold for ProcessCellPair: " << optimalSortingThreshold << std::endl;


    //Calculate optimal SortingThreshold for ProcessCell:

    optimalSortingThreshold = 0;

    for (int run = 1; run <= runs; run++) {
      std::cout << "-------------------" << std::endl;
      int priorDifference = 0;
      int currentDifference = 0;

      for (numParticles = 2; numParticles <= 20 /*configuration.getParticles().size();*/; numParticles+=1) {

        //SortingThreshold = 0
        milliseconds msStart = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
       configuration.sortingThreshold.value = 0;
        for (int ii = 0; ii < repeat; ii++) {
            simulation.processCell(numParticles, /*sortingThreshold =*/ 0);
        }
        milliseconds msEnd = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
        milliseconds msDifferenceThreshold0 = msEnd - msStart;

        //SortingThreshold = Max
        msStart = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
        configuration.sortingThreshold.value = configuration.getParticles().size();
        for (int ii = 0; ii < repeat; ii++) {
          simulation.processCell(numParticles, /*sortingThreshold =*/ INT_MAX);
        }
        msEnd = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
         milliseconds msDifferenceThresholdMax = msEnd - msStart;

        currentDifference = msDifferenceThreshold0.count() - msDifferenceThresholdMax.count();
        std::cout   << "Number of Particles: " << numParticles << ", diff =" << currentDifference << std::endl;

        //When the optimal SortingThreshold lies between two numbers, calculate which one it is closer to
        if (currentDifference < 0) {
          if (abs(currentDifference) < priorDifference) {
            optimalSortingThreshold += numParticles;
          }
          else {
            optimalSortingThreshold += numParticles-1;
          }
          break;
        }
        priorDifference = currentDifference;
      }
      if (numParticles == 20) {
        optimalSortingThreshold += numParticles;
      }
    }
    optimalSortingThreshold = (optimalSortingThreshold + runs/2) / runs;
    std::cout << std::endl << "Optimal SortingThreshold for ProcessCell: " << optimalSortingThreshold << std::endl;
  }

  autopas::AutoPas_MPI_Finalize();
  return EXIT_SUCCESS;
}

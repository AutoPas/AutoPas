/**
 *
 * @file main.cpp
 * @date 01.07.19
 * @author M. Geitner
 */


#include <chrono>
#include <fstream>
#include <iostream>
#include "autopas/AutoPas.h"


using namespace std;
using namespace autopas;

int main(){



  std::chrono::high_resolution_clock::time_point startTotal, stopTotal;

  startTotal = std::chrono::high_resolution_clock::now();
  string logFileName = "log_kokkos.txt";

  // select either std::out or a logfile for autopas log output.
  // This does not affect md-flex output.
  std::ofstream logFile;
  std::streambuf *streamBuf;
  if (logFileName.empty()) {
    streamBuf = std::cout.rdbuf();
  } else {
    logFile.open(logFileName);
    streamBuf = logFile.rdbuf();
  }
  std::ostream outputStream(streamBuf);
  // Initialization
  autopas::AutoPas<KokkosParticle, FullParticleCell<KokkosParticle>> autopas(outputStream);
  //autopas::Logger::get()->set_level(logLevel);
  auto functor = KokkosLJFunctor<KokkosParticle, FullParticleCell<KokkosParticle>>();

  double cutoff = 1.0;
  double sigma = 1.0;
  double epsilon = 1.0;

  autopas.setCutoff(cutoff);
  /*autopas.setVerletSkin(verletSkinRadius);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setTuningInterval(tuningInterval);
  autopas.setTuningStrategyOption(tuningStrategy);
  autopas.setNumSamples(tuningSamples);
  autopas.setSelectorStrategy(selectorStrategy);*/
  std::set<autopas::ContainerOption> containerOptions{};
  containerOptions.insert(autopas::ContainerOption::directSum);
  std::set<autopas::TraversalOption> traversalOptions{};
  traversalOptions.insert(autopas::TraversalOption::kokkosDirectSumTraversal);
  std::set<autopas::DataLayoutOption > dataLayoutOptions{};
  dataLayoutOptions.insert(autopas::DataLayoutOption::aos);
  std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled};


  autopas.setAllowedContainers(containerOptions);
  autopas.setAllowedTraversals(traversalOptions);
  autopas.setAllowedDataLayouts(dataLayoutOptions);
  autopas.setAllowedNewton3Options(newton3Options);
  //autopas.setAllowedCellSizeFactors(cellSizeFactors);

  //set boxMin and boxMax
  double boxLength = 10;
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();//init autopas before particles are added

  std::array<KokkosParticle, 5> arrParticles{};
  for (int i = 0; i < 5; ++i) {
    arrParticles[i] = KokkosParticle({0.1 * i , 0.2 * i, 0.3 * i},{0.0, 0.0, 0.0}, i);
  }
  for (int i = 0; i < 5; ++i) {
    autopas.addParticle(arrParticles[i]);
  }
  autopas.iteratePairwise(&functor);//iterate
  for (int i = 0; i < 5; ++i) {
    std::cout <<arrParticles[i].toString() << "\n";
  }
  return EXIT_SUCCESS;
}
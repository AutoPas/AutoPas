/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include <mpi.h>

#include "MDFlexSimulation.h"
#include "src/domainDecomposition/RegularGrid.h"

class MDFlexMPI : public MDFlexSimulation {
 public:
  MDFlexMPI(int dimensionCount, int argc, char **argv);
  ~MDFlexMPI() = default;

  void run() override;
  void initializeDomainDecomposition(int &dimensionCount) override;

  DomainDecomposition *getDomainDecomposition() override {
    return static_cast<DomainDecomposition *>(&(*_domainDecomposition));
  }

 private:
  std::shared_ptr<RegularGrid> _domainDecomposition;

  void updateParticles();
  void executeSuperstep(const int iterationsPerSuperstep);

  void sendParticles(std::vector<ParticleType> &particles, int &receiver);
  void receiveParticles(std::vector<ParticleType> &receivedParticles, int &source);
};

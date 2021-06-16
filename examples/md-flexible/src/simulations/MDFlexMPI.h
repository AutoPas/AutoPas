/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"
#include "src/domainDecomposition/RegularGrid.h"

class MDFlexMPI : public MDFlexSimulation {
 public:
  /**
   * Constructor.
   * @param dimensionCount The number of dimensions to use for the domain decomposition.
   * @param argc The argument count provided to the main function.
   * @param argv The argument vector provided to the main function.
   */
  MDFlexMPI(int dimensionCount, int argc, char **argv);

  /**
   * Destructor
   */
  ~MDFlexMPI() = default;

  /**
   * Runs the simulation.
   */
  void run() override;

  /**
   * Initializes the domain decomposition for this simulation.
   */
  void initializeDomainDecomposition(int &dimensionCount) override;

  /**
   * Returns the domain decomposition for this simulation.
   */
  DomainDecomposition *getDomainDecomposition() override {
    return static_cast<DomainDecomposition *>(&(*_domainDecomposition));
  }

 private:
  /**
   * This simulation's domain decomposition.
   */
  std::shared_ptr<RegularGrid> _domainDecomposition;

  /**
   * Updates the particles in the local AutoPas container.
   */
  void updateParticles();

  /**
   * Executes a superstep of the simulation.
   */
  void executeSuperstep(const int iterationsPerSuperstep);

  /**
   * Sends particles of type ParticleType to a specific receiver.
   * @param particles The particles to be sent to the receiver.
   * @param receiver The recipient of the particles.
   */
  void sendParticles(std::vector<ParticleType> &particles, int &receiver);

  /**
   * Receives particels of type ParticleType which have been send by a specific sender.
   * @param receivedParticels The container where the received particles will be stored.
   * @param source The sender of the particles.
   */
  void receiveParticles(std::vector<ParticleType> &receivedParticles, int &source);
};

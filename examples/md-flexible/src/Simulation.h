/**
 * @file Simulation.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>

#include "autopas/AutoPas.h"
#include "src/ParallelVtkWriter.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/domainDecomposition/DomainDecomposition.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

/**
 * Handles minimal initialization requriements for MD-Flexible simulations.
 * Derivce this class to create custom simulations.
 */
class Simulation {
 public:
  /**
   * Initializes the simulation on a domain according to the arguments passed to the main function.
   * @param dimensionCount The number of dimensions in the simulation domain.
   * @param argc The number of arguments passed in argv.
   * @param argv The arguments passed to the program.
   */
  Simulation(const MDFlexConfig &configuration, RegularGridDecomposition &domainDecomposition);

  /**
   * Destructor.
   */
  ~Simulation() = default;

  /**
   * Runs the simulation
   */
  void run();

 protected:
  /**
   * Stores the configuration used for the simulation.
   * The configuration is defined by the .yaml file passed to the application  with the '--yaml-file' argument.
   */
  MDFlexConfig _configuration;

  /**
   * The the nodes' autopas container used for simulation.
   * This member will not be initialized by the constructor and therefore has to be initialized by the deriving class.
   */
  std::shared_ptr<autopas::AutoPas<ParticleType>> _autoPasContainer;

  /**
   * Shared pointer to the logfile.
   */
  std::shared_ptr<std::ofstream> _logFile;

  /**
   * Pointer to the output stream.
   */
  std::ostream *_outputStream;

  /**
   * Number of completed iterations. Aka. number of current iteration.
   */
  size_t _iteration = 0;

  /**
   * Counts completed iterations that were used for tuning
   */
  size_t _numTuningIterations = 0;

  /**
   * Counts completed tuning phases.
   */
  size_t _numTuningPhasesCompleted = 0;

  /**
   * Indicator if the previous iteration was used for tuning.
   */
  bool _previousIterationWasTuningIteration = false;

  /**
   * Precision of floating point numbers printed.
   */
  constexpr static auto _floatStringPrecision = 3;

  /**
   * Homogeneity of the scenario, calculated by the standard deviation of the density.
   */
  double _homogeneity = 0;

  /**
   * Struct containing all timers used for the simulation.
   */
  struct Timers {
    autopas::utils::Timer positionUpdate;
    autopas::utils::Timer forceUpdateTotal;
    autopas::utils::Timer forceUpdatePairwise;
    autopas::utils::Timer forceUpdateGlobal;
    autopas::utils::Timer forceUpdateTuning;
    autopas::utils::Timer forceUpdateNonTuning;
    autopas::utils::Timer velocityUpdate;
    autopas::utils::Timer simulate;
    autopas::utils::Timer vtk;
    autopas::utils::Timer initialization;
    autopas::utils::Timer total;
    autopas::utils::Timer thermostat;
    autopas::utils::Timer boundaries;
  } _timers;

  /**
   * Parallel VTK file writer.
   */
  std::shared_ptr<ParallelVtkWriter> _vtkWriter;

  /**
   * Defines, if vtk files should be created or not.
   */
  bool _createVtkFiles;

  /**
   * Initializes the simulation.
   * This function needs to be called in the constructor of the deriving class, because initializeDomainDecomposition
   * can not be called by the constructor of Simulation, because it is a pure virtual function.
   */
  void initialize(int dimensionCount, int argc, char **argv);

  /**
   * Executes a superstep of the simulation.
   */
  void executeSuperstep(const int iterationsPerSuperstep);

  /**
   * Checks if there are any iterations left to compute.
   */
  bool needsMoreIterations();

  /**
   * Estimates the number of tuning iterations which ocurred during the simulation so far.
   */
  std::tuple<size_t, bool> estimateNumberOfIterations() const;

  /**
   * Prints a progress bar to the terminal.
   */
  void printProgress(size_t iterationProgress, size_t maxIterations, bool maxIsPrecise);

  /**
   * Writes the current simulation state to a vtk-file.
   */
  void writeVTKFile();

  /**
   * Returns an 'mpi_rank_<rank>_', where <rank> is the rank of the current MPI process.
   */
  std::string getMPISuffix();

  /**
   * Turns the timers into a human readable string.
   */
  std::string timerToString(const std::string &name, long timeNS, size_t numberWidth, long maxTime);

  /**
   * Updates the position of particles in the local AutoPas container.
   */
  void updatePositions();

  /**
   * Updates the forces of particles in the local AutoPas container.
   */
  void updateForces();

  /**
   * Updates the velocities of particles in the local AutoPas container.
   */
  void updateVelocities();
  /**
   * Updates the thermostat of for the local domain.
   * @todo The thermostat shoud act globally and therefore needs to be communicated to all processes.
   */
  void updateThermostat();

 private:
  /**
   * This simulation's domain decomposition.
   */
  RegularGridDecomposition _domainDecomposition;

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

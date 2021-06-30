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

#include "autopas/AutoPasDecl.h"
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
   * @param configuration: The configuration of this simulation.
   * @param domainDecomposition: The domain decomposition used for this simulation
   */
  Simulation(const MDFlexConfig &configuration, RegularGridDecomposition &domainDecomposition);

  /**
   * Destructor.
   */
  ~Simulation();

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
    /**
     * Records the time used for the position updates of all particles.
     */
    autopas::utils::Timer positionUpdate;

    /**
     * Records the time used for the total force update of all particles.
     */
    autopas::utils::Timer forceUpdateTotal;

    /**
     * Records the time used for the pairwise force update of all particles.
     */
    autopas::utils::Timer forceUpdatePairwise;

    /**
     * Records the time used for the global force update of all particles.
     */
    autopas::utils::Timer forceUpdateGlobal;

    /**
     * Records the time used for the force update of all particles during the tuning iterations.
     */
    autopas::utils::Timer forceUpdateTuning;

    /**
     * Records the time used for force updates of all particles during the non tuning iterations.
     */
    autopas::utils::Timer forceUpdateNonTuning;

    /**
     * Records the time used for the velocity updates of all particles.
     */
    autopas::utils::Timer velocityUpdate;

    /**
     * Records the time used for actively simulating the provided scenario.
     * This excludes initialization time, among others.
     */
    autopas::utils::Timer simulate;

    /**
     * Records the time used for the creation of the timestep records.
     */
    autopas::utils::Timer vtk;

    /**
     * Records the time used for the initialization of the simulation.
     */
    autopas::utils::Timer initialization;

    /**
     * Records the total time required for the simulation.
     */
    autopas::utils::Timer total;

    /**
     * Records the time required for the thermostat updates.
     */
    autopas::utils::Timer thermostat;
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
   * Executes a sequence of  supersteps for the simulation.
   * @param iterations: The number of iterations which will be simulated during the excution of this function.
   */
  void executeSupersteps(const int iterations);

  /**
   * Estimates the number of tuning iterations which ocurred during the simulation so far.
   * @return an estimation of the number of tuning iterations which occured so far.
   */
  std::tuple<size_t, bool> estimateNumberOfIterations() const;

  /**
   * Prints a progress bar to the terminal.
   * @param iterationProgress: the number of already computed iterations.
   * @param maxIterations: the number of maximum iterations.
   * @param maxIsPrecise: Decides if the "~" symbol will be printed before the max iterations.
   */
  void printProgress(size_t iterationProgress, size_t maxIterations, bool maxIsPrecise);

  /**
   * Turns the timers into a human readable string.
   * @param name: The timer's name.
   * @param timeNS: The time in nano seconds.
   * @param numberWidth: The precision of the printed number.
   * @param maxTime: The simulation's total execution time.
   * @return All information of the timer in a human readable string.
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

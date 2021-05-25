/**
 * @file MDFlexSimulation.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <tuple>

#include "../TypeDefinitions.h"
#include "../configuration/MDFlexConfig.h"
#include "../domainDecomposition/DomainDecomposition.h"
#include "autopas/AutoPas.h"

/**
 * Handles minimal initialization requriements for MD-Flexible simulations.
 * Derivce this class to create custom simulations.
 */
class MDFlexSimulation {
 public:
  /**
   * Runs the simulation
   */
  virtual void run() = 0;
  
  /**
   * Initializes the domain decomposition for this simulation.
   */
  virtual void initializeDomainDecomposition(int &dimensionCount) = 0;

  /**
   * Returns the domain decomposition for this simulation.
   */
  virtual DomainDecomposition *getDomainDecomposition() = 0;

 protected:
  /**
   * Initializes the simulation on a domain according to the arguments passed to the main function.
   * @param dimensionCount The number of dimensions in the simulation domain.
   * @param argc The number of arguments passed in argv.
   * @param argv The arguments passed to the program.
   */
  MDFlexSimulation() = default;

  /**
   * Destructor.
   */
  virtual ~MDFlexSimulation();

  /**
   * Stores the argument count passed to the constructor for later reuse.
   */
  int _argc;

  /**
   * Stores the arguments passed to the constructor for later reuse.
   */
  char **_argv;

  /**
   * Stores the configuration used for the simulation.
   * The configuration is defined by the .yaml file passed to the application  with the '--yaml-file' argument.
   */
  std::shared_ptr<MDFlexConfig> _configuration;

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
   * Initializes the simulation.
   * Call this function in the constructor of derived classes.
   */
  void initialize(int dimensionCount, int argc, char **argv);

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

 private:
  /**
   * Initializes the local AutoPas container.
   */
  void initializeAutoPasContainer();
};

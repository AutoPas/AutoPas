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

#include "GlobalVariableLogger.h"
#include "GlobalVariableLogger.h"
#include "TimeDiscretization.h"
#include "autopas/AutoPasDecl.h"
#include "src/ParallelVtkWriter.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/domainDecomposition/DomainDecomposition.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

/**
 * Handles minimal initialization requirements for MD-Flexible simulations.
 * Derive this class to create custom simulations.
 */
class Simulation {
 public:
  /**
   * Initializes the simulation on a domain according to the arguments passed to the main function.
   * @param configuration: The configuration of this simulation.
   * @param domainDecomposition: The domain decomposition used for this simulation
   */
  Simulation(const MDFlexConfig &configuration, std::shared_ptr<RegularGridDecomposition> &domainDecomposition);

  /**
   * Destructor.
   */
  ~Simulation() = default;

  /**
   * Runs the simulation, implementing velocity verlet.
   *
   * If md-flexible is compiled for multi-site molecules, rotational integration is done with an implementation of the
   * the quaternion approach (A) as described in Rozmanov, 2010, Robust rotational-velocity-Verlet integration methods.
   *
   * @note: An initial force calculation should actually take place before the first position update, according to
   * Störmer Verlet the literature (Griebel et al., Numerical simulation in molecular dynamics: numerics, algorithms,
   * parallelization, applications, 2007). In this implementation, the first position update is calculated with zero
   * forces, if none are given initially, and only the velocity is used. This means that, essentially, the molecules do
   * not start the simulation at the given positions but are already moved by one timestep in the direction of their
   * initial velocity. If you want to observe the exact movement of particles, you should bear this in mind. The initial
   * force calculation is omitted here for the following reasons:
   * 1. The primary focus for md-flexible is to be a test application and demonstrator of how AutoPas can be used. Not
   * to be a perfect simulator.
   * 2. md-flexible is often used during development for performance measurements or memory analysis. An additional
   * force calculation step would unnecessarily complicate the simulation loop and distort the statistics.
   * 3. We also do not need an initial force calculation or the oldForce for checkpoint loading. Reason: In the last
   * simulation step before saving the checkpoint, F was calculated, and a velocity update was performed with this F. In
   * the first iteration after loading the checkpoint file, the position update is calculated with this F (loaded from
   * the checkpoint file) from the last simulation step, then F_old is set to F and the new force is calculated. This
   * automatically continues the Störmer Verlet algorithm seamlessly.
   */
  void run();

  /**
   * Finalizes the simulation.
   * Stops remaining timers and logs the result of all the timers.
   * This needs to be called before MPI_Finalize if MPI is enabled.
   */
  void finalize();

 protected:
  /**
   * Stores the configuration used for the simulation.
   * The configuration is defined by the .yaml file passed to the application  with the '--yaml-file' argument.
   */
  MDFlexConfig _configuration;

  /**
   * The the nodes' AutoPas container used for simulation.
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
   * Variable to store total potential energy for the current iteration from both 2B and 3B interactions.
   */
  double _totalPotentialEnergy{};

  /**
   * Variable to store total virial sum for the current iteration from both 2B and 3B interactions.
   */
  double _totalVirialSum{};

  /**
   * Logger for global variables calculated by the Functor, eg. Virial, Potential Energy, etc.
   */
  std::unique_ptr<GlobalVariableLogger> _globalLogger;

  /**
   * Number of completed iterations. Aka. number of current iteration.
   * The first iteration has number 0.
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
   * Indicator if the current iteration was a tuning iteration.
   */
  bool _currentIterationIsTuningIteration = false;

  /**
   * Indicator if the simulation is paused due to the PauseDuringTuning option.
   */
  bool _simulationIsPaused = false;

  /**
   * Indicator if the previous iteration was used for tuning.
   */
  bool _previousIterationWasTuningIteration = false;

  /**
   * Precision of floating point numbers printed.
   */
  constexpr static auto _floatStringPrecision = 3;

  /**
   * Struct containing all timers used for the simulation.
   */
  struct Timers {
    /**
     * Records the time used for the position updates of all particles.
     */
    autopas::utils::Timer positionUpdate;

    /**
     * Records the time used for the quaternion updates of all particles.
     */
    autopas::utils::Timer quaternionUpdate;

    /**
     * Records the time used for the total force update of all particles.
     */
    autopas::utils::Timer forceUpdateTotal;

    /**
     * Records the time used for the pairwise force update of all particles.
     */
    autopas::utils::Timer forceUpdatePairwise;

    /**
     * Records the time used for the triwise force update of all particles.
     */
    autopas::utils::Timer forceUpdateTriwise;

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
     * Records the time used for the angular velocity updates of all particles.
     */
    autopas::utils::Timer angularVelocityUpdate;

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

    /**
     * Records the time required to exchange the halo particles.
     */
    autopas::utils::Timer haloParticleExchange;

    /**
     * Records the time required to reflect particles.
     */
    autopas::utils::Timer reflectParticlesAtBoundaries;

    /**
     * Records the time required to exchange migrating particles.
     */
    autopas::utils::Timer migratingParticleExchange;

    /**
     * Records the time required for load balancing.
     */
    autopas::utils::Timer loadBalancing;

    /**
     * Used for the diffuse load balancing as the metric to determine the imbalance.
     */
    autopas::utils::Timer computationalLoad;

    /**
     * Records the time required for the update of the AutoPas container.
     */
    autopas::utils::Timer updateContainer;
  };

  /**
   * Sensor for energy measurement
   */
  autopas::utils::EnergySensor _totalEnergySensor;

  /**
   * The timers used during the simulation.
   */
  struct Timers _timers;

  /**
   * Parallel VTK file writer.
   */
  std::shared_ptr<ParallelVtkWriter> _vtkWriter;

  /**
   * Defines, if vtk files should be created or not.
   */
  bool _createVtkFiles;

 private:
  /**
   * Load particles from this object's config into this object's AutoPas container.
   *
   * @note This also clears all particles from the config file!
   */
  void loadParticles();

  /**
   * Returns the number of expected maximum number of iterations of the Simulation.
   * This is exact if the number of iterations was specified, or an estimate based on the number of tuning iterations if
   * the number of tuning phases is specified.
   * @return <size_t, bool> Max number of iterations, bool whether it is an estimate.
   */
  [[nodiscard]] std::tuple<size_t, bool> estimateNumberOfIterations() const;

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
   * @param numberWidth: The minimal field width of the printed number.
   * @param maxTime: The simulation's total execution time.
   * @return All information of the timer in a human readable string.
   */
  [[nodiscard]] static std::string timerToString(const std::string &name, long timeNS, int numberWidth = 0,
                                                 long maxTime = 0ul);

  /**
   * Updates the position of particles in the local AutoPas container. In addition, the oldForce is set to the value of
   * the current forces and the force buffers of the particles are reset to the global force.
   */
  void updatePositionsAndResetForces();

  /**
   * Update the quaternion orientation of the particles in the local AutoPas container.
   */
  void updateQuaternions();

  /**
   * Updates the forces of particles in the local AutoPas container. Includes torque updates (if an appropriate functor
   * is used).
   */
  void updateInteractionForces();

  /**
   * Updates the velocities of particles in the local AutoPas container.
   */
  void updateVelocities();

  /**
   * Updates the angular velocities of the particles in the local AutoPas container.
   */
  void updateAngularVelocities();

  /**
   * Updates the thermostat of for the local domain.
   * @todo The thermostat should act globally and therefore needs to be communicated to all processes.
   */
  void updateThermostat();

  /**
   * This simulation's domain decomposition.
   */
  std::shared_ptr<RegularGridDecomposition> _domainDecomposition;
  /**
   * If MPI is enabled, accumulates the times of all ranks on rank 0.
   * Otherwise, this function does nothing.
   * @param time: the time to accumulate.
   * @return the accumulated time of all ranks.
   */
  [[nodiscard]] static long accumulateTime(const long &time);

  /**
   * Handles the pausing of the simulation and updates the _simulationIsPaused flag.
   */
  void updateSimulationPauseState();

  /**
   * Logs the number of total/owned/halo particles in the simulation, aswell as the standard deviation of Homogeneity.
   */
  void logSimulationState();

  /**
   * Logs the times recorded by the timers.
   * When MPI is enabled it acumulates the times (user time) of all ranks. In this case, the total
   * time will exceed the wall-clock time.
   */
  void logMeasurements();

  /**
   * Calculates the pairwise forces between particles in the autopas container.
   * @return Tells the user if the current iteration of force calculations was a tuning iteration.
   */
  bool calculatePairwiseForces();

  /**
   * Calculates the triwise forces between particles in the autopas container.
   * @return Tells the user if the current iteration of force calculations was a tuning iteration.
   */
  bool calculateTriwiseForces();

  /**
   * Adds global forces to the particles in the container.
   * @param globalForce The global force which will be applied to each particle in the container.
   */
  void calculateGlobalForces(const std::array<double, 3> &globalForce);

  /**
   * Indicates if enough iterations were completed yet.
   * Uses class member variables.
   * @return
   */
  [[nodiscard]] bool needsMoreIterations() const;

  /**
   * Checks if the global number of particles is the expected value. If not an exception is thrown.
   * If the simulation contains e.g. outflow boundaries this function does nothing!
   * @note This function is primarily for debugging purposes as it triggers global communication.
   * @param expectedNumParticlesGlobal Expected global value.
   * @param numParticlesCurrentlyMigratingLocal Number of particles that are currently not inserted but should be
   * re-inserted. E.g. immigrants and / or emigrants.
   * @param lineNumber Will be shown in the Exception so that it is easier to find the offending call. Pass __LINE__.
   */
  [[maybe_unused]] void checkNumParticles(size_t expectedNumParticlesGlobal, size_t numParticlesCurrentlyMigratingLocal,
                                          int lineNumber);

  /**
   *
   * Apply the functor chosen and configured via _configuration to the given lambda function f(auto functor).
   * @note This templated function is private and hence implemented in the .cpp
   *
   * @tparam ReturnType Return type of f.
   * @tparam FunctionType Function type ReturnType f(auto functor).
   * @param f lambda function.
   * @return Return value of f.
   */
  template <class ReturnType, class FunctionType>
  ReturnType applyWithChosenFunctor(FunctionType f);

  /**
   *
   * Apply the functor chosen and configured via _configuration to the given lambda function f(auto functor).
   * @note This templated function is private and hence implemented in the .cpp
   *
   * @tparam ReturnType Return type of f.
   * @tparam FunctionType Function type ReturnType f(auto functor).
   * @param f lambda function.
   * @return Return value of f.
   */
  template <class ReturnType, class FunctionType>
  ReturnType applyWithChosenFunctor3B(FunctionType f);
};

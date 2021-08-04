/**
 * @file Simulation.h
 * @author N. Fottner
 * @date 12.05.2019
 */
#pragma once

#include "TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/Timer.h"
#include "parsing/MDFlexConfig.h"
#ifdef AUTOPAS_INTERNODE_TUNING
#include <mpi.h>
#endif

/**
 * The main simulation class.
 */
class Simulation {
 public:
  /**
   * Particle type used for the simulation.
   */
  using ParticleType = ::ParticleType;

  /**
   * Constructor.
   *
   * This starts the timer for total simulation time.
   */
  explicit Simulation() { _timers.total.start(); };

  /**
   * Destructor.
   *
   * Closes the log file if applicable.
   */
  ~Simulation() = default;

  /**
   * Writes a VTK file for the current state of the AutoPas object.
   * @param autopas
   */
  void writeVTKFile(autopas::AutoPas<ParticleType> &autopas);

  /**
   * Initializes the ParticlePropertiesLibrary with properties from _config.
   */
  void initializeParticlePropertiesLibrary();

  /**
   * Initializes the AutoPas Object with the given config and initializes the simulation domain with the Object
   * Generators.
   * @param mdFlexConfig
   * @param autopas
   */
  void initialize(const MDFlexConfig &mdFlexConfig, autopas::AutoPas<ParticleType> &autopas);

  /**
   * Calculates the pairwise forces in the system and measures the runtime.
   * @tparam Force Calculation Functor
   * @param autopas
   */
  template <class FunctorType>
  void calculateForces(autopas::AutoPas<ParticleType> &autopas);

  /**
   * Calculate influences from global, non pairwise forces, e.g. gravity.
   * @param autopas
   */
  void globalForces(autopas::AutoPas<ParticleType> &autopas);

  /**
   * This function processes the main simulation loop
   * -calls the time discretization class(calculate fores, etc ...)
   * -do the output each timestep
   * -collects the duration of every Calculation(Position,Force,Velocity)
   * @param autopas
   */
  void simulate(autopas::AutoPas<ParticleType> &autopas);

  /**
   * Indicates if enough iterations were completed yet.
   * Uses class member variables.
   * @return
   */
  [[nodiscard]] bool needsMoreIterations() const;

  /**
   * Returns Suffix for the mpi rank the process is running on.
   * Otherwise returns empty string.
   * @return suffix
   */
  [[nodiscard]] std::string getMPISuffix() const;

  /**
   * Gives an estimate of how many iterations the simulation will do in total.
   * @return The estimate and true iff this is the correct number and not only an estimate.
   */
  [[nodiscard]] std::tuple<size_t, bool> estimateNumIterations() const;

  /**
   * Prints statistics like duration of calculation etc of the Simulation.
   * @param autopas
   */
  void printStatistics(autopas::AutoPas<ParticleType> &autopas);

  /**
   * Getter for ParticlePropertiesLibrary of Simulation.
   * @note Used for testing.
   * @return unique_prt(ParticlePropertiesLibrary)
   */
  [[nodiscard]] const std::unique_ptr<ParticlePropertiesLibrary<double, size_t>> &getPpl() const;

 private:
  /**
   * Print a progressbar and progress information to the console.
   * @note Calling this function deletes the whole current line in the terminal.
   *
   * @param iterationProgress
   * @param maxIterations
   * @param maxIsPrecise Indicate whether maxIterations is precise (true) or an estimate (false).
   */
  void printProgress(size_t iterationProgress, size_t maxIterations, bool maxIsPrecise);

  /**
   * Use the variant of the LJFunctor that uses the shifted Lennard-Jones potential.
   */
  constexpr static bool _shifting = true;
  /**
   * Use the variant of the LJFunctor that supports mixing of particle types.
   */
  constexpr static bool _mixing = true;

  /**
   * Configuration of the scenario.
   */
  std::shared_ptr<MDFlexConfig> _config;

  /**
   * Container for properties (e.g. mass) per particle type.
   */
  std::unique_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;

  /**
   * Struct containing all timers used for the simulation.
   */
  struct timers {
    autopas::utils::Timer positionUpdate, forceUpdateTotal, forceUpdatePairwise, forceUpdateGlobal, forceUpdateTuning,
        forceUpdateNonTuning, velocityUpdate, simulate, vtk, init, total, thermostat, boundaries;
  } _timers;

  /**
   * Convert a time and a name to a properly formatted string.
   * @param name incl. offset.
   * @param timeNS in nanoseconds.
   * @param numberWidth Width to which the time should be offset.
   * @param maxTime if passed the percentage of timeNS of maxTime is appended.
   * @return formatted std::string
   */
  static std::string timerToString(const std::string &name, long timeNS, size_t numberWidth = 0, long maxTime = 0);

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

//  std::string _homoName;

  //autopas::AutoPas_MPI_Comm _comm{AUTOPAS_MPI_COMM_NULL};

};

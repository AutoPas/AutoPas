/**
 * @file MDFlexSimulation.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "../domainDecomposition/DomainDecomposition.h"
#include "../configuration/MDFlexConfig.h"
#include "../TypeDefinitions.h"
#include "autopas/AutoPas.h"

#include <iostream>
#include <memory>
#include <string>
#include <tuple>

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

	virtual void initializeDomainDecomposition(int &dimensionCount) = 0;

	virtual DomainDecomposition* getDomainDecomposition() = 0;

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
  ~MDFlexSimulation();

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

	//std::shared_ptr<DomainDecomposition> _domainDecomposition;

	/**
	* The the nodes' autopas container used for simulation.
	* This member will not be initialized by the constructor and therefore has to be initialized by the deriving class.
	*/
	std::shared_ptr<autopas::AutoPas<ParticleType>> _autoPasContainer;

  std::shared_ptr<std::ofstream> _logFile;

  std::ostream* _outputStream;

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
	
	void initialize(int dimensionCount, int argc, char **argv);
	bool needsMoreIterations();
	std::tuple<size_t, bool> estimateNumberOfIterations() const;
	void printProgress(size_t iterationProgress, size_t maxIterations, bool maxIsPrecise);
	void printStatistics();
	void writeVTKFile();
	std::string getMPISuffix();
	std::string timerToString(const std::string &name, long timeNS, size_t numberWidth, long maxTime);
	void updatePositions();
	void updateForces();
	void updateVelocities();

	private:
		void initializeAutoPasContainer();
		void initializeObjects();
};

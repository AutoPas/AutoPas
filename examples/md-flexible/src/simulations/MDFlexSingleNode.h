/**
 * @file MDFlexSingleNode.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

/**
 * Runs the MD-Flex simulation on a single node.
 * This is the default demonstration of AutoPas.
 */
class MDFlexSingleNode : MDFlexSimulation {
 public:
  MDFlexSingleNode(int argc, char **argv);
  ~MDFlexSingleNode();

  /**
   * Runs the simulation
   */
  void run() override;

 private:
  /**
   * Stores the argument count passed to the constructor for later reuse.
   */
  int _argc;

  /**
   * Stores the arguments passed to the constructor for later reuse.
   */
  char **_argv;

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
   * Use the variant of the LJFunctor that uses the shifted Lennard-Jones potential.
   */
  constexpr static bool _shifting = true;
  /**
   * Use the variant of the LJFunctor that supports mixing of particle types.
   */
  constexpr static bool _mixing = true;


  /**
   * Print a progressbar and progress information to the console.
   * @note Calling this function deletes the whole current line in the terminal.
   * @param iterationProgress
   * @param maxIterations
   * @param maxIsPrecise Indicate whether maxIterations is precise (true) or an estimate (false).
   */
  void printProgress(size_t iterationProgress, size_t maxIterations, bool maxIsPrecise);

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
   * Writes a VTK file for the current state of the AutoPas object.
   */
  void writeVTKFile();

  /**
   * Calculates the pairwise forces in the system and measures the runtime.
   * @tparam Force Calculation Functor
   */
  template <class FunctorType>
  void calculateForces();

  /**
   * Calculate influences from global, non pairwise forces, e.g. gravity.
   */
  void globalForces();

  /**
   * Indicates if enough iterations were completed yet.
   * Uses class member variables.
   * @return
   */
  [[nodiscard]] bool needsMoreIterations() const;

  /**
   * Gives an estimate of how many iterations the simulation will do in total.
   * @return The estimate and true iff this is the correct number and not only an estimate.
   */
  [[nodiscard]] std::tuple<size_t, bool> estimateNumIterations() const;

  /**
   * Prints statistics like duration of calculation etc of the Simulation.
   */
	void printStatistics();
};

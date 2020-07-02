/**
 * @file MDFlexConfig.h
 * @author F. Gratl
 * @date 18.10.2019
 */

#pragma once

#include <getopt.h>

#include <map>
#include <set>
#include <utility>

#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/utils/NumberSet.h"
#include "src/Objects/CubeGauss.h"
#include "src/Objects/CubeGrid.h"
#include "src/Objects/CubeUniform.h"
#include "src/Objects/Sphere.h"

/**
 * Class containing all necessary parameters for configuring a md-flexible simulation.
 */
class MDFlexConfig {
 public:
  MDFlexConfig() = default;

  /**
   * Struct to bundle information for options.
   * @tparam T Datatype of the option
   * @tparam getOptChar int for the switch case that is used during cli argument parsing with getOpt.
   * @note ints should be unique so they can be used for a switch case.
   * @note getOptChar should never be -1 because getopt uses this value to indicate that there are no more cli arguments
   * @note use __LINE__ as a cheap generator for unique ints.
   * @todo c++20: With support for non-type template parameters replace the template getOptChar with the name string
   * Then the getOptChar can be generated from that (e.g. with a hash) at compile time. Discussed here:
   * https://github.com/AutoPas/AutoPas/pull/469#discussion_r431270944
   */
  template <class T, int getOptChar>
  struct MDFlexOption {
    /**
     * Value of this option.
     */
    T value;

    /**
     * Indicate whether this option is a flag or takes arguments.
     */
    bool requiresArgument;

    /**
     * String representation of the option name.
     */
    std::string name;

    /**
     * String describing this option. This is displayed when md-flexible is invoked with --help.
     */
    std::string description;

    /**
     * Member to access the template parameter.
     */
    constexpr static int getoptChar{getOptChar};

    /**
     * Constructor
     * @param value Default value for this option.
     * @param newName String representation of the option name.
     * @param requiresArgument Indicate whether this option is a flag or takes arguments.
     * @param newDescription String describing this option. This is displayed when md-flexible is invoked with --help.
     */
    MDFlexOption(T value, std::string newName, bool requiresArgument, std::string newDescription)
        : requiresArgument(requiresArgument),
          name(std::move(newName)),
          description(std::move(newDescription)),
          value(std::move(value)) {}

    /**
     * Returns a getopt option struct for this object.
     * @return
     */
    [[nodiscard]] struct option toGetoptOption() const {
      struct option retStruct {
        name.c_str(), requiresArgument, nullptr, getOptChar
      };
      return retStruct;
    }
  };

  /**
   * Convert the content of the config to a string representation.
   * @return
   */
  [[nodiscard]] std::string to_string() const;

  /**
   * Checks parsed Objects and determines the necessary size of the simulation box.
   */
  void calcSimulationBox();

  /**
   * Adds parameters to all relevant particle property attributes and checks if the type already exists.
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   */
  void addParticleType(unsigned long typeId, double epsilon, double sigma, double mass);

  /**
   * Choice of the functor
   */
  enum class FunctorOption { lj12_6, lj12_6_AVX, lj12_6_Globals };

  /**
   * Choice of the particle generators specified in the command line
   */
  enum class GeneratorOption { grid, uniform, gaussian, sphere };

  //  All options in the config
  //  Make sure that the description is parsable by `CLIParser::createZSHCompletionFile()`!

  /**
   * yamlFilename
   */
  MDFlexOption<std::string, __LINE__> yamlFilename{"", "yaml-filename", true, "Path to a .yaml input file."};

  // AutoPas options:
  /**
   * containerOptions
   */
  MDFlexOption<std::set<autopas::ContainerOption>, __LINE__> containerOptions{
      autopas::ContainerOption::getMostOptions(), "container", true,
      "List of container options to use. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::ContainerOption::getAllOptions(), " ", {"(", ")"})};
  /**
   * dataLayoutOptions
   */
  MDFlexOption<std::set<autopas::DataLayoutOption>, __LINE__> dataLayoutOptions{
      autopas::DataLayoutOption::getMostOptions(), "data-layout", true,
      "List of data layout options to use. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::DataLayoutOption::getAllOptions(), " ", {"(", ")"})};
  /**
   * selectorStrategy
   */
  MDFlexOption<autopas::SelectorStrategyOption, __LINE__> selectorStrategy{
      autopas::SelectorStrategyOption::fastestAbs, "selector-strategy", true,
      "Strategy how to reduce the sample measurements to a single value. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::SelectorStrategyOption::getAllOptions(), " ", {"(", ")"})};
  /**
   * traversalOptions
   */
  MDFlexOption<std::set<autopas::TraversalOption>, __LINE__> traversalOptions{
      autopas::TraversalOption::getMostOptions(), "traversal", true,
      "List of traversal options to use. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::TraversalOption::getAllOptions(), " ", {"(", ")"})};
  /**
   * traversalOptions
   */
  MDFlexOption<std::set<autopas::LoadEstimatorOption>, __LINE__> loadEstimatorOptions{
      autopas::LoadEstimatorOption::getMostOptions(), "load-estimator", true,
      "List of load estimator function choices for traversals that do heuristic load balancing. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::LoadEstimatorOption::getAllOptions(), " ", {"(", ")"})};
  /**
   * newton3Options
   */
  MDFlexOption<std::set<autopas::Newton3Option>, __LINE__> newton3Options{
      autopas::Newton3Option::getMostOptions(), "newton3", true,
      "List of newton3 options to use. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::Newton3Option::getAllOptions(), " ", {"(", ")"})};
  /**
   * cellSizeFactors
   */
  MDFlexOption<std::shared_ptr<autopas::NumberSet<double>>, __LINE__> cellSizeFactors{
      std::make_shared<autopas::NumberSetFinite<double>>(std::set<double>{1.}), "cell-size", true,
      "Factor for the interaction length to determine the cell size."};
  /**
   * logFileName
   */
  MDFlexOption<std::string, __LINE__> logFileName{"", "log-file", true,
                                                  "Path to an .out file to store the log output."};
  /**
   * logLevel
   */
  MDFlexOption<autopas::Logger::LogLevel, __LINE__> logLevel{
      autopas::Logger::LogLevel::info, "log-level", true,
      "Log level for AutoPas. Set to debug for tuning information. "
      "Possible Values: (trace debug info warn error critical off)"};
  /**
   * tuningStrategyOption
   */
  MDFlexOption<autopas::TuningStrategyOption, __LINE__> tuningStrategyOption{
      autopas::TuningStrategyOption::fullSearch, "tuning-strategy", true,
      "Strategy how to reduce the sample measurements to a single value. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::TuningStrategyOption::getAllOptions(), " ", {"(", ")"})};
  /**
   * tuningInterval
   */
  MDFlexOption<unsigned int, __LINE__> tuningInterval{5000, "tuning-interval", true,
                                                      "Number of iterations between two tuning phases."};
  /**
   * tuningSamples
   */
  MDFlexOption<unsigned int, __LINE__> tuningSamples{3, "tuning-samples", true,
                                                     "Number of samples to collect per configuration."};
  /**
   * tuningMaxEvidence
   */
  MDFlexOption<unsigned int, __LINE__> tuningMaxEvidence{
      10, "tuning-max-evidence", true,
      "For Bayesian based tuning strategies: Maximum number of evidences "
      "tuning strategies that have no finishing indicator take."};
  /**
   * relativeOptimumRange
   */
  MDFlexOption<double, __LINE__> relativeOptimumRange{
      1.2, "relative-optimum-range", true,
      "For predictive based tuning strategies: Configurations whose predicted performance lies within this range of "
      "the predicted optimal performance will be tested."};
  /**
   * maxTuningPhasesWithoutTest
   */
  MDFlexOption<unsigned int, __LINE__> maxTuningPhasesWithoutTest{
      10, "max-tuning-phases-without-test", true,
      "For predictive based tuning strategies: Maximal number of "
      "tuning phases a configurations can be excluded from testing."};
  /**
   * relativeRangeforBlacklist
   */
  MDFlexOption<unsigned int, __LINE__> relativeRangeForBlacklist{
      0, "relative-range-for-blacklist", true,
      "For predictive based tuning strategies: Relative range to the optimum in which the first evidence of a "
      "configuration needs to be ot not get put on the blacklist. If the blacklist should not be used give 0 as an "
      "argument"};
  /**
   * evidenceFirstPrediction
   */
  MDFlexOption<unsigned int, __LINE__> evidenceFirstPrediction{
      3, "evidence-for-prediction", true,
      "For predictive based tuning strategies: The number of evidence for a configuration that needs to be gathered "
      "before the first prediction is made. This number also determines how much evidence is used for the calculation "
      "of the prediction and for a polynomial extrapolation method, this is also the degree of the polynomial."};
  /**
   * extrapolationMethodOption
   */
  MDFlexOption<autopas::ExtrapolationMethodOption, __LINE__> extrapolationMethodOption{
      autopas::ExtrapolationMethodOption::linearRegression, "extrapolation-method", true,
      "For predictive based tuning strategies: The extrapolation method that calculates the prediction. Possible "
      "Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::ExtrapolationMethodOption::getAllOptions(), " ", {"(", ")"})};
  /**
   * vtkFileName
   */
  MDFlexOption<std::string, __LINE__> vtkFileName{"", "vtk-filename", true, "Basename for all VTK output files."};
  /**
   * vtkWriteFrequency
   */
  MDFlexOption<size_t, __LINE__> vtkWriteFrequency{100, "vtk-write-frequency", true,
                                                   "Number of iterations after which a VTK file is written."};
  /**
   * verletClusterSize
   */
  MDFlexOption<unsigned int, __LINE__> verletClusterSize{4, "verlet-cluster-size", true,
                                                         "Number of particles in Verlet clusters."};
  /**
   * verletRebuildFrequency
   */
  MDFlexOption<unsigned int, __LINE__> verletRebuildFrequency{
      20, "verlet-rebuild-frequency", true, "Number of iterations after which containers are rebuilt."};
  /**
   * verletSkinRadius
   */
  MDFlexOption<double, __LINE__> verletSkinRadius{.2, "verlet-skin-radius", true,
                                                  "Skin added to the cutoff to form the interaction length."};
  /**
   * boxMin
   */
  MDFlexOption<std::array<double, 3>, 0> boxMin{
      {0, 0, 0}, "box-min", true, "Lower front left corner of the simulation box."};
  /**
   * boxMax
   */
  MDFlexOption<std::array<double, 3>, 0> boxMax{
      {1, 1, 1}, "box-max", true, "Upper back right corner of the simulation box."};
  /**
   * acquisitionFunctionOption
   */
  MDFlexOption<autopas::AcquisitionFunctionOption, __LINE__> acquisitionFunctionOption{
      autopas::AcquisitionFunctionOption::upperConfidenceBound, "tuning-acquisition-function", true,
      "For Bayesian based tuning strategies: Function to determine the predicted knowledge gain when testing a given "
      "configuration. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::AcquisitionFunctionOption::getAllOptions(), " ", {"(", ")"})};

  // Simulation Options:
  /**
   * cutoff
   */
  MDFlexOption<double, __LINE__> cutoff{2., "cutoff", true, "Lennard-Jones force cutoff."};
  /**
   * functorOption
   */
  MDFlexOption<FunctorOption, __LINE__> functorOption{
      FunctorOption::lj12_6, "functor", true,
      "Force functor to use. Possible Values: (lennard-jones lennard-jones-AVX2 lennard-jones-globals)"};
  /**
   * iterations
   */
  MDFlexOption<size_t, __LINE__> iterations{10, "iterations", true, "Number of iterations to simulate."};
  /**
   * tuningPhases
   */
  MDFlexOption<size_t, __LINE__> tuningPhases{
      0, "tuning-phases", true, "Number of tuning phases to simulate. This option overwrites --iterations."};
  /**
   * Periodic boundaries.
   * This starts with a "not" such that it can be used as a flag with a sane default.
   */
  MDFlexOption<bool, __LINE__> periodic{
      true, "periodic-boundaries", true,
      "(De)Activate periodic boundaries. Possible Values: (true false) Default: true."};
  /**
   * dontMeasureFlops
   */
  MDFlexOption<bool, __LINE__> dontMeasureFlops{true, "no-flops", false, "Set to omit the calculation of flops."};
  /**
   * Omit the creation of a config file at the end of the Simulation.
   * This starts with a "not" such that it can be used as a flag with a sane default.
   */
  MDFlexOption<bool, __LINE__> dontCreateEndConfig{
      true, "no-end-config", false, "Set to omit the creation of a yaml file at the end of a simulation."};
  /**
   * deltaT
   */
  MDFlexOption<double, __LINE__> deltaT{0.001, "deltaT", true,
                                        "Length of a timestep. Set to 0 to deactivate time integration."};
  /**
   * epsilonMap
   */
  MDFlexOption<std::map<unsigned long, double>, 0> epsilonMap{
      {{0ul, 1.}}, "particle-epsilon", true, "Mapping from particle type to an epsilon value."};
  /**
   * sigmaMap
   */
  MDFlexOption<std::map<unsigned long, double>, 0> sigmaMap{
      {{0ul, 1.}}, "particle-sigma", true, "Mapping from particle type to a sigma value."};
  /**
   * massMap
   */
  MDFlexOption<std::map<unsigned long, double>, 0> massMap{
      {{0ul, 1.}}, "particle-mass", true, "Mapping from particle type to a mass value."};

  // Options for additional Object Generation on command line
  /**
   * boxLength
   */
  MDFlexOption<double, __LINE__> boxLength{10, "box-length", true, "Length of the simulation box as a cuboid."};
  /**
   * distributionMean
   */
  MDFlexOption<std::array<double, 3>, __LINE__> distributionMean{
      {5., 5., 5.}, "distribution-mean", true, "Mean of the gaussian distribution for random particle initialization."};
  /**
   * distributionStdDev
   */
  MDFlexOption<std::array<double, 3>, __LINE__> distributionStdDev{
      {2., 2., 2.},
      "distribution-stddeviation",
      true,
      "Standard deviation of the gaussian distribution for random particle initialization."};
  /**
   * particlesPerDim
   */
  MDFlexOption<size_t, __LINE__> particlesPerDim{10, "particles-per-dimension", true,
                                                 "Size of the scenario for the grid generator."};
  /**
   * particlesTotal
   */
  MDFlexOption<size_t, __LINE__> particlesTotal{
      1000, "particles-total", true, "Total number of particles for the random distribution based generators."};
  /**
   * particleSpacing
   * For a stable grid initialize this as 2^(1/6) sigma
   */
  MDFlexOption<double, __LINE__> particleSpacing{1.1225 * 1, "particle-spacing", true,
                                                 "Space between two particles for the grid generator."};
  /**
   * generatorOption
   */
  MDFlexOption<GeneratorOption, __LINE__> generatorOption{
      GeneratorOption::grid, "particle-generator", true,
      "Scenario generator. Possible Values: (grid uniform gaussian sphere) Default: grid"};

  // Object Generation:
  /**
   * objectsStr
   */
  static inline const char *objectsStr{"Objects"};
  /**
   * bottomLeftBackCornerStr
   */
  static inline const char *bottomLeftBackCornerStr{"bottomLeftCorner"};
  /**
   * velocityStr
   */
  static inline const char *velocityStr{"velocity"};
  /**
   * particleTypeStr
   */
  static inline const char *particleTypeStr{"particle-type"};
  /**
   * particlesPerObjectStr
   */
  static inline const char *particlesPerObjectStr{"numberOfParticles"};
  /**
   * cubeGridObjectsStr
   */
  static inline const char *cubeGridObjectsStr{"CubeGrid"};
  /**
   * cubeGridObjects
   */
  std::vector<CubeGrid> cubeGridObjects{};
  /**
   * cubeGaussObjectsStr
   */
  static inline const char *cubeGaussObjectsStr{"CubeGauss"};
  /**
   * cubeGaussObjects
   */
  std::vector<CubeGauss> cubeGaussObjects{};
  /**
   * cubeUniformObjectsStr
   */
  static inline const char *cubeUniformObjectsStr{"CubeUniform"};
  /**
   * cubeUniformObjects
   */
  std::vector<CubeUniform> cubeUniformObjects{};
  /**
   * sphereObjectsStr
   */
  static inline const char *sphereObjectsStr{"Sphere"};
  /**
   * sphereCenterStr
   */
  static inline const char *sphereCenterStr{"center"};
  /**
   * sphereRadiusStr
   */
  static inline const char *sphereRadiusStr{"radius"};
  /**
   * sphereObjects
   */
  std::vector<Sphere> sphereObjects{};

  // Thermostat Options
  /**
   * useThermostat
   */
  MDFlexOption<bool, __LINE__> useThermostat{
      false, "thermostat", true,
      "(De)Activate the thermostat. Only useful when used to overwrite a yaml file. "
      "Possible Values: (true false) Default: false"};
  /**
   * initTemperature
   */
  MDFlexOption<double, 0> initTemperature{0., "initialTemperature", true,
                                          "Thermostat option. Initial temperature of the system."};
  /**
   * targetTemperature
   */
  MDFlexOption<double, 0> targetTemperature{0., "targetTemperature", true,
                                            "Thermostat option. Target temperature of the system."};
  /**
   * deltaTemp
   */
  MDFlexOption<double, 0> deltaTemp{0., "deltaTemperature", true,
                                    "Thermostat option. Maximal temperature jump the thermostat is allowed to apply."};
  /**
   * thermostatInterval
   */
  MDFlexOption<size_t, 0> thermostatInterval{
      0, "thermostatInterval", true,
      "Thermostat option. Number of Iterations between two applications of the thermostat."};
  /**
   * addBrownianMotion
   */
  MDFlexOption<bool, 0> addBrownianMotion{
      true, "addBrownianMotion", true,
      "Thermostat option. Whether the particle velocities should be initialized using "
      "Brownian motion. Possible Values: (true false) Default: true"};

  // Checkpoint Options
  /**
   * checkpointfile
   */
  MDFlexOption<std::string, __LINE__> checkpointfile{"", "checkpoint", true,
                                                     "Path to a .vtk File to load as a checkpoint."};

  /**
   * valueOffset used for cli-output alignment
   */
  static constexpr size_t valueOffset{33};
};

/**
 * Stream insertion operator for MDFlexConfig.
 * @param os
 * @param config
 * @return
 */
inline std::ostream &operator<<(std::ostream &os, const MDFlexConfig &config) { return os << config.to_string(); }

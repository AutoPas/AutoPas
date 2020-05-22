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
   * @tparam T
   */
  struct MDFlexOptionInterface {
    /**
     * Indicate whether this option is a flag or takes arguments.
     */
    bool requiresArgument;
    /**
     * Char for the switch case that is used during cli argument parsing with getOpt.
     * @note Set to 0 (not '0') if not intended to be used from command line.
     */
    char getoptSwitchChar;
    /**
     * String representation of the option name.
     */
    std::string name;
    /**
     * String describing this option. This is displayed when md-flexible is invoked with --help.
     */
    std::string description;

    /**
     * Constructor
     * @param value Default value for this option.
     * @param newName String representation of the option name.
     * @param requiresArgument Indicate whether this option is a flag or takes arguments.
     * @param getOptSwitchChar Char for the switch case that is used during cli argument parsing with getOpt.
     * @param newDescription String describing this option. This is displayed when md-flexible is invoked with --help.
     */
    MDFlexOptionInterface(std::string newName, bool requiresArgument, char getOptSwitchChar,
                 std::string newDescription)
        : requiresArgument(requiresArgument), getoptSwitchChar(getOptSwitchChar), name(std::move(newName)), description(std::move(newDescription)) {
    }

    /**
     * Returns a getopt option struct for this object.
     * @return
     */
    [[nodiscard]] struct option toGetoptOption() const {
      struct option retStruct{name.c_str(), requiresArgument, nullptr, getoptSwitchChar};
      return retStruct;
    }
  };

  /**
   * Struct to bundle information for options.
   * @tparam T
   */
  template <class T>
  struct MDFlexOption : MDFlexOptionInterface {
    /**
     * Value of this option.
     */
    T value;

    /**
     * Constructor
     * @param value Default value for this option.
     * @param newName String representation of the option name.
     * @param requiresArgument Indicate whether this option is a flag or takes arguments.
     * @param getOptSwitchChar Char for the switch case that is used during cli argument parsing with getOpt.
     * @param newDescription String describing this option. This is displayed when md-flexible is invoked with --help.
     */
    MDFlexOption(T value, const std::string& newName, bool requiresArgument, char getOptSwitchChar,
                 std::string newDescription)
        : MDFlexOptionInterface(newName, requiresArgument, getOptSwitchChar, newDescription), value(std::move(value)){
    }

  };

  /**
   * Convert the content of the config to a string representation.
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

  MDFlexOption<std::string> yamlFilename{"", "yaml-filename", true, 'z', "Path to input file."};

  // AutoPas options:
  MDFlexOption<std::set<autopas::ContainerOption>> containerOptions{
      autopas::ContainerOption::getAllOptions(), "container", true, 'c',
      "List of container options to use. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::ContainerOption::getAllOptions())};
  MDFlexOption<std::set<autopas::DataLayoutOption>> dataLayoutOptions{
      autopas::DataLayoutOption::getAllOptions(), "data-layout", true, 'd',
      "List of data layout options to use. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::DataLayoutOption::getAllOptions())};
  MDFlexOption<autopas::SelectorStrategyOption> selectorStrategy{
      autopas::SelectorStrategyOption::fastestAbs, "selector-strategy", true, 'y',
      "Strategy how to reduce the sample measurements to a single value. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::SelectorStrategyOption::getAllOptions())};
  MDFlexOption<std::set<autopas::TraversalOption>> traversalOptions{
      autopas::TraversalOption::getAllOptions(), "traversal", true, 't',
      "List of traversal options to use. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::TraversalOption::getAllOptions())};
  MDFlexOption<std::set<autopas::Newton3Option>> newton3Options{
      autopas::Newton3Option::getAllOptions(), "newton3", true, '3',
      "List of newton3 options to use. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::Newton3Option::getAllOptions())};
  MDFlexOption<std::shared_ptr<autopas::NumberSet<double>>> cellSizeFactors{
      std::make_shared<autopas::NumberSetFinite<double>>(std::set<double>{1.}), "cell-size", true, 'a',
      "Factor for the interaction length to determine the cell size."};
  MDFlexOption<std::string> logFileName{"", "log-file", true, 'L',
                                        "Path to a file to store the log output. Possible Values: "};
  MDFlexOption<autopas::Logger::LogLevel> logLevel{autopas::Logger::LogLevel::info, "log-level", true, 'l',
                                                   "Log level for AutoPas. Set to debug for tuning information."};
  MDFlexOption<autopas::TuningStrategyOption> tuningStrategyOption{
      autopas::TuningStrategyOption::fullSearch, "tuning-strategy", true, 'T',
      "Strategy how to reduce the sample measurements to a single value."};
  MDFlexOption<unsigned int> tuningInterval{100, "tuning-interval", true, 'I',
                                            "Number of iterations between two tuning phases."};
  MDFlexOption<unsigned int> tuningSamples{3, "tuning-samples", true, 'S',
                                           "Number of samples to collect per configuration."};
  MDFlexOption<unsigned int> tuningMaxEvidence{10, "tuning-max-evidence", true, 'E',
                                               "For Bayesian based tuning strategies: Maximum number of evidences "
                                               "tuning strategies that have no finishing indicator take."};
  MDFlexOption<double> relativeOptimumRange{
      1.2, "relative-optimum-range", true, 'o',
      "For predictive based tuning strategies: Configurations whose predicted performance lies within this range of "
      "the predicted optimal performance will be tested."};
  MDFlexOption<unsigned int> maxTuningPhasesWithoutTest{5, "max-tuning-phases-without-test", true, 'M',
                                                        "For predictive based tuning strategies: Maximal number of "
                                                        "tuning phases a configurations can be excluded from testing."};
  MDFlexOption<std::string> vtkFileName{"", "vtk-filename", true, 'w', "Basename for all VTK output files."};
  MDFlexOption<size_t> vtkWriteFrequency{100, "vtk-write-frequency", true, 'W',
                                         "Number of iterations after which a VTK file is written."};
  MDFlexOption<unsigned int> verletClusterSize{4, "verlet-cluster-size", true, 'q',
                                               "Number of particles in Verlet clusters."};
  MDFlexOption<unsigned int> verletRebuildFrequency{1, "verlet-rebuild-frequency", true, 'v',
                                                    "Number of iterations after which containers are rebuilt."};
  MDFlexOption<double> verletSkinRadius{.2, "verlet-skin-radius", true, 'r',
                                        "Skin added to the cutoff to form the interaction length."};
  MDFlexOption<std::array<double, 3>> boxMin{
      {0, 0, 0}, "box-min", true, 0, "Lower front left corner of the simulation box."};
  MDFlexOption<std::array<double, 3>> boxMax{
      {5, 5, 5}, "box-max", true, 0, "Upper back right corner of the simulation box."};
  MDFlexOption<autopas::AcquisitionFunctionOption> acquisitionFunctionOption{
      autopas::AcquisitionFunctionOption::lowerConfidenceBound, "tuning-acquisition-function", true, 'A',
      "For Bayesian based tuning strategies: Function to determine the predicted knowledge gain when testing a given "
      "configuration. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::AcquisitionFunctionOption::getAllOptions())};

  // Simulation Options:
  MDFlexOption<double> cutoff{1., "cutoff", true, 'C', "Lennard-Jones force cutoff."};
  MDFlexOption<FunctorOption> functorOption{
      FunctorOption::lj12_6, "functor", true, 'f',
      "Force functor to use. Possible Values: [lennard-jones, lennard-jones-AVX2, lennard-jones-globals]"};
  MDFlexOption<size_t> iterations{10, "iterations", true, 'i', "Number of iterations to simulate."};
  MDFlexOption<size_t> tuningPhases{0, "tuning-phases", true, 'P',
                                    "Number of tuning phases to simulate. This option overwrites --iterations."};
  MDFlexOption<bool> periodic{true, "periodic-boundaries", true, 'p',
                              "(De)Activate periodic boundaries. Possible Values: [true, false] Default: true."};
  // this starts with a "not" such that it can be used as a flag with a sane default.
  MDFlexOption<bool> dontMeasureFlops{true, "no-flops", false, 'F', "Set to omit the calculation of flops."};
  // this starts with a "not" such that it can be used as a flag with a sane default.
  MDFlexOption<bool> dontCreateEndConfig{true, "no-end-config", false, 'e',
                                         "Set to omit the creation of a yaml file at the end of a simulation."};
  MDFlexOption<double> deltaT{0.001, "deltaT", true, 'D',
                              "Length of a timestep. Set to 0 to deactivate time integration."};
  MDFlexOption<std::map<unsigned long, double>> epsilonMap{
      {{0ul, 1.}}, "particle-epsilon", true, 0, "Mapping from particle type to an epsilon value."};
  MDFlexOption<std::map<unsigned long, double>> sigmaMap{
      {{0ul, 1.}}, "particle-sigma", true, 0, "Mapping from particle type to a sigma value."};
  MDFlexOption<std::map<unsigned long, double>> massMap{
      {{0ul, 1.}}, "particle-mass", true, 0, "Mapping from particle type to a mass value."};

  // Options for additional Object Generation on command line
  MDFlexOption<double> boxLength{10, "box-length", true, 'b', "Length of the simulation box as a cuboid."};
  MDFlexOption<std::array<double, 3>> distributionMean{
      {5., 5., 5.},
      "distribution-mean",
      true,
      'm',
      "Mean of the gaussian distribution for random particle initialization."};
  MDFlexOption<std::array<double, 3>> distributionStdDev{
      {2., 2., 2.},
      "distribution-stddeviation",
      true,
      'z',
      "Standard deviation of the gaussian distribution for random particle initialization."};
  MDFlexOption<size_t> particlesPerDim{10, "particles-per-dimension", true, 'n',
                                       "Size of the scenario for the grid generator."};
  MDFlexOption<size_t> particlesTotal{1000, "particles-total", true, 'N',
                                      "Total number of particles for the random distribution based generators."};
  MDFlexOption<double> particleSpacing{.5, "particle-spacing", true, 's',
                                       "Space between two particles for the grid generator."};
  MDFlexOption<GeneratorOption> generatorOption{
      GeneratorOption::grid, "particle-generator", true, 'g',
      "Scenario generator. Possible Values: [grid, uniform, gaussian, sphere] Default: grid"};

  // Object Generation:
  static inline const char *objectsStr{"Objects"};
  static inline const char *bottomLeftBackCornerStr{"bottomLeftCorner"};
  static inline const char *velocityStr{"velocity"};
  static inline const char *particleTypeStr{"particle-type"};
  static inline const char *particlesPerObjectStr{"numberOfParticles"};
  static inline const char *cubeGridObjectsStr{"CubeGrid"};
  std::vector<CubeGrid> cubeGridObjects{};
  static inline const char *cubeGaussObjectsStr{"CubeGauss"};
  std::vector<CubeGauss> cubeGaussObjects{};
  static inline const char *cubeUniformObjectsStr{"CubeUniform"};
  std::vector<CubeUniform> cubeUniformObjects{};
  static inline const char *sphereObjectsStr{"Sphere"};
  static inline const char *sphereCenterStr{"center"};
  static inline const char *sphereRadiusStr{"radius"};
  std::vector<Sphere> sphereObjects{};

  // Thermostat Options
  MDFlexOption<bool> useThermostat{false, "thermostat", true, 'u',
                                   "(De)Activate the thermostat. Only useful when used to overwrite a yaml file. "
                                   "Possible Values: [true, false] Default: false"};
  MDFlexOption<double> initTemperature{0., "initialTemperature", true, 0,
                                       "Thermostat option. Initial temperature of the system."};
  MDFlexOption<double> targetTemperature{0., "targetTemperature", true, 0,
                                         "Thermostat option. Target temperature of the system."};
  MDFlexOption<double> deltaTemp{0., "deltaTemperature", true, 0,
                                 "Thermostat option. Maximal temperature jump the thermostat is allowed to apply."};
  MDFlexOption<size_t> thermostatInterval{
      0, "thermostatInterval", true, 0,
      "Thermostat option. Number of Iterations between two applications of the thermostat."};
  MDFlexOption<bool> addBrownianMotion{true, "addBrownianMotion", true, 0,
                                       "Thermostat option. Whether the particle velocities should be initialized using "
                                       "Brownian motion. Possible Options: [true, false] Default: true"};

  // Checkpoint Options
  MDFlexOption<std::string> checkpointfile{"", "checkpoint", true, '4', "Path to a VTK File to load as a checkpoint."};

  // used for cli-output alignment
  static constexpr size_t valueOffset{33};
};

/**
 * Stream insertion operator for MDFlexConfig.
 * @param os
 * @param config
 * @return
 */
inline std::ostream &operator<<(std::ostream &os, const MDFlexConfig &config) { return os << config.to_string(); }
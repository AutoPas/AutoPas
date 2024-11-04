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
#include "autopas/options/EnergySensorOption.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningMetricOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/NumberSet.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/objects/CubeClosestPacked.h"
#include "src/configuration/objects/CubeGauss.h"
#include "src/configuration/objects/CubeGrid.h"
#include "src/configuration/objects/CubeUniform.h"
#include "src/configuration/objects/Sphere.h"
#include "src/domainDecomposition/LoadBalancerOption.h"
#include "src/options/BoundaryTypeOption.h"

/**
 * Class containing all necessary parameters for configuring a md-flexible simulation.
 */
class MDFlexConfig {
 public:
  /**
   * Constructor that initializes the configuration from the CLI arguments (incl. yaml file argument).
   * @param argc: the argument count of the arguments passed to the main function.
   * @param argv: the argument vector passed to the main function.
   */
  MDFlexConfig(int argc, char **argv);

  /**
   * Constructor using only default values.
   * Useful for testing but might require setting some values before this is valid.
   */
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
          value(std::move(value)),
          description(std::move(newDescription)) {}

    /**
     * Returns a getopt option struct for this object.
     * @return
     */
    [[nodiscard]] auto toGetoptOption() const {
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
   * Returns the ParticlePropertiesLibrary containing the properties of the particle types used in this simulation.
   * @return the ParticlePropertiesLibrary
   */
  std::shared_ptr<ParticlePropertiesLibraryType> getParticlePropertiesLibrary() { return _particlePropertiesLibrary; }

  /**
   * Adds parameters of a LJ site and checks if the siteId already exists.
   *
   * For single site simulations, the molecule's molId is used to look up the site with siteId = molId.
   *
   * @param siteId unique site type id
   * @param epsilon
   * @param sigma
   * @param mass
   */
  void addSiteType(unsigned long siteId, double epsilon, double sigma, double mass);

  /**
   * Adds site positions and types for a given molecule type and checks if the molId already exists
   *
   * @note When md-flexible is compiled for only single-site molecules, calls to this function return errors.
   *
   * @param molId unique mol type id
   * @param siteIds vector of ids of site types
   * @param relSitePos vector of relative site positions
   * @param momentOfInertia diagonalized moment of inertia as a length 3 array of double representing diagonal elements
   */
  void addMolType(unsigned long molId, const std::vector<unsigned long> &siteIds,
                  const std::vector<std::array<double, 3>> &relSitePos, std::array<double, 3> momentOfInertia);

  /**
   * Flushes the particles.
   */
  void flushParticles();

  /**
   * Loads the particles from the checkpoint file defined in the configuration file.
   * If the checkpoint has been recorded using multiple processes, the rank of the current process needs to be passed.
   * The provided rank also needs to respect the domain decomposition. E. g. if the regular grid decomposition is
   * used, don't pass the MPI_COMM_WORLD rank, as it might differ from the grid rank derived in the decomposition
   * scheme. The wrong rank might result in a very bad network topology and therefore increase communication costs.
   * @param rank: The MPI rank of the current process.
   * @param numRanks: The size of the MPI communicator used for the simulation.
   */
  void loadParticlesFromCheckpoint(const size_t &rank, const size_t &numRanks);

  /**
   * Choice of the functor
   */
  enum class FunctorOption { lj12_6, lj12_6_AVX, lj12_6_SVE, lj12_6_Globals, lj12_6_XSIMD, lj12_6_MIPP, lj12_6_SIMDe, lj12_6_HWY };

  /**
   * Choice of the particle generators specified in the command line
   */
  enum class GeneratorOption { grid, uniform, gaussian, sphere, closestPacked };

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
   * loadEstimatorOptions
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

  MDFlexOption<std::set<autopas::VectorizationPatternOption>, __LINE__> vecPatternOptions {
      autopas::VectorizationPatternOption::getMostOptions(), "vectorizationPattern", true,
      "Vectorization Pattern for HWY Functor."
  };
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
   * tuningStrategyOptions
   */
  MDFlexOption<std::vector<autopas::TuningStrategyOption>, __LINE__> tuningStrategyOptions{
      {},
      "tuning-strategies",
      true,
      "Ordered pipeline of strategies to find the optimal algorithmic configuration. "
      "Leave empty to perform an exhaustive search. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(
              []() {
                auto options = autopas::TuningStrategyOption::getAllOptions();
                options.erase(autopas::TuningStrategyOption::fullSearch);
                return options;
              }(),
              " ", {"(", ")"})};
  /**
   * tuningMetricOption
   */
  MDFlexOption<autopas::TuningMetricOption, __LINE__> tuningMetricOption{
      autopas::TuningMetricOption::time, "tuning-metric", true,
      "Metric to use for tuning. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::TuningMetricOption::getAllOptions(), " ", {"(", ")"})};

  /**
   * enerySensorOption
   */
  MDFlexOption<autopas::EnergySensorOption, __LINE__> energySensorOption{
      autopas::EnergySensorOption::none, "energy-sensor", true,
      "Sensor, used for energy consumption measurement. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(autopas::EnergySensorOption::getAllOptions(), " ", {"(", ")"})};
  /**
   * ruleFilename
   */
  MDFlexOption<std::string, __LINE__> ruleFilename{
      "", "rule-filename", true, "Path to a .rule file containing rules for the rule-based tuning method."};

  /**
   * fuzzyRuleFilename
   */
  MDFlexOption<std::string, __LINE__> fuzzyRuleFilename{
      "", "fuzzy-rule-filename", true, "Path to a .frule file containing rules for the fuzzy-based tuning method."};

  /**
   * MPITuningMaxDifferenceForBucket
   */
  MDFlexOption<double, __LINE__> MPITuningMaxDifferenceForBucket{
      0.2, "mpi-tuning-max-difference-for-bucket", true,
      "For MPI-tuning: Maximum of the relative difference in the comparison metric for two ranks which exchange their "
      "tuning information."};

  /**
   * MPITuningWeightForMaxDensity
   */
  MDFlexOption<double, __LINE__> MPITuningWeightForMaxDensity{
      0.1, "mpi-tuning-weight-for-max-density", true,
      "For MPI-tuning: Weight for maxDensity in the calculation for bucket distribution."};

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
      5, "max-tuning-phases-without-test", true,
      "For predictive based tuning strategies: Maximal number of "
      "tuning phases a configurations can be excluded from testing."};
  /**
   * relativeBlacklistRange
   */
  MDFlexOption<double, __LINE__> relativeBlacklistRange{
      3., "relative-blacklist-range", true,
      "For Slow Config Filter: When the first evidence of a configuration is further away from the "
      "optimum than this relative range, the configuration is ignored for the rest of the simulation."};
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
   * vtkOutputFolder
   */
  MDFlexOption<std::string, __LINE__> vtkOutputFolder{"output", "vtk-output-folder", true,
                                                      "The location where the vtk output will be created."};
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
      15, "verlet-rebuild-frequency", true, "Number of iterations after which containers are rebuilt."};
  /**
   * verletSkinRadiusPerTimeStep
   */
  MDFlexOption<double, __LINE__> verletSkinRadiusPerTimestep{
      .2, "verlet-skin-radius-per-timestep", true,
      "Skin added to the cutoff to form the interaction length. The total skin width is this number times "
      "verletRebuildFrequency."};

  /**
   * fastParticlesThrow
   */
  MDFlexOption<bool, __LINE__> fastParticlesThrow{false, "fastParticlesThrow", false,
                                                  "Decide if particles that move farther than skin/2/rebuildFrequency "
                                                  "will throw an exception during the position update or not."};
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
   * loadBalancingInterval
   */
  MDFlexOption<unsigned int, __LINE__> loadBalancingInterval{
      100, "load-balancing-interval", true, "Defines the iteration interval at which load balancing should occur."};

  /**
   * subdivideDimension
   */
  MDFlexOption<std::array<bool, 3>, 0> subdivideDimension{
      {true, true, true},
      "subdivide-dimension",
      true,
      "Indicates in which dimensions the global domain can be subdivided by the MPI decomposition"};

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
  MDFlexOption<FunctorOption, __LINE__> functorOption {
    // choose a reasonable default depending on what is available at compile time
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
    FunctorOption::lj12_6_AVX,
#elif defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
    FunctorOption::lj12_6_SVE,
#else
    FunctorOption::lj12_6,
#endif
        "functor", true,
        "Force functor to use. Possible Values: (lennard-jones "
        "lennard-jones-AVX lennard-jones-SVE lennard-jones-globals)"
  };
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
   * Boundary types.
   */
  MDFlexOption<std::array<options::BoundaryTypeOption, 3>, __LINE__> boundaryOption{
      {options::BoundaryTypeOption::periodic, options::BoundaryTypeOption::periodic,
       options::BoundaryTypeOption::periodic},
      "boundary-type",
      true,
      "Boundary condition types for each of the three dimensions. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(options::BoundaryTypeOption::getAllOptions(), " ", {"(", ")"}) +
          " Default: {periodic, periodic, periodic}"};
  /**
   * Omit the creation of a config file at the end of the Simulation.
   * This starts with a "not" such that it can be used as a flag with a sane default.
   */
  MDFlexOption<bool, __LINE__> dontCreateEndConfig{
      true, "no-end-config", false, "Set to omit the creation of a yaml file at the end of a simulation."};
  /**
   * Omit the output of the progress bar and progress information. This might be useful if this upsets your output.
   */
  MDFlexOption<bool, __LINE__> dontShowProgressBar{false, "no-progress-bar", false,
                                                   "Set to omit printing the progress bar."};
  /**
   * deltaT
   */
  MDFlexOption<double, __LINE__> deltaT{0.001, "deltaT", true,
                                        "Length of a timestep. Set to 0 to deactivate time integration."};

  /**
   * pauseSimulationDuringTuning
   */
  MDFlexOption<bool, __LINE__> pauseSimulationDuringTuning{false, "pause-simulation-during-tuning", false,
                                                           "Pauses the update of the simulation during tuning phases."};
  /**
   * sortingThreshold
   * This value is used in traversal that use the CellFunctor. If the sum of the number of particles in two cells is
   * greater or equal to that value, the CellFunctor creates a sorted view of the particles to avoid unnecessary
   * distance checks.
   */
  MDFlexOption<size_t, __LINE__> sortingThreshold{
      8, "sorting-threshold", true,
      "Threshold for traversals that use the CellFunctor to start sorting. If the sum of the number of particles in "
      "two cells is greater or equal to that value, the CellFunctor creates a sorted view of the particles to avoid "
      "unnecessary distance checks."};

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
      "Scenario generator. Possible Values: (grid uniform gaussian sphere closestPacking) Default: grid"};

  // Site Type Generation
  /**
   * siteStr
   */
  static inline const char *const siteStr{"Sites"};
  /**
   * epsilonMap
   */
  MDFlexOption<std::map<unsigned long, double>, 0> epsilonMap{
      {{0ul, 1.}}, "epsilon", true, "Mapping from site type to an epsilon value."};
  /**
   * sigmaMap
   */
  MDFlexOption<std::map<unsigned long, double>, 0> sigmaMap{
      {{0ul, 1.}}, "sigma", true, "Mapping from site type to a sigma value."};
  /**
   * massMap
   */
  MDFlexOption<std::map<unsigned long, double>, 0> massMap{
      {{0ul, 1.}}, "mass", true, "Mapping from site type to a mass value."};
  // Molecule Type Generation
  // Strings for parsing yaml files.
  /**
   * moleculesStr
   */
  static inline const char *const moleculesStr{"Molecules"};
  /**
   * moleculeToSiteIdStr
   */
  static inline const char *const moleculeToSiteIdStr{"site-types"};
  /**
   * moleculeToSitePosStr
   */
  static inline const char *const moleculeToSitePosStr{"relative-site-positions"};
  /**
   * momentOfInertiaStr
   */
  static inline const char *const momentOfInertiaStr{"moment-of-inertia"};
  // Maps where the molecule type information is actually stored
  /**
   * molToSiteIdMap
   */
  std::map<unsigned long, std::vector<unsigned long>> molToSiteIdMap{};
  /**
   * molToSitePosMap
   */
  std::map<unsigned long, std::vector<std::array<double, 3>>> molToSitePosMap{};
  /**
   * momentOfInertiaMap
   */
  std::map<unsigned long, std::array<double, 3>> momentOfInertiaMap{};
  // Object Generation:
  /**
   * objectsStr
   */
  static inline const char *const objectsStr{"Objects"};
  /**
   * particleTypeStr. A md-flex mode blind string name for the particle type of the object's particles. E.g. this is
   * site-type-id in a single-site simulation and molecule-type-id in a multi-site simulation.
   */
  static inline const char *const particleTypeStr{"particle-type-id"};
  /**
   * bottomLeftBackCornerStr
   */
  static inline const char *const bottomLeftBackCornerStr{"bottomLeftCorner"};
  /**
   * velocityStr
   */
  static inline const char *const velocityStr{"velocity"};
  /**
   * particlesPerObjectStr
   */
  static inline const char *const particlesPerObjectStr{"numberOfParticles"};
  /**
   * cubeGridObjectsStr
   */
  static inline const char *const cubeGridObjectsStr{"CubeGrid"};
  /**
   * cubeGridObjects
   */
  std::vector<CubeGrid> cubeGridObjects{};
  /**
   * cubeGaussObjectsStr
   */
  static inline const char *const cubeGaussObjectsStr{"CubeGauss"};
  /**
   * cubeGaussObjects
   */
  std::vector<CubeGauss> cubeGaussObjects{};
  /**
   * cubeUniformObjectsStr
   */
  static inline const char *const cubeUniformObjectsStr{"CubeUniform"};
  /**
   * cubeUniformObjects
   */
  std::vector<CubeUniform> cubeUniformObjects{};
  /**
   * sphereObjectsStr
   */
  static inline const char *const sphereObjectsStr{"Sphere"};
  /**
   * sphereCenterStr
   */
  static inline const char *const sphereCenterStr{"center"};
  /**
   * sphereRadiusStr
   */
  static inline const char *const sphereRadiusStr{"radius"};
  /**
   * sphereObjects
   */
  std::vector<Sphere> sphereObjects{};
  /**
   * cubeClosestPackedObjects
   */
  std::vector<CubeClosestPacked> cubeClosestPackedObjects{};
  /**
   * cubeClosestPackedObjectsStr
   */
  static inline const char *const cubeClosestPackedObjectsStr{"CubeClosestPacked"};

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
  MDFlexOption<double, __LINE__> initTemperature{0., "initialTemperature", true,
                                                 "Thermostat option. Initial temperature of the system."};
  /**
   * targetTemperature
   */
  MDFlexOption<double, __LINE__> targetTemperature{0., "targetTemperature", true,
                                                   "Thermostat option. Target temperature of the system."};
  /**
   * deltaTemp
   */
  MDFlexOption<double, __LINE__> deltaTemp{
      0., "deltaTemperature", true, "Thermostat option. Maximal temperature jump the thermostat is allowed to apply."};
  /**
   * thermostatInterval
   */
  MDFlexOption<size_t, __LINE__> thermostatInterval{
      0, "thermostatInterval", true,
      "Thermostat option. Number of Iterations between two applications of the thermostat."};
  /**
   * addBrownianMotion
   */
  MDFlexOption<bool, __LINE__> addBrownianMotion{
      true, "addBrownianMotion", true,
      "Thermostat option. Whether the particle velocities should be initialized using "
      "Brownian motion. Possible Values: (true false) Default: true"};

  /**
   * Global external force like e.g. gravity
   */
  MDFlexOption<std::array<double, 3>, __LINE__> globalForce{
      {0, 0, 0},
      "globalForce",
      true,
      "Global force applied on every particle. Useful to model e.g. gravity. Default: {0,0,0}"};

  /**
   * Convenience function testing if the global force contains only 0 entries.
   * @return false if any entry in globalForce.value is != 0.
   */
  [[nodiscard]] bool globalForceIsZero() const {
    bool isZero = true;
    for (auto gf : globalForce.value) {
      isZero &= gf == 0;
    }
    return isZero;
  }

  // Checkpoint Options
  /**
   * checkpointfile
   */
  MDFlexOption<std::string, __LINE__> checkpointfile{"", "checkpoint", true,
                                                     "Path to a .pvtu File to load as a checkpoint."};

  /**
   * loadBalancer
   */
  MDFlexOption<LoadBalancerOption, __LINE__> loadBalancer{
      LoadBalancerOption::invertedPressure, "load-balancer", true,
      "Defines which load balancing approach will be used with the adaptive grid decomposition. If ALL is chosen as "
      "load balancer, MD-Flexible uses ALL's TENSOR method. Possible Values: " +
          autopas::utils::ArrayUtils::to_string(LoadBalancerOption::getAllOptions(), " ", {"(", ")"})};

  /**
   * Whether to use the tuning logger or not.
   *
   * @see TuningStrategyLoggerProxy
   */
  MDFlexOption<bool, __LINE__> useTuningLogger{false, "use-tuning-logger", true,
                                               "If tuning information should be logged. Possible Values: (true false)"};

  /**
   * The suffix for files created by the tuning logger.
   */
  MDFlexOption<std::string, __LINE__> outputSuffix{"", "output-suffix", true,
                                                   "An identifier that is contained in the filename of all log files."};

  /**
   * valueOffset used for cli-output alignment
   */
  static constexpr int valueOffset{33};

  /**
   * Stores the particles generated based on the provided configuration file
   * These particles can be added to the respective autopas container,
   * but have to be converted to the respective particle type, first.
   */
  std::vector<ParticleType> particles{};

 private:
  /**
   * Stores the physical properties of the particles used in the an MDFlexSimulation
   */
  std::shared_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary;

  /**
   * Initializes the ParticlePropertiesLibrary
   */
  void initializeParticlePropertiesLibrary();

  /**
   * Initializes all particles present at the start of the simulation.
   */
  void initializeObjects();
};

/**
 * Stream insertion operator for MDFlexConfig.
 * @param os
 * @param config
 * @return
 */
inline std::ostream &operator<<(std::ostream &os, const MDFlexConfig &config) { return os << config.to_string(); }

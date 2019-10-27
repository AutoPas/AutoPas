/**
 * @file MDFlexConfig.h
 * @author F. Gratl
 * @date 10/18/19
 */

#pragma once

#include <map>
#include <set>

#include "Objects.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/utils/NumberSet.h"

/**
 * Class containing all necessary parameters for configuring a md-flexible simulation.
 */
class MDFlexConfig {
 public:
  /**
   * Print the content of the config to std::out.
   */
  void print();

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
  enum class FunctorOption { lj12_6, lj12_6_AVX };

  /**
   * Choice of the particle generators specified in the command line
   */
  enum class GeneratorOption { grid, uniform, gaussian, sphere };

  static inline const char *yamlFilenameStr = "yaml-filename";
  std::string yamlFilename;  // default configuration = CubeGrid

  // AutoPas options:
  static inline const char *containerOptionsStr = "container";
  std::set<autopas::ContainerOption> containerOptions = autopas::allContainerOptions;
  static inline const char *dataLayoutOptionsStr = "data-layout";
  std::set<autopas::DataLayoutOption> dataLayoutOptions = autopas::allDataLayoutOptions;
  static inline const char *selectorStrategyStr = "selector-strategy";
  autopas::SelectorStrategyOption selectorStrategy = autopas::SelectorStrategyOption::fastestAbs;
  static inline const char *traversalOptionsStr = "traversal";
  std::set<autopas::TraversalOption> traversalOptions = autopas::allTraversalOptions;
  static inline const char *newton3OptionsStr = "newton3";
  std::set<autopas::Newton3Option> newton3Options = autopas::allNewton3Options;
  static inline const char *cellSizeFactorsStr = "cell-size";
  std::shared_ptr<autopas::NumberSet<double>> cellSizeFactors =
      std::make_shared<autopas::NumberSetFinite<double>>(std::set<double>{1.});
  static inline const char *logFileNameStr = "log-file";
  std::string logFileName;
  static inline const char *logLevelStr = "log-level";
  spdlog::level::level_enum logLevel = spdlog::level::info;
  static inline const char *tuningStrategyOptionsStr = "tuning-strategy";
  autopas::TuningStrategyOption tuningStrategyOption = autopas::TuningStrategyOption::fullSearch;
  static inline const char *tuningIntervalStr = "tuning-interval";
  unsigned int tuningInterval = 100;
  static inline const char *tuningSamplesStr = "tuning-samples";
  unsigned int tuningSamples = 3;
  static inline const char *tuningMaxEvidenceStr = "tuning-max-evidence";
  unsigned int tuningMaxEvidence = 10;
  static inline const char *vtkFileNameStr = "vtk-filename";
  std::string vtkFileName;
  static inline const char *vtkWriteFrequencyStr = "vtk-write-frequency";
  size_t vtkWriteFrequency = 100;
  static inline const char *verletRebuildFrequencyStr = "verlet-rebuild-frequency";
  unsigned int verletRebuildFrequency = 1;
  static inline const char *verletSkinRadiusStr = "verlet-skin-radius";
  double verletSkinRadius = .2;
  static inline const char *boxMinStr = "box-min";
  std::array<double, 3> boxMin = {0, 0, 0};
  static inline const char *boxMaxStr = "box-max";
  std::array<double, 3> boxMax = {5, 5, 5};

  // Simulation Options:
  static inline const char *cutoffStr = "cutoff";
  double cutoff = 1.;
  static inline const char *functorOptionStr = "functor";
  FunctorOption functorOption = FunctorOption::lj12_6;
  static inline const char *iterationsStr = "iterations";
  size_t iterations = 10;
  static inline const char *periodicStr = "periodic-boundaries";
  bool periodic = true;
  static inline const char *measureFlopsStr = "no-flops";
  bool measureFlops = true;
  static inline const char *deltaTStr = "deltaT";
  double deltaT = 0.001;
  static inline const char *epsilonStr = "particle-epsilon";
  std::map<unsigned long, double> epsilonMap = {{0, 1}};
  static inline const char *sigmaStr = "particle-sigma";
  std::map<unsigned long, double> sigmaMap = {{0, 1}};
  static inline const char *massStr = "particle-mass";
  std::map<unsigned long, double> massMap = {{0, 1}};

  // Options for additional Object Generation on command line
  static inline const char *boxLengthStr = "box-length";
  double boxLength = 10;
  static inline const char *distributionMeanStr = "distribution-mean";
  double distributionMean = 5.;
  static inline const char *distributionStdDevStr = "distribution-stddeviation";
  double distributionStdDev = 2.;
  static inline const char *particlesPerDimStr = "particles-per-dimension";
  size_t particlesPerDim = 10;
  static inline const char *particlesTotalStr = "particles-total";
  size_t particlesTotal = 1000;
  static inline const char *particlesSpacingStr = "particle-spacing";
  double particleSpacing = .5;
  static inline const char *generatorOptionStr = "particle-generator";
  GeneratorOption generatorOption = GeneratorOption::grid;

  // Object Generation:
  static inline const char *objectsStr = "Objects";
  static inline const char *bottomLeftBackCornerStr = "bottomLeftCorner";
  static inline const char *velocityStr = "velocity";
  static inline const char *particleTypeStr = "particle-type";
  static inline const char *particlesPerObjectStr = "numberOfParticles";
  static inline const char *cubeGridObjectsStr = "CubeGrid";
  std::vector<CubeGrid> cubeGridObjects = {};
  static inline const char *cubeGaussObjectsStr = "CubeGauss";
  std::vector<CubeGauss> cubeGaussObjects = {};
  static inline const char *cubeUniformObjectsStr = "CubeUniform";
  std::vector<CubeUniform> cubeUniformObjects = {};
  static inline const char *sphereObjectsStr = "Sphere";
  static inline const char *sphereCenterStr = "center";
  static inline const char *sphereRadiusStr = "radius";
  std::vector<Sphere> sphereObjects = {};

  // Thermostat Options
  static inline const char *thermostatStr = "thermostat";
  bool useThermostat = false;
  static inline const char *initTemperatureStr = "initialTemperature";
  double initTemperature = 0.;
  static inline const char *targetTemperatureStr = "targetTemperature";
  double targetTemperature = 0.;
  static inline const char *deltaTempStr = "deltaTemperature";
  double deltaTemp = 0.;
  static inline const char *thermostatIntervalStr = "thermostatInterval";
  size_t thermostatInterval = 0.;
  static inline const char *useCurrentTempForBrownianMotionStr = "useCurrentTempForBrownianMotion";
  bool useCurrentTempForBrownianMotion = true;

  // used for cli-output alignment
  static constexpr size_t valueOffset = 32;
};

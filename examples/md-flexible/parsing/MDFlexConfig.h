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
  enum class GeneratorOption { empty, grid, uniform, gaussian };

  std::string yamlFilename;  // default configuration = CubeGrid

  // AutoPas options:
  std::set<autopas::ContainerOption> containerOptions = autopas::allContainerOptions;
  std::set<autopas::DataLayoutOption> dataLayoutOptions = autopas::allDataLayoutOptions;
  autopas::SelectorStrategyOption selectorStrategy = autopas::SelectorStrategyOption::fastestAbs;
  std::set<autopas::TraversalOption> traversalOptions = autopas::allTraversalOptions;
  autopas::TuningStrategyOption tuningStrategyOption = autopas::TuningStrategyOption::fullSearch;
  std::set<autopas::Newton3Option> newton3Options = autopas::allNewton3Options;
  std::shared_ptr<autopas::NumberSet<double>> cellSizeFactors =
      std::make_shared<autopas::NumberSetFinite<double>>(std::set<double>{1.});
  spdlog::level::level_enum logLevel = spdlog::level::info;
  unsigned int tuningInterval = 100;
  unsigned int tuningSamples = 3;
  unsigned int tuningMaxEvidence = 10;
  std::string VTKFileName;
  size_t vtkWriteFrequency = 100;
  std::string logFileName;
  unsigned int verletRebuildFrequency = 20;
  double verletSkinRadius = .2;
  std::array<double, 3> boxMin = {0,0,0};
  std::array<double, 3> boxMax = {1,1,1};

  // Simulation Options:
  double cutoff = 1.;
  FunctorOption functorOption = FunctorOption::lj12_6;
  size_t iterations = 10;
  bool periodic = false;
  bool measureFlops = true;
  double deltaT = 0.001;
  std::map<unsigned long, double> epsilonMap = {{0,1}};
  std::map<unsigned long, double> sigmaMap = {{0,1}};
  std::map<unsigned long, double> massMap = {{0,1}};

  // Options for additional Object Generation on command line
  double boxLength = 10;
  double distributionMean = 5.;
  double distributionStdDev = 2.;
  size_t particlesPerDim = 10;
  size_t defaultParticlesTotal = 1000;
  double particleSpacing = .5;
  GeneratorOption generatorOption = GeneratorOption::grid;

  // Object Generation:
  std::vector<CubeGrid> cubeGridObjects = {
      {{particlesPerDim, particlesPerDim, particlesPerDim}, particleSpacing, {0., 0., 0.}, {0., 0., 0.}, 0, 1, 1, 1}};
  std::vector<CubeGauss> cubeGaussObjects = {};
  std::vector<CubeUniform> cubeUniformObjects = {};
  std::vector<Sphere> sphereObjects = {};

  // Thermostat Options default to false
  bool thermostat = false;
  bool initializeThermostat = false;
  double initTemperature = 0.;
  size_t numberOfTimesteps = 0.;
  bool thermoTarget = false;
  double targetTemperature = 0.;
  double deltaTemp = 0.;

  // used for cli-output alignment
  static constexpr size_t valueOffset = 32;
};

#pragma once

#include <getopt.h>
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include "Objects.h"
#include "Thermostat.h"
#include "autopas/autopasIncludes.h"
#include "autopas/utils/NumberSet.h"
class YamlParser {
  /**
   * @file MDFlexParser.h
   * @date 23.02.2018
   * @author F. Gratl
   */
 public:
  /**
   * Choice of the functor
   */
  enum FunctorOption { lj12_6, lj12_6_AVX };

  /**
   * Choice of the particle generators specified in the command line
   */
  enum GeneratorOption { empty, grid, uniform, gaussian };

  /**Constructor für YAMl Parser:
   * */
  YamlParser() = default;
  /**Copy Contructor
   * */
  YamlParser(const YamlParser &parser) = default;
  /**Copy assignment operator
   * */
  YamlParser &operator=(const YamlParser &parser) = default;

  /**Getter for BoxMin for Autopas Object, needed in autopas initialization
   * @return BoxMin
   * */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const;
  /**Getter for BoxMax for Autopas Object, needed in autopas initialization
   * @return BoxMax
   * */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const;

  /** Parse the Input from the command line
   *  If a Yaml File is specified, it will be parsed first
   *  If additional Simulation Options come after the .yaml Argument, these Options override the .yaml ones
   * */
  bool parseInput(int argc, char **argv);

  /**Parses the Input for the simulation from the Yaml File
   * @param filename
   * */
  void parseYamlFile();

  /**Prints Configuration of Simulation:
   * */
  void printConfig();
  /**Calculates the total number of Particles generated
   * @return particlestotal
   * */
  size_t particlesTotal();

  /**Calculate the required Box for the AutoPas Object
   * */
  void calcAutopasBox();

  /**Set up the data Structure for Particle types with their properties to be read by the Simulation
   * */
  void addType(unsigned long typeId, double epsilon, double sigma, double mass);

  [[nodiscard]] const std::set<autopas::ContainerOption> &getContainerOptions() const;

  [[nodiscard]] const std::set<autopas::DataLayoutOption> &getDataLayoutOptions() const;

  [[nodiscard]] autopas::SelectorStrategyOption getSelectorStrategy() const;

  [[nodiscard]] const std::set<autopas::TraversalOption> &getTraversalOptions() const;

  [[nodiscard]] autopas::TuningStrategyOption getTuningStrategyOption() const;

  [[nodiscard]] const std::set<autopas::Newton3Option> &getNewton3Options() const;

  [[nodiscard]] const autopas::NumberSet<double> &getCellSizeFactors() const;

  [[nodiscard]] double getCutoff() const;

  [[nodiscard]] FunctorOption getFunctorOption() const;

  [[nodiscard]] size_t getIterations() const;

  [[nodiscard]] spdlog::level::level_enum getLogLevel() const;

  [[nodiscard]] bool getMeasureFlops() const;

  [[nodiscard]] unsigned int getTuningInterval() const;

  [[nodiscard]] unsigned int getTuningSamples() const;

  [[nodiscard]] unsigned int getTuningMaxEvidence() const;

  [[nodiscard]] const std::string &getVTKFileName() const;

  [[nodiscard]] const std::string &getLogFileName() const;

  [[nodiscard]] unsigned int getVerletRebuildFrequency() const;

  [[nodiscard]] double getVerletSkinRadius() const;

  [[nodiscard]] double getDeltaT() const;

  [[nodiscard]] const std::vector<CubeGrid> &getCubeGrid() const;

  [[nodiscard]] const std::vector<CubeGauss> &getCubeGauss() const;

  [[nodiscard]] const std::vector<CubeUniform> &getCubeUniform() const;

  [[nodiscard]] const std::vector<Sphere> &getSphere() const;

  [[nodiscard]] const std::map<unsigned long, double> &getEpsilonMap() const;

  [[nodiscard]] const std::map<unsigned long, double> &getSigmaMap() const;

  [[nodiscard]] const std::map<unsigned long, double> &getMassMap() const;

  void setFilename(const std::string &inputFilename);
  [[nodiscard]] size_t getVtkWriteFrequency() const;

  void setVtkWriteFrequency(size_t vtkWriteFrequency);
  void setVtkFileName(const std::string &vtkFileName);

    bool isThermostat() const;

    double getInitTemperature() const;

    size_t getNumberOfTimesteps() const;

    double getTargetTemperature() const;

    double getDeltaTemp() const;

    bool isThermoTarget() const;

private:
  static constexpr size_t valueOffset = 32;
  // defaults:
  std::string filename;  // default configuration = CubeGrid
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
  std::string VTKFileName = "";
  size_t vtkWriteFrequency = 100;
  std::string logFileName = "";
  unsigned int verletRebuildFrequency = 20;
  double verletSkinRadius = .2;
  std::array<double, 3> BoxMin = {};
  std::array<double, 3> BoxMax = {};

  // Simulation Options:
  double cutoff = 1.;
  FunctorOption functorOption = FunctorOption::lj12_6;
  size_t iterations = 10;
  bool periodic = true;
  bool measureFlops = true;
  double delta_t = 0.001;
  std::map<unsigned long, double> epsilonMap;
  std::map<unsigned long, double> sigmaMap;
  std::map<unsigned long, double> massMap;

  // Options for additional Object Generation on command line
  double boxLength = 10;
  double distributionMean = 5.;
  double distributionStdDev = 2.;
  size_t particlesPerDim = 10;
  size_t defaultParticlesTotal = 1000;
  double particleSpacing = .5;
  GeneratorOption generatorOption = GeneratorOption::empty;

  // Object Generation:
  std::vector<CubeGrid> CubeGridObjects = {};
  std::vector<CubeGauss> CubeGaussObjects = {};
  std::vector<CubeUniform> CubeUniformObjects = {};
  std::vector<Sphere> SphereObjects = {};

  //Thermostat Options default to false
  bool thermostat = false;
  double initTemperature=0.;
  size_t numberOfTimesteps=0.;
  bool ThermoTarget = false;
  double targetTemperature=0.;
  double delta_temp=0.;

};
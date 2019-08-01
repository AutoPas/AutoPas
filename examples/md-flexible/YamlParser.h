#pragma once

#include <getopt.h>
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include "Objects.h"
#include "autopas/autopasIncludes.h"
#include "autopas/utils/NumberSet.h"
class YamlParser {
  /**
   * @file MDFlexParser.h
   * @date 23.02.2018
   * @author F. Gratl
   */
 public:
  enum FunctorOption { lj12_6, lj12_6_AVX };

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

  [[nodiscard]] const std::string &getWriteVtk() const;

  [[nodiscard]] const std::string &getLogFileName() const;

  [[nodiscard]] unsigned int getVerletRebuildFrequency() const;

  [[nodiscard]] double getVerletSkinRadius() const;

  [[nodiscard]] double getDeltaT() const;

  [[nodiscard]] const std::vector<CubeGrid> &getCubeGrid() const;

  [[nodiscard]] const std::vector<CubeGauss> &getCubeGauss() const;

  [[nodiscard]] const std::vector<CubeUniform> &getCubeUniform() const;

  [[nodiscard]] const std::vector<Sphere> &getSphere() const;

  const std::map<unsigned long, double> &getEpsilonMap() const;

  const std::map<unsigned long, double> &getSigmaMap() const;

  const std::map<unsigned long, double> &getMassMap() const;

 private:
 public:
  void setFilename(const std::string &inputFilename);

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
  std::string writeVTK = "";
  std::string logFileName = "";
  unsigned int verletRebuildFrequency = 5;
  double verletSkinRadius = .2;
  std::array<double, 3> BoxMin = {0., 0., 0.};
  std::array<double, 3> BoxMax = {10., 10., 10.};

  // Simulation Options:
  double cutoff = 1.;
  FunctorOption functorOption = FunctorOption::lj12_6;
  size_t iterations = 10;
  bool measureFlops = true;
  double delta_t = 0.001;
  std::map<unsigned long, double> epsilonMap;
  std::map<unsigned long, double> sigmaMap;
  std::map<unsigned long, double> massMap;

  // Object Generation:
  std::vector<CubeGrid> CubeGridObjects = {CubeGrid({10, 10, 10}, 1., {0., 0., 0.}, {5., 5., 5.}, 0, 1.0, 1.0, 1.0)};
  std::vector<CubeGauss> CubeGaussObjects = {};
  std::vector<CubeUniform> CubeUniformObjects = {};
  std::vector<Sphere> SphereObjects = {};
};
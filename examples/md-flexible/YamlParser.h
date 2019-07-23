#pragma once

#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include "Objects.h"
#include "autopas/autopasIncludes.h"
#include "autopas/utils/NumberSet.h"
#include <limits>
#include <array>
#include <algorithm>
using namespace std;
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

    const array<double, 3> &getBoxMin() const;

    const array<double, 3> &getBoxMax() const;

    /**Parses the Input for the simulation
     * @param filename
     * */
  void parseInput(std::string &filename);

  /**Prints Configuration of Simulation:
   * */
   //@todo output besser strukturieren(object generation am ende, ...)
  void printConfig();
  /**Calculates the total number of Particles generated
   * @return particlestotal
   * */
  size_t particlesTotal();

  /**Calculate the required Box for the AutoPas Object
   * */
  void calcAutopasBox();

  const set<ContainerOption> &getContainerOptions() const;

  const set<DataLayoutOption> &getDataLayoutOptions() const;

  SelectorStrategyOption getSelectorStrategy() const;

  const set<TraversalOption> &getTraversalOptions() const;

  TuningStrategyOption getTuningStrategyOption() const;

  const set<Newton3Option> &getNewton3Options() const;

  const NumberSet<double> &getCellSizeFactors() const;

  double getCutoff() const;

  FunctorOption getFunctorOption() const;

  size_t getIterations() const;

  spdlog::level::level_enum getLogLevel() const;

  bool getMeasureFlops() const;

  unsigned int getTuningInterval() const;

  unsigned int getTuningSamples() const;

  unsigned int getTuningMaxEvidence() const;

  const string &getWriteVtk() const;

  const string &getLogFileName() const;

  unsigned int getVerletRebuildFrequency() const;

  double getVerletSkinRadius() const;

  double getEpsilon() const;

  double getSigma() const;

  double getDeltaT() const;

  double getMass() const;

  const vector<CubeGrid> &getCubeGrid() const;

  const vector<CubeGauss> &getCubeGauss() const;

  const vector<CubeUniform> &getCubeUniform() const;

  const vector<Sphere> &getSphere() const;

 private:
  static constexpr size_t valueOffset = 32;
  // defaults:

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
  std::array<double,3> BoxMin={0.,0.,0.};
  std::array<double,3> BoxMax={0.,0.,0.};
  // Simulation Options:
  double cutoff = 1.;
  FunctorOption functorOption = FunctorOption::lj12_6;
  size_t iterations = 10;
  bool measureFlops = true;
  double epsilon = 5.0;
  double sigma = 1.0;
  double delta_t = 0.001;
  double mass = 1.0;
  // Object Generation:
  std::vector<CubeGrid> CubeGridObjects = {};
  std::vector<CubeGauss> CubeGaussObjects = {};
  std::vector<CubeUniform> CubeUniformObjects = {};
  std::vector<Sphere> SphereObjects = {};
};
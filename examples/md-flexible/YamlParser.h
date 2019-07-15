#pragma once

#include <yaml/yaml.h>
#include "autopas/autopasIncludes.h"
using namespace std;
class YamlParser {
    /**
 * @file MDFlexParser.h
 * @date 23.02.2018
 * @author F. Gratl
 */
public:
    enum FunctorOption { lj12_6, lj12_6_AVX };
    enum GeneratorOption { grid, uniform, gaussian };
    /**Constructor für YAMl Parser:
     * */
    YamlParser(std::string filename);
    /**Destructor für YAML Parser
     * */
    YAMLParser() = default;

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

    // Simulation Options:
    double boxLength = -1;
    double cutoff = 1.;
    double distributionMean = 5.;
    double distributionStdDev = 2.;
    FunctorOption functorOption = FunctorOption::lj12_6;
    GeneratorOption generatorOption = GeneratorOption::grid;
    size_t iterations = 10;
    spdlog::level::level_enum logLevel = spdlog::level::info;
    bool measureFlops = true;
    size_t particlesPerDim = 20;
    size_t particlesTotal = 1000;
    double particleSpacing = .4;
    unsigned int tuningInterval = 100;
    unsigned int tuningSamples = 3;
    unsigned int tuningMaxEvidence = 10;
    std::string writeVTK = "";
    std::string logFileName = "";
    unsigned int verletRebuildFrequency = 5;
    double verletSkinRadius = .2;
    double epsilon = 5.0;
    double sigma = 1.0;
    double delta_t = 0.001;
    double mass = 1.0;
    //Object Generation Options:




};




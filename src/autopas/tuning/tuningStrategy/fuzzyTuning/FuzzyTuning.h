/**
 * @file FuzzyTuning.h
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#pragma once

#include <string>

#include "OutputMapper.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyControlSystem.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/LinguisticVariable.h"

namespace autopas {

using namespace autopas::FuzzyLogic;

/**
 * A tuning strategy that uses fuzzy logic to make predictions about the performance of configurations.
 * The goal of this tuning strategy is to allow expert users to encode their knowledge about the performance of
 * different configurations in a rule file and use this knowledge to let a FuzzyControlSystem combine all this
 * information to make predictions about the optimal configuration for the current simulation state.
 *
 * Similar to the RuleBasedTuning strategy, the "knowledge" is encoded in a rule file consisting of 4 parts:
 * - @see FuzzyControlSettings: Defines the settings of all FuzzyControlSystems.
 * - @see LinguisticVariable: Define the input and output variables of the fuzzy control system. Can be chosen mostly
 *   arbitrarily the only limitation is that currently all input variables of the Systems must be defined over
 * components of the @see LiveInfo struct. All the currently implemented membership functions can be looked up at @see
 * FuzzySetFactory.
 * - @see OutputMapper: Maps the output of the fuzzy control systems to a configuration suitable for autopas. This is
 *   necessary because the fuzzy control systems are functions $f: R^n -> R$ and we need to map the output value to a
 *   configuration. At the moment, this is done by specifying a bunch of possible @see Configuration at specific output
 *   values. The FuzzyControlSystem will then choose the configuration lying closest to the predicted output value.
 * - @see FuzzyRule: Define the rules that the fuzzy control system should follow.
 *
 * All those options are defined by the user in a domain specific language, formally described in FuzzyLanguage.g4,
 * which is parsed with the help of antlr4.
 * The rules are loaded once at the start of the simulation.
 *
 *
 * <b>Summary of the Fuzzy Language:</b>
 *
 * This is a small example of how a typical rule file could look like:
 * ```
 * # Define the settings of the fuzzy control system
 * FuzzySystemSettings:
 *      defuzzificationMethod: "MeanOfMaximum"
 *
 * # Define all of the linguistic variables together with their linguistic terms
 * FuzzyVariable: domain: "homogeneity" range: (-0.009, 0.1486)
 *      "lower than 0.049":     SigmoidFinite(0.0914, 0.049, 0.0065)
 *      "lower than 0.041":     SigmoidFinite(0.0834, 0.041, -0.001)
 *      "higher than 0.049":    SigmoidFinite(0.0065, 0.049, 0.0914)
 *      "higher than 0.041":    SigmoidFinite(-0.001, 0.041, 0.0834)
 *
 * FuzzyVariable: domain: "threadCount" range: (-19.938, 48.938)
 *      "lower than 18.0":      SigmoidFinite(38.938, 18.0,  -2.938)
 *      "lower than 26.0":      SigmoidFinite(46.938, 26.0,   5.061)
 *      "lower than 8.0":       SigmoidFinite(28.938,  8.0, -12.938)
 *      "higher than 18.0":     SigmoidFinite(-2.938, 18.0,  38.938)
 *      "higher than 26.0":     SigmoidFinite(5.0617, 26.0,  46.938)
 *      "higher than 8.0":      SigmoidFinite(-12.93,  8.0,  28.938)
 *
 * FuzzyVariable: domain: "particlesPerCellStdDev" range: (-0.017, 0.072)
 *      "lower than 0.013":     SigmoidFinite(0.0639, 0.038,  0.012)
 *      "lower than 0.014":     SigmoidFinite(0.0399, 0.014, -0.011)
 *      "higher than 0.013":    SigmoidFinite(0.012,  0.013,  0.0639)
 *      "higher than 0.014":    SigmoidFinite(-0.011, 0.014,  0.0399)
 *
 * FuzzyVariable: domain: "Newton 3" range: (0, 1)
 *       "disabled, enabled":   Gaussian(0.3333, 0.1667)
 *       "enabled":             Gaussian(0.6667, 0.1667)
 *
 * # Define how the result of the output variables should be interpreted in the context of autopas
 * OutputMapping:
 *  "Newton 3":
 *      0.333 => [newton3="disabled"], [newton3="enabled"]
 *      0.666 => [newton3="enabled"]
 *
 * # Define a bunch of rules connecting the input variables to the output variables
 * if ("threadCount" == "lower than 18.0") && ("threadCount" == "higher than 8.0")
 *     && ("homogeneity" == "lower than 0.041")
 *   then ("Newton 3" == "enabled")
 * if ("threadCount" == "higher than 26.0") && ("particlesPerCellStdDev" == "lower than 0.013")
 *   then ("Newton 3" == "disabled, enabled")
 * ...
 * ```
 *
 * A larger example file is stored in /examples/md-flexible/input/fuzzyRules.frule.
 *
 *
 * Due to the compilation cost of ANTLR and issues with compiling the bundled dependency uuid on some machines, this
 * tuning strategy can be disabled with the CMake option AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING=OFF.
 */
class FuzzyTuning : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the FuzzyTuning strategy.
   * @param fuzzyRuleFileName The path of the fuzzy rule file.
   */
  explicit FuzzyTuning(std::string fuzzyRuleFileName);

  [[nodiscard]] TuningStrategyOption getOptionType() const override;

  bool needsLiveInfo() const override;

  void receiveLiveInfo(const LiveInfo &value) override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  /**
   * Getter for the fuzzy control settings.
   * @return The fuzzy control settings.
   */
  [[nodiscard]] std::shared_ptr<FuzzyControlSettings> getFuzzyControlSettings() const;

  /**
   * Getter for the fuzzy control systems.
   * @return The fuzzy control systems.
   */
  [[nodiscard]] const std::map<std::string, std::shared_ptr<FuzzyControlSystem>> &getFuzzyControlSystems() const;

  /**
   * Getter for the output mappings.
   * @return The output mappings.
   */
  [[nodiscard]] const std::map<std::string, std::shared_ptr<OutputMapper>> &getOutputMappings() const;

 private:
  /**
   * Parses the given fuzzy rule file and returns a map of the form {consequent_domain: FuzzyControlSystem}.
   * @param fuzzyRuleFilename The name of the fuzzy rule file.
   * @return A tuple containing the linguistic variables, the output mappings and the fuzzy control systems.
   */
  [[nodiscard]] static std::tuple<
      std::shared_ptr<FuzzyControlSettings>, std::vector<std::shared_ptr<LinguisticVariable>>,
      std::map<std::string, std::shared_ptr<OutputMapper>>, std::map<std::string, std::shared_ptr<FuzzyControlSystem>>>
  parse(const std::string &fuzzyRuleFilename);

  /**
   * Interprets the output of the fuzzy control systems as individual predictors and sets the configuration queue to the
   * cartesian product of all the predictions.
   *
   * This mode is intended for rule-bases designed for predicting individual parameters of the configuration and
   * combines the predictions of all the fuzzy control systems to a configuration queue adhering to the constraints of
   * the systems.
   *
   * Lets say we have 3 FuzzyControlSystems which resulted in the the following configuration patterns:
   * System 1:
   *   [container="LinkedCells", dataLayout="SoA"], [container="VerletClusterLists", dataLayout="SoA"],
   *   [container="VerletListsCells", dataLayout="AoS"]
   * System 2:
   *   [traversal="lc_c04"], [traversal="lc_c04_HCP"], [traversal="lc_c08"]
   * System 3:
   *   [newton3="enabled"]
   *
   * This method combines all the configurations similar to a cartesian product and would result in a configuration
   * queue with 3*3*1=9 configurations. However, as there are some invalid parameter combinations the resulting queue is
   * generally smaller. In every case however, all elements of the resulting queue fulfill the constraints of *all* the
   * fuzzy control systems.
   *
   * The resulting configuration queue could therefore include configurations like:
   * [ {Container: LinkedCells , CellSizeFactor: 1.000000 , Traversal: lc_c08 , Load Estimator: none , Data Layout: SoA
   * , Newton 3: enabled}, ... ]
   *
   * If there is no constraint about certain parameters all possible values for those parameters are allowed.
   *
   * @param configQueue The queue of configurations to be tested.
   */
  void updateQueueInterpretOutputAsIndividualSystems(std::vector<Configuration> &configQueue);

  /**
   * Interprets the output of the fuzzy control systems as a bunch of suitability values for each configuration and
   * updates the configuration queue to contain only the configurations that are most suitable according to the
   * predictions.
   *
   * This mode is intended for rule-bases designed to predict the suitability of **every** possible configuration.
   *
   * This mode assumes there is a LinguisticVariable for each configuration that can be applied. Such variables could
   * look like:
   *
   * FuzzyVariable: domain: "Suitability LinkedCells_AoS_lc_c01_disabled" range: (0, 1)
   *   "terrible":    Gaussian(0.125, 0.017)
   *   "medium":      Gaussian(0.5, 0.017)
   *   "excellent":   Gaussian(0.875, 0.017)
   *  ...
   *
   * The FuzzyControlSystem will then predict the suitability of each configuration and sort them according to the
   * predicted suitability (the higher the better). The resulting configuration-queue will then consist of the
   * configurations with the highest suitability values.
   *
   * It is possible to specify a band around the highest value such that more good configurations are included.
   * The width of this band can be chosen with the "suitabilityThreshold" parameter of the FuzzyControlSettings. (A
   * value of 0.1 means that all configurations within 10% of the highest suitability value are picked)
   *
   * @param configQueue The queue of configurations to be tested.
   */
  void updateQueueInterpretOutputAsSuitability(std::vector<Configuration> &configQueue);

  /**
   * The name of the fuzzy rule file.
   */
  std::string _fuzzyRuleFileName;

  /**
   * The current live info used to make predictions.
   */
  std::map<std::string, double> _currentLiveInfo;

  /**
   * The settings for the Fuzzy System
   */
  std::shared_ptr<FuzzyControlSettings> _fuzzyControlSettings;

  /**
   * The fuzzy control systems parsed from the fuzzy rule file.
   */
  std::map<std::string, std::shared_ptr<FuzzyControlSystem>> _fuzzyControlSystems;

  /**
   * The output mappings parsed from the fuzzy rule file.
   */
  std::map<std::string, std::shared_ptr<OutputMapper>> _outputMappings;
};

};  // namespace autopas

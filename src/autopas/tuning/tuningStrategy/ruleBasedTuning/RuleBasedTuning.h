/**
 * @file RuleBasedTuning.h
 * @author humig
 * @date 28.06.2021
 */

#pragma once

#include <fstream>
#include <list>
#include <unordered_map>

#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
#include "RuleBasedProgramParser.h"
#include "RuleBasedProgramTree.h"
#include "RuleVM.h"
#include "autopas/utils/WrapMPI.h"
#endif

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {
/**
 * A tuning strategy that uses information collected live from the domain to exclude configurations that knowingly
 * perform worse than other configurations in the current simulation state. The remaining configurations are tested
 * consecutively.
 *
 * This "knowledge" is encoded as rules in a rule file. The rules are defined by the user in a domain specific language,
 * formally described in RuleLanguage.g4, which is parsed by the RuleBasedProgramParser with the help of antlr4.
 * The rules are dynamically loaded and executed as a program for the RuleVM in the beginning of each tuning phase.
 *
 *
 * <b>Summary of the Rule Language:</b>
 *
 * The heart of this language are so called configuration orders. Here is one example:
 * ```
 * [container="LinkedCells", dataLayout="SoA", newton3="enabled"] >= [container="LinkedCells", dataLayout="AoS",
 * newton3="enabled"] with same traversal;
 * ```
 * This is a rule that states that *all* configurations that match the pattern to the left of `>=` are better than *all*
 * configurations that match the pattern to the right, provided they have the same `traversal`. The pattern to the left
 * e.g. matches all configuration with container type `LinkedCells`, data layout `SoA` and newton3 enabled.
 *
 * If we want to have this same rule with another container type, e.g. `LinkedCellsReferences`, we can combine these two
 * using a `define_list` statement and use that as container type in the rule:
 *
 * ```
 * define_list LinkedCellsContainer = "LinkedCells", "LinkedCellsReferences";
 * [container=LinkedCellsContainer, dataLayout="SoA", newton3="enabled"] >= [container="LinkedCells", dataLayout="AoS",
 * newton3="enabled"] with same traversal;
 * ```
 *
 * Because always applying the same configuration orders would not be very interesting, we can apply them conditionally.
 * For this, all live info collected in the LiveInfo class is available as a variable in the rule program. One example
 * is `numParticles`. Now let us define a rule which disables all `SoA` configurations if there are only very few
 * particles, because their `AoS` counterpart is better:
 *
 * ```
 * define lowNumParticlesThreshold = 500;
 * if numParticles < lowNumParticlesThreshold:
 *    [dataLayout="AoS"] >= [dataLayout="SoA"] with same container, newton3, traversal, loadEstimator;
 * endif
 * ```
 * At first, we define our threshold as a variable. Variables are always constant and have global visibility. Then, we
 * check if the threshold is not surpassed, and if yes, we apply the configuration order.
 *
 * A larger example file is stored in /examples/md-flexible/input/tuningRules.rule.
 *
 * Additionally, the class supports a so called "verify mode" where full search is performed and the given rules are
 * checked for correctness.
 *
 *
 * Due to the compilation cost of ANTLR and issues with compiling the bundled dependency uuid on some machines, this
 * tuning strategy can be disabled with the CMake option AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING=OFF.
 *
 */
class RuleBasedTuning : public TuningStrategyInterface {
 public:
  /**
   * A function type used to print errors found in verify mode.
   */
  using PrintTuningErrorFunType =
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
      std::function<void(const RuleSyntax::ConfigurationOrder &order, const Configuration &actualBetterConfig,
                         unsigned long betterRuntime, const Configuration &shouldBeBetterConfig,
                         unsigned long shouldBeBetterRuntime, const LiveInfo &liveInfo)>;
#else
      int;  // This is simply a dummy type and will never be used.
#endif

  /**
   * Constructs a RuleBasedTuning strategy.
   * @param searchSpace Set of all allowed configurations.
   * @param verifyModeEnabled If verify mode should be enabled. False by default.
   * @param ruleFileName The name of the file where the rules are stored.
   * @param tuningErrorPrinter The function to call in verify mode if errors are found.
   */
  explicit RuleBasedTuning(const std::set<Configuration> &searchSpace, bool verifyModeEnabled = false,
                           std::string ruleFileName = "tuningRules.rule",
                           PrintTuningErrorFunType tuningErrorPrinter = {});

  TuningStrategyOption getOptionType() const override;

  bool needsLiveInfo() const override;

  void receiveLiveInfo(const LiveInfo &info) override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  /**
   * @returns the total time that would have been skipped if verify mode was disabled and thus
   * some configurations would have been skipped due to the rules.
   */
  [[nodiscard]] long getLifetimeWouldHaveSkippedTuningTime() const;

  /**
   * @returns the total time spent trialling configurations suggested by this tuning strategy.
   */
  [[nodiscard]] long getLifetimeTuningTime() const;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

 private:
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
  /**
   * Creates a multi-line string representation of all rules contained in the given rule file.
   *
   * Works by printing all non-empty, non-comment lines in the file.
   *
   * @param filePath Path to the file whose rules should be printed.
   * @return String representation of all rules.
   */
  std::string rulesToString(const std::string &filePath) const;

  /**
   * Determines if the result of this configuration, compared with the results of other configurations tested previously
   * during this tuning phase, contradicts any of the applicable configuration orders and, if so, calls
   * tuningErrorPrinter.
   *
   * @param configuration
   */
  void verifyCurrentConfigTime(const Configuration &configuration) const;

  /**
   * Executes the rule file for the current simulation state.
   * Puts all known live info as "defines" (=definition of variables) in front of the program.
   */
  std::vector<RuleSyntax::ConfigurationOrder> applyRules(const std::vector<Configuration> &searchSpace);

  std::vector<RuleSyntax::ConfigurationOrder> _lastApplicableConfigurationOrders{};

  // The following member variables are only conditionally compiled to avoid warnings about unused variables.

  bool _verifyModeEnabled{};
  /**
   * Sum of all evidence since the last reset. This, times the number of samples, is approximately the time spent
   * trying out configurations.
   */
  long _tuningTime{0};
  long _wouldHaveSkippedTuningTime{0};
  /**
   * If the rules would immediately remove all options deactivate the whole strategy until the next reset.
   */
  bool _rulesTooHarsh{false};
  PrintTuningErrorFunType _tuningErrorPrinter{};
#endif

  std::set<Configuration> _searchSpace{};
  std::list<Configuration> _originalSearchSpace{};
  std::set<Configuration> _removedConfigurations{};

  std::string _ruleFileName{};

  std::unordered_map<Configuration, long, ConfigHash> _traversalTimes{};

  long _tuningTimeLifetime{0};
  long _wouldHaveSkippedTuningTimeLifetime{0};

  LiveInfo _currentLiveInfo{};
};
}  // namespace autopas

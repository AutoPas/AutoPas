/**
 * @file RuleBasedTuning.h
 * @author humig
 * @date 28.06.2021
 */

#pragma once

#include <fstream>
#include <iostream>
#include <list>
#include <unordered_set>
#include <variant>

#include "RuleBasedProgramParser.h"
#include "RuleBasedProgramTree.h"
#include "RuleVM.h"
#include "autopas/tuning/tuningStrategy/FullSearch.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {
/**
 * A tuning strategy that uses information collected live from the domain to exclude configurations that knowingly
 * perform worse than other configurations in the current simulation state. The remaining configurations are tested
 * using FullSearch.
 *
 * This "knowledge" is encoded as rules in a rule file. The rules are expected to be in their own little language,
 * formally described in RuleLanguage.g4. The rules are dynamically loaded and executed in the beginning of each tuning
 * phase.
 *
 * Here is a quick summary of this language:
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
 */
class RuleBasedTuning : public FullSearch {
 public:
  /**
   * A function type used to print errors found in verify mode.
   */
  using PrintTuningErrorFunType =
      std::function<void(const rule_syntax::ConfigurationOrder &order, const Configuration &actualBetterConfig,
                         unsigned long betterRuntime, const Configuration &shouldBeBetterConfig,
                         unsigned long shouldBeBetterRuntime, const LiveInfo &liveInfo)>;

  /**
   * Constructs a RuleBasedTuning strategy.
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param verifyModeEnabled If verify mode should be enabled. False by default.
   * @param ruleFileName The name of the file where the rules are stored.
   * @param tuningErrorPrinter The function to call in verify mode if errors are found.
   */
  RuleBasedTuning(const std::set<ContainerOption> &allowedContainerOptions,
                  const std::set<double> &allowedCellSizeFactors,
                  const std::set<TraversalOption> &allowedTraversalOptions,
                  const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                  const std::set<Newton3Option> &allowedNewton3Options, bool verifyModeEnabled = false,
                  std::string ruleFileName = "tuningRules.rule", PrintTuningErrorFunType tuningErrorPrinter = {})
      : FullSearch(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                   allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options),
        _originalSearchSpace(this->_searchSpace.begin(), this->_searchSpace.end()),
        _verifyModeEnabled(verifyModeEnabled),
        _ruleFileName(std::move(ruleFileName)),
        _tuningErrorPrinter(std::move(tuningErrorPrinter)) {}

  bool needsLiveInfo() const override { return true; }

  void receiveLiveInfo(const LiveInfo &info) override {
    _currentLiveInfo = info;
    _lastApplicableConfigurationOrders = applyRules();
  }

  void addEvidence(long time, size_t iteration) override {
    FullSearch::addEvidence(time, iteration);
    tuningTime += time;
    tuningTimeLifetime += time;
    if (_verifyModeEnabled) {
      verifyCurrentConfigTime();
      if (_removedConfigurations.find(this->getCurrentConfiguration()) != _removedConfigurations.end()) {
        wouldHaveSkippedTuningTime += time;
        wouldHaveSkippedTuningTimeLifetime += time;
      }
    }
  }

  void reset(size_t iteration) override {
    FullSearch::reset(iteration);

    if (_verifyModeEnabled and tuningTime > 0) {
      AutoPasLog(INFO, "Rules would have saved {} ns removing {}/{} configurations. ({}% of total tuning time)",
                 wouldHaveSkippedTuningTime, _removedConfigurations.size(), _originalSearchSpace.size(),
                 // TODO: This lambda could probably be replaced by some formatting parameters in the fmt string.
                 [&]() {
                   const auto percent =
                       ((static_cast<double>(wouldHaveSkippedTuningTime) / static_cast<double>(tuningTime)) * 100);
                   const auto percentRounded = std::round(percent * 100) / 100;
                   return percentRounded;
                 }());
    }

    tuningTime = 0;
    wouldHaveSkippedTuningTime = 0;
  }

  /**
   * @returns in verify mode the summed up time that would have been skipped if verify mode was disabled and
   * configurations would have been skipped due to the rules.
   */
  [[nodiscard]] auto getLifetimeWouldHaveSkippedTuningTime() const { return wouldHaveSkippedTuningTimeLifetime; }

  /**
   * @returns the summed up time of all configurations that have been tested by this tuning strategy.
   */
  [[nodiscard]] auto getLifetimeTuningTime() const { return tuningTimeLifetime; }

 private:
  /**
   * Goes through all applicable configuration orders and checks if the result of the current configuration contradicts
   * any rules when comparing with previously tested configurations in this tuning phase. If yes, calls
   * tuningErrorPrinter.
   */
  void verifyCurrentConfigTime() const {
    for (const auto &order : _lastApplicableConfigurationOrders) {
      bool shouldBeBetter;
      if (order.smaller.matches(*_currentConfig)) {
        shouldBeBetter = true;
      } else if (order.greater.matches(*_currentConfig)) {
        shouldBeBetter = false;
      } else {
        continue;
      }

      const auto currentConfigTime = _traversalTimes.at(*_currentConfig);
      const auto &comparePattern = shouldBeBetter ? order.greater : order.smaller;
      for (const auto &[otherConfig, time] : _traversalTimes) {
        bool error = false;
        if (comparePattern.matches(otherConfig) and order.haveEqualSameProperties(*_currentConfig, otherConfig)) {
          if ((shouldBeBetter and time > currentConfigTime) or (not shouldBeBetter and time < currentConfigTime)) {
            error = true;
          }
        }
        if (error) {
          if (shouldBeBetter) {
            _tuningErrorPrinter(order, *_currentConfig, currentConfigTime, otherConfig, time, _currentLiveInfo);
          } else {
            _tuningErrorPrinter(order, otherConfig, time, *_currentConfig, currentConfigTime, _currentLiveInfo);
          }
        }
      }
    }
  }

  /**
   * Executes the rule file for the current simulation state. Puts all known live info as defines in front of the
   * program.
   */
  std::vector<rule_syntax::ConfigurationOrder> applyRules() {
    AutoPasLog(DEBUG, _currentLiveInfo.toString());

    std::vector<RuleVM::MemoryCell> initialStack;
    std::vector<std::pair<std::string, rule_syntax::Define>> defines{};
    for (const auto &[name, value] : _currentLiveInfo.get()) {
      initialStack.emplace_back(value);
      defines.push_back({name, {name, std::make_shared<rule_syntax::Literal>(value)}});
    }

    std::ifstream ifs{_ruleFileName};
    const std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));

    rule_syntax::RuleBasedProgramParser parser{defines};
    auto [programTree, context] = parser.parse(content);

    const auto generatedProgram = programTree.generateCode(context);

    RuleVM vm{};
    const auto removePatterns = vm.execute(generatedProgram, initialStack);

    AutoPasLog(DEBUG, "Remove patterns (Count {}):", removePatterns.size());
    std::vector<ConfigurationPattern> toRemovePatterns{};
    std::vector<rule_syntax::ConfigurationOrder> applicableConfigurationOrders{};
    for (const auto &patternIdx : removePatterns) {
      const auto pattern = context.smallerConfigurationPatternByIndex(patternIdx);
      toRemovePatterns.push_back(pattern);
      AutoPasLog(DEBUG, "Remove {}", pattern.toString());

      applicableConfigurationOrders.push_back(context.getConfigurationOrders().at(patternIdx));
    }

    auto newSearchSpace = _originalSearchSpace;
    _removedConfigurations.clear();
    auto &removedConfigsLocal = _removedConfigurations;
    newSearchSpace.remove_if([&toRemovePatterns, &removedConfigsLocal](const Configuration &configuration) {
      const bool remove = std::any_of(toRemovePatterns.begin(), toRemovePatterns.end(),
                                      [&configuration](const auto &pattern) { return pattern.matches(configuration); });
      if (remove) {
        removedConfigsLocal.insert(configuration);
      }
      return remove;
    });
    AutoPasLog(DEBUG, "Rules remove {} out of {} configurations", _originalSearchSpace.size() - newSearchSpace.size(),
               _originalSearchSpace.size());
    if (not _verifyModeEnabled) {
      this->_searchSpace = {newSearchSpace.begin(), newSearchSpace.end()};
    }

    return applicableConfigurationOrders;
  }

  const std::list<Configuration> _originalSearchSpace;
  std::unordered_set<Configuration, ConfigHash> _removedConfigurations;
  std::vector<rule_syntax::ConfigurationOrder> _lastApplicableConfigurationOrders;
  bool _verifyModeEnabled;
  std::string _ruleFileName;

  long tuningTime = 0;
  long wouldHaveSkippedTuningTime = 0;
  long tuningTimeLifetime = 0;
  long wouldHaveSkippedTuningTimeLifetime = 0;

  PrintTuningErrorFunType _tuningErrorPrinter;

  LiveInfo _currentLiveInfo;
};
}  // namespace autopas

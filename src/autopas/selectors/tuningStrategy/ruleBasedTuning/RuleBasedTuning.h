/**
 * @file RuleBasedTuning.h
 * @author humig
 * @date 28.06.2021
 */

#pragma once

#include <fstream>
#include <iostream>
#include <list>
#include <variant>

#include "RuleBasedProgramParser.h"
#include "RuleBasedProgramTree.h"
#include "RuleVM.h"
#include "autopas/selectors/tuningStrategy/FullSearch.h"

namespace autopas {

class RuleBasedTuning : public FullSearch {
 public:
  RuleBasedTuning(const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
             const std::set<TraversalOption> &allowedTraversalOptions,
             const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
             const std::set<DataLayoutOption> &allowedDataLayoutOptions,
             const std::set<Newton3Option> &allowedNewton3Options, bool verifyModeEnabled = true)
      : FullSearch(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions, allowedLoadEstimatorOptions,
                   allowedDataLayoutOptions, allowedNewton3Options), _originalSearchSpace(this->_searchSpace.begin(),
                             this->_searchSpace.end()), _verifyModeEnabled(verifyModeEnabled)
      {}

  bool needsLiveInfo() const override {
    return true;
  }

  void receiveLiveInfo(LiveInfo info) override {
    _lastApplicableConfigurationOrders = applyRules(std::move(info));
  }

  void addEvidence(long time, size_t iteration) override {
    FullSearch::addEvidence(time, iteration);
    if(_verifyModeEnabled) {
      verifyCurrentConfigTime();
    }
  }

 private:
  void verifyCurrentConfigTime() const {
    for(const auto& order : _lastApplicableConfigurationOrders) {
      bool shouldBeBetter;
      if(order.smaller.matches(*_currentConfig)) {
        shouldBeBetter = true;
      } else if(order.greater.matches(*_currentConfig)) {
        shouldBeBetter = false;
      } else {
        continue;
      }

      const auto& comparePattern = shouldBeBetter ? order.greater : order.smaller;
      for(const auto& [otherConfig, time] : _traversalTimes) {
        bool error = false;
        if(comparePattern.matches(otherConfig) and order.haveEqualSameProperties(*_currentConfig, otherConfig)) {
          if((shouldBeBetter and time > _traversalTimes.at(*_currentConfig)) or
              (not shouldBeBetter and time < _traversalTimes.at(*_currentConfig))) {
            error = true;
          }
        }
        if (error) {
          std::cerr << "Error in ConfigurationOrder " << order.toString() << ":" << std::endl;
          std::cerr << "\tConfig " << _currentConfig->toShortString() << ": " << _traversalTimes.at(*_currentConfig)
                    << std::endl;
          std::cerr << "\tConfig " << otherConfig.toShortString() << ": " << time << std::endl << std::endl;
        }
      }
    }
  }

  std::vector<rule_syntax::ConfigurationOrder> applyRules(LiveInfo info) {
    std::cout << info.toString() << std::endl;

    std::vector<RuleVM::MemoryCell> initialStack;
    std::vector<std::pair<std::string, rule_syntax::Define>> defines{};
    for(const auto& [name, value] : info.get()) {
      initialStack.emplace_back(value);
      defines.push_back({name, {name, std::make_shared<rule_syntax::Literal>(value)}});
    }

    std::ifstream ifs("/home/tobias/tmp/myfile.txt");
    std::string content( (std::istreambuf_iterator<char>(ifs) ),
                         (std::istreambuf_iterator<char>()) );

    rule_syntax::RuleBasedProgramParser parser{defines};
    auto [programTree, context] = parser.parse(content);

    auto generatedProgram = programTree.generateCode(context);

    RuleVM vm{};
    auto removePatterns = vm.execute(generatedProgram, initialStack);

    std::cout << "Remove patterns (Count " << removePatterns.size() << "):" << std::endl;
    std::vector<ConfigurationPattern> toRemovePatterns{};
    std::vector<rule_syntax::ConfigurationOrder> applicableConfigurationOrders{};
    for(const auto& patternIdx : removePatterns) {
      auto pattern = context.smallerConfigurationPatternByIndex(patternIdx);
      toRemovePatterns.push_back(pattern);
      auto str = pattern.toString();
      std::cout << "Remove " << str << std::endl;

      applicableConfigurationOrders.push_back(context.getConfigurationOrders().at(patternIdx));
    }

    if(not _verifyModeEnabled) {
      auto newSearchSpace = _originalSearchSpace;
      newSearchSpace.remove_if(
          [&toRemovePatterns](const Configuration& configuration) {
            return std::any_of(toRemovePatterns.begin(), toRemovePatterns.end(),
                            [&configuration](const auto& pattern) {
                              return pattern.matches(configuration);
                            });
          });

      this->_searchSpace = {newSearchSpace.begin(), newSearchSpace.end()};
    }

    return applicableConfigurationOrders;
  }

 private:
  const std::list<Configuration> _originalSearchSpace;
  std::vector<rule_syntax::ConfigurationOrder> _lastApplicableConfigurationOrders;
  bool _verifyModeEnabled;
};
}  // namespace autopas

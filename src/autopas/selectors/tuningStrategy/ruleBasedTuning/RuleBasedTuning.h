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
             const std::set<Newton3Option> &allowedNewton3Options)
      : FullSearch(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions, allowedLoadEstimatorOptions,
                   allowedDataLayoutOptions, allowedNewton3Options), _originalSearchSpace(this->_searchSpace.begin(),
                             this->_searchSpace.end())
      {}

  bool needsLiveInfo() const override {
    return true;
  }

  void receiveLiveInfo(LiveInfo info) override {
    applyRules(std::move(info));
  }

 private:
  void applyRules(LiveInfo info) {
    std::cout << info.toString() << std::endl;

    auto newSearchSpace = _originalSearchSpace;

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
    for(const auto& patternIdx : removePatterns) {
      auto pattern = context.getConfigurationPatterns().at(patternIdx);
      toRemovePatterns.push_back(pattern);
      auto str = pattern.toString();
      std::cout << "Remove " << str << std::endl;
    }

    newSearchSpace.remove_if(
        [&toRemovePatterns](const Configuration& configuration) {
          bool res =
              std::any_of(toRemovePatterns.begin(), toRemovePatterns.end(),
                          [&configuration](const auto& pattern) {
                            return pattern.matches(configuration);
                         });
          if (res) {
            //std::cout << "!!!! Remove " << configuration.toString() << std::endl;
          }
          return res;
        });

    this->_searchSpace = {newSearchSpace.begin(), newSearchSpace.end()};
  }

 private:
  const std::list<Configuration> _originalSearchSpace;
};
}  // namespace autopas

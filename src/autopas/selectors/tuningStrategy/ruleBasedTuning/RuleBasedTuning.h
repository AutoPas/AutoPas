/**
 * @file RuleBasedTuning.h
 * @author humig
 * @date 28.06.2021
 */

#pragma once

#include <iostream>
#include <list>
#include <variant>

#include "autopas/selectors/tuningStrategy/FullSearch.h"

#include "RuleVM.h"
#include "RuleBasedProgramTree.h"
#include "RuleBasedProgramParser.h"

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
    auto cellsPerDim = utils::ArrayMath::mulScalar(info.domainSize, 1.0 / info.cutoff);
    auto numCells = static_cast<unsigned long>(std::ceil(cellsPerDim[0]) * std::ceil(cellsPerDim[1]) *
                                               std::ceil(cellsPerDim[2]));

    auto newSearchSpace = _originalSearchSpace;

    std::vector<ConfigurationPattern> patterns;
    patterns.push_back({{ContainerOption::linkedCells, ContainerOption::verletListsCells}});
    patterns.push_back({{ContainerOption::pairwiseVerletLists}});

    std::vector<RuleVM::MemoryCell> initialStack;
    initialStack.emplace_back(numCells);
    initialStack.emplace_back(info.numParticles);

    using I = RuleVM::Instruction;
    RuleVM::Program program{{
                                I{RuleVM::LOADA, 0ul},
                                I{RuleVM::LOADA, 1ul},
                                I{RuleVM::GREATER, 0ul},
                                I{RuleVM::CONDOUTPUTC, 0ul},
                                I{RuleVM::POP, 0ul},
                                I{RuleVM::LOADA, 1ul},
                                I{RuleVM::LOADC, 1000000ul},
                                I{RuleVM::GREATER, 0ul},
                                I{RuleVM::LOADA, 0ul},
                                I{RuleVM::LOADA, 1ul},
                                I{RuleVM::LESS, 0ul},
                                I{RuleVM::AND, 0ul},
                                I{RuleVM::CONDOUTPUTC, 1ul},
                                I{RuleVM::POP, 0ul},
                                I{RuleVM::HALT, 0ul}}, 100};


    using namespace rule_syntax;
    Define numParticlesDef{"numParticles", std::make_shared<Literal>(info.numParticles)};
    Define numCellsDef{"numCells", std::make_shared<Literal>(numCells)};
    CodeGenerationContext context{{}};
    context.addVariable(numCellsDef.variable);
    context.addVariable(numParticlesDef.variable);

    auto verletLists = std::make_shared<DefineList>("VerletListsContainer",
                                                    std::vector<std::shared_ptr<Literal>>{
                                                        std::make_shared<Literal>(ContainerOption::linkedCells),
                                                        std::make_shared<Literal>(ContainerOption::verletListsCells)});
    auto threshold = std::make_shared<Define>("LC_VL_THRESHOLD", std::make_shared<Literal>(2ul << 20));
    RuleBasedProgramTree programTree{
        {
            verletLists,
            threshold,
            std::make_shared<If>(
                std::make_shared<BinaryOperator>(BinaryOperator::AND,
                                                 std::make_shared<BinaryOperator>(BinaryOperator::GREATER,
                                                                                  std::make_shared<Variable>(&numParticlesDef),
                                                                                      std::make_shared<Variable>(&*threshold)),
                                                 std::make_shared<BinaryOperator>(BinaryOperator::LESS,
                                                                                  std::make_shared<Variable>(&numCellsDef),
                                                                                      std::make_shared<Variable>(&numParticlesDef))),
                std::vector<std::shared_ptr<Statement>>{std::make_shared<ConfigurationOrder>(patterns[0], patterns[0])})
    }};
    auto generatedProgram = programTree.generateCode(context);

    RuleVM vm{};
    auto removePatternsOld = vm.execute(program, initialStack);
    auto removePatterns = vm.execute(generatedProgram, initialStack);

    std::cout << "Old: ";
    for(const auto& val : removePatternsOld) {
      std::cout << val << " ";
    }
    std::cout << std::endl << "New: ";
    for(const auto& val : removePatterns) {
      std::cout << val;
    }
    std::cout << std::endl;

    std::vector<ConfigurationPattern> toRemovePatterns{};
    for(const auto& patternIdx : removePatterns) {
      toRemovePatterns.push_back(patterns.at(patternIdx));
    }

    newSearchSpace.remove_if([&toRemovePatterns](const Configuration& configuration) {
                    bool res = std::any_of(toRemovePatterns.begin(), toRemovePatterns.end(),
                         [&configuration](const auto& pattern) {
                            return pattern.matches(configuration);
                         });
                    if (res) {
                      std::cout << "!!!! Remove " << configuration.toString() << std::endl;
                    }
                    return res;
                  });

    /*if (numCells > info.numParticles) {
      newSearchSpace.remove_if([](const Configuration& configuration){
        return configuration.container != ContainerOption::verletClusterLists;
      });
      std::cout << "only allow cluster lists" << std::endl;
    } else if (info.numParticles > 1000000 and numCells < info.numParticles * 3) {
      std::cout << "remove pairwiseVerletLists" << std::endl;
      newSearchSpace.remove_if([](const Configuration& configuration){
        return configuration.container == ContainerOption::pairwiseVerletLists;
      });
    }*/



    this->_searchSpace = {newSearchSpace.begin(), newSearchSpace.end()};
  }

 private:
  const std::list<Configuration> _originalSearchSpace;
};
}  // namespace autopas

/**
 * @file FuzzyTuning.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "FuzzyTuning.h"

#include <sys/stat.h>

#include <numeric>
#include <utility>

#include "antlr4-runtime.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyControlSystem.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/FuzzyRuleErrorListener.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/TranslationVisitor.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageLexer.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageParser.h"

namespace autopas {

FuzzyTuning::FuzzyTuning(std::string fuzzyRuleFileName) : _fuzzyRuleFileName(std::move(fuzzyRuleFileName)) {
  // Check if the given rule file exists and throw if not
  struct stat buffer {};
  if (stat(_fuzzyRuleFileName.c_str(), &buffer) != 0) {
    utils::ExceptionHandler::exception("Rule file {} does not exist!", _fuzzyRuleFileName);
  }

  auto [linguisticVariables, fuzzyControlSystems] = parse(_fuzzyRuleFileName);

  // By default, dump the rules for reproducibility reasons.
  std::string linguisticVariablesStr =
      std::accumulate(linguisticVariables.begin(), linguisticVariables.end(), std::string("\n"),
                      [](const std::string &acc, const std::shared_ptr<fuzzy_logic::LinguisticVariable> &b) {
                        return acc + std::string(*b);
                      });

  std::string fuzzyControlSystemsStr =
      std::accumulate(fuzzyControlSystems.begin(), fuzzyControlSystems.end(), std::string("\n"),
                      [](const std::string &acc,
                         const std::pair<const std::string, std::shared_ptr<fuzzy_logic::FuzzyControlSystem>> &b) {
                        return acc + std::string(*b.second);
                      });

  AutoPasLog(INFO, "{}", linguisticVariablesStr);
  AutoPasLog(INFO, "{}", fuzzyControlSystemsStr);

  _fuzzyControlSystems = fuzzyControlSystems;
}

bool FuzzyTuning::needsLiveInfo() const { return true; }

void FuzzyTuning::receiveLiveInfo(const LiveInfo &info) {
  // clear the current live info
  _currentLiveInfo.clear();

  // only filter out the numeric values
  for (const auto &infoEntry : info.get()) {
    const auto &name = infoEntry.first;
    const auto &value = infoEntry.second;
    std::visit(
        [&](auto type) {
          if constexpr (std::is_same_v<decltype(type), size_t> or std::is_same_v<decltype(type), double>) {
            _currentLiveInfo[name] = type;
          }
        },
        value);
  }
}

void FuzzyTuning::addEvidence(const Configuration &configuration, const Evidence &evidence) {}

void FuzzyTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                        const EvidenceCollection &evidenceCollection) {}

void FuzzyTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                      const EvidenceCollection &evidenceCollection) {
  using namespace autopas::fuzzy_logic;

  for (const auto &[name, fcs] : _fuzzyControlSystems) {
    auto prediction = fcs->predict(_currentLiveInfo);

    std::cout << "Prediction for " << name << ": " << prediction << std::endl;
  }

  // set new search space
  // std::copy(newSearchSpace.begin(), newSearchSpace.end(), std::back_inserter(configQueue));
}

TuningStrategyOption FuzzyTuning::getOptionType() { return TuningStrategyOption::fuzzyTuning; }

std::pair<std::vector<std::shared_ptr<LinguisticVariable>>, std::map<std::string, std::shared_ptr<FuzzyControlSystem>>>
FuzzyTuning::parse(const std::string &fuzzyRuleFilename) {
  using namespace antlr4;
  using namespace autopas::fuzzy_logic;

  try {
    // Read file
    std::ifstream stream(fuzzyRuleFilename);
    ANTLRInputStream input(stream);
    FuzzyLanguageLexer lexer(&input);

    // Lexer
    FuzzyRuleErrorListener ErrorListener;
    lexer.removeErrorListeners();
    lexer.addErrorListener(&ErrorListener);

    CommonTokenStream tokens(&lexer);

    // Parser
    FuzzyLanguageParser parser(&tokens);
    parser.removeErrorListeners();
    parser.addErrorListener(&ErrorListener);

    tree::ParseTree *tree = parser.rule_file();

    // Translation
    TranslationVisitor visitor;
    auto fuzzy_rule_program =
        visitor.visit(tree)
            .as<std::pair<std::vector<std::shared_ptr<LinguisticVariable>>,
                          std::map<std::string, std::shared_ptr<fuzzy_logic::FuzzyControlSystem>>>>();

    return fuzzy_rule_program;
  } catch (const std::exception &e) {
    utils::ExceptionHandler::exception("Error while parsing fuzzy rule file: \"{}\". Error: {}", fuzzyRuleFilename,
                                       e.what());
    throw std::runtime_error("This should never be reached");
  }
}

}  // namespace autopas
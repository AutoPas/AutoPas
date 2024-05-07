/**
 * @file FuzzyTuning.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "FuzzyTuning.h"

#include <sys/stat.h>

#include "antlr4-runtime.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyControlSystem.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyRule.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzySetFactory.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/LinguisticVariable.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/xmlParser/XMLParser.h"
#include "parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageBaseVisitor.h"
#include "parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageLexer.h"
#include "parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageParser.h"

namespace autopas {

FuzzyTuning::FuzzyTuning(const std::string &fuzzyRuleFileName) : _fuzzyRuleFileName(fuzzyRuleFileName) {
  // Check if the given rule file exists and throw if not
  struct stat buffer;
  if (stat(_fuzzyRuleFileName.c_str(), &buffer) != 0) {
    utils::ExceptionHandler::exception("Rule file {} does not exist!", _fuzzyRuleFileName);
  }

  std::ifstream file(_fuzzyRuleFileName);

  std::string programCode((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

  antlr4::ANTLRInputStream input{programCode};
  autopas_generated_fuzzy_rule_syntax::FuzzyLanguageLexer lexer(&input);
  antlr4::CommonTokenStream tokens(&lexer);

  tokens.fill();

  autopas_generated_fuzzy_rule_syntax::FuzzyLanguageParser parser{&tokens};
  antlr4::tree::ParseTree *tree = parser.rule_file();

  auto visitor = std::make_unique<autopas_generated_fuzzy_rule_syntax::FuzzyLanguageBaseVisitor>();
  visitor->visit(tree);

  // parse(_ruleFileName);

  // By default, dump the rules for reproducibility reasons.
  // AutoPasLog(INFO, "Rule File {}:\n{}", _fuzzyRuleFileName, rulesToString(_fuzzyRuleFileName));
}

bool FuzzyTuning::needsLiveInfo() const { return true; }

void FuzzyTuning::receiveLiveInfo(const LiveInfo &info) { _currentLiveInfo = info; }

void FuzzyTuning::addEvidence(const Configuration &configuration, const Evidence &evidence) {}

void FuzzyTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                        const EvidenceCollection &evidenceCollection) {}

void FuzzyTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                      const EvidenceCollection &evidenceCollection) {
  using namespace autopas::fuzzy_logic;

  const std::map<std::string, LiveInfo::InfoType> live_info = _currentLiveInfo.get();

  auto spreadParticlesPerCell = LinguisticVariable("spreadParticlesPerCell", std::pair(0, 10));
  auto avgParticles = LinguisticVariable("avgParticlesPerCell", std::pair(0, 10));

  auto density = LinguisticVariable("density", std::pair(0, 30));

  spreadParticlesPerCell.addLinguisticTerm(makeGaussian("low", 0, 2));
  spreadParticlesPerCell.addLinguisticTerm(makeGaussian("medium", 5, 2));
  spreadParticlesPerCell.addLinguisticTerm(makeGaussian("high", 10, 2));

  avgParticles.addLinguisticTerm(makeGaussian("low", 0, 2));
  avgParticles.addLinguisticTerm(makeGaussian("medium", 5, 2));
  avgParticles.addLinguisticTerm(makeGaussian("high", 10, 2));

  density.addLinguisticTerm(makeGaussian("low", 0, 10));
  density.addLinguisticTerm(makeGaussian("medium", 15, 10));
  density.addLinguisticTerm(makeGaussian("high", 30, 10));

  auto fcs = FuzzyControlSystem();

  fcs.addRule(FuzzyRule(avgParticles == "low" && spreadParticlesPerCell == "low", density == "low"));
  fcs.addRule(FuzzyRule(avgParticles == "low" && spreadParticlesPerCell == "medium", density == "medium"));

  size_t spread =
      std::get<size_t>(live_info.at("maxParticlesPerCell")) - std::get<size_t>(live_info.at("minParticlesPerCell"));

  FuzzySet::Data sample = {{"avgParticlesPerCell", std::get<double>(live_info.at("avgParticlesPerCell"))},
                           {"spreadParticlesPerCell", spread}};

  auto x2 = fcs.predict(sample);

  std::cout << "Predicted tip: " << x2 << std::endl;

  // set new search space
  // std::copy(newSearchSpace.begin(), newSearchSpace.end(), std::back_inserter(configQueue));
}

TuningStrategyOption FuzzyTuning::getOptionType() { return TuningStrategyOption::fuzzyTuning; }

}  // namespace autopas
/**
 * @file TranslationVisitor.cpp
 * @author Manuel Lerchner
 * @date 09.05.24
 */

#include "TranslationVisitor.h"

#include "autopas/tuning/tuningStrategy/fuzzyTuning/OutputMapper.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzySetFactory.h"
#include "autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedProgramTree.h"
#include "autopas/utils/ExceptionHandler.h"

using namespace autopas;
using std::any_cast;

antlrcpp::Any TranslationVisitor::visitRule_file(FuzzyLanguageParser::Rule_fileContext *context) {
  std::vector<std::shared_ptr<LinguisticVariable>> linguisticVariables;
  std::vector<FuzzyRule> fuzzyRules;

  // get the settings
  const auto settings = any_cast<std::shared_ptr<FuzzyControlSettings>>(visit(context->settings()));

  // get the linguistic variables
  for (auto *fuzzy_variableContext : context->linguistic_variable()) {
    const auto lv = any_cast<std::shared_ptr<LinguisticVariable>>(visit(fuzzy_variableContext));
    linguisticVariables.push_back(lv);

    _state._linguisticVariables[lv->getName()] = lv;
  }

  // get the configuration mappings
  const auto outputMapping =
      any_cast<std::map<std::string, std::shared_ptr<OutputMapper>>>(visit(context->output_mapping()));

  // get the fuzzy rules
  for (auto *fuzzy_ruleContext : context->fuzzy_rule()) {
    const auto rule = any_cast<FuzzyRule>(visit(fuzzy_ruleContext));
    fuzzyRules.push_back(rule);
  }

  // construct the fuzzy control systems
  std::map<std::string, std::shared_ptr<FuzzyControlSystem>> fuzzyControlSystems;

  for (auto &rule : fuzzyRules) {
    const auto dimensions = rule.getConsequent()->getCrispSet()->getDimensions();
    if (dimensions.size() != 1) {
      utils::ExceptionHandler::exception("Only rules with one dimensional output are supported! Rule: " +
                                         std::string(rule));
    }
    const auto output_domain = dimensions.begin()->first;

    if (fuzzyControlSystems.find(output_domain) == fuzzyControlSystems.end()) {
      fuzzyControlSystems[output_domain] = std::make_shared<FuzzyControlSystem>(settings);
    }

    fuzzyControlSystems[output_domain]->addRule(rule);
  }

  return std::make_tuple(settings, linguisticVariables, outputMapping, fuzzyControlSystems);
};

antlrcpp::Any TranslationVisitor::visitSettings(FuzzyLanguageParser::SettingsContext *ctx) {
  auto settings = std::make_shared<FuzzyControlSettings>();

  for (size_t i = 0; i < ctx->STRING().size(); ++i) {
    const auto key = ctx->IDENTIFIER(i)->getText();
    const auto value = ctx->STRING(i)->getText();

    settings->insert({key, value});
  }

  return settings;
}

antlrcpp::Any TranslationVisitor::visitLinguistic_variable(FuzzyLanguageParser::Linguistic_variableContext *context) {
  // get the name and range of the linguistic variable
  const auto linguisticTerm = context->STRING()->getText();
  const std::pair<double, double> range = {std::stod(context->NUMBER(0)->getText()),
                                           std::stod(context->NUMBER(1)->getText())};

  auto linguisticVariable = std::make_shared<LinguisticVariable>(linguisticTerm, range);

  // add all linguistic terms
  for (auto *fuzzy_termContext : context->fuzzy_term()) {
    const auto term = any_cast<std::shared_ptr<FuzzySet>>(visit(fuzzy_termContext));
    linguisticVariable->addLinguisticTerm(term);
  }

  return linguisticVariable;
};

antlrcpp::Any TranslationVisitor::visitFuzzy_term(FuzzyLanguageParser::Fuzzy_termContext *context) {
  const std::string linguisticTerm = context->STRING()->getText();

  const auto function = visit(context->function());
  const auto [functionName, params] = any_cast<std::pair<std::string, std::vector<double>>>(function);

  return FuzzySetFactory::makeFuzzySet(linguisticTerm, functionName, params);
};

antlrcpp::Any TranslationVisitor::visitFunction(FuzzyLanguageParser::FunctionContext *context) {
  const std::string function = context->IDENTIFIER()->getText();
  std::vector<double> params;
  for (auto *number : context->NUMBER()) {
    params.push_back(std::stod(number->getText()));
  }

  return std::pair(function, params);
};

antlrcpp::Any TranslationVisitor::visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *context) {
  const auto antecedent = any_cast<std::shared_ptr<FuzzySet>>(visit(context->fuzzy_set(0)));
  const auto consequent = any_cast<std::shared_ptr<FuzzySet>>(visit(context->fuzzy_set(1)));

  return FuzzyRule(antecedent, consequent);
};

antlrcpp::Any TranslationVisitor::visitOr(FuzzyLanguageParser::OrContext *context) {
  const auto left = any_cast<std::shared_ptr<FuzzySet>>(visit(context->fuzzy_set(0)));
  const auto right = any_cast<std::shared_ptr<FuzzySet>>(visit(context->fuzzy_set(1)));

  return left || right;
};

antlrcpp::Any TranslationVisitor::visitBrackets(FuzzyLanguageParser::BracketsContext *context) {
  return visit(context->fuzzy_set());
};

antlrcpp::Any TranslationVisitor::visitAnd(FuzzyLanguageParser::AndContext *context) {
  const auto left = any_cast<std::shared_ptr<FuzzySet>>(visit(context->fuzzy_set(0)));
  const auto right = any_cast<std::shared_ptr<FuzzySet>>(visit(context->fuzzy_set(1)));

  return left && right;
};

antlrcpp::Any TranslationVisitor::visitSelect(FuzzyLanguageParser::SelectContext *context) {
  const auto lvName = context->STRING(0)->getText();
  const auto termName = context->STRING(1)->getText();

  const auto lv = _state._linguisticVariables[lvName];

  return lv->operator==(termName);
};

antlrcpp::Any TranslationVisitor::visitNegate(FuzzyLanguageParser::NegateContext *context) {
  const auto fuzzySet = any_cast<std::shared_ptr<FuzzySet>>(visit(context->fuzzy_set()));
  return !fuzzySet;
};

antlrcpp::Any TranslationVisitor::visitOutput_mapping(FuzzyLanguageParser::Output_mappingContext *context) {
  std::map<std::string, std::shared_ptr<OutputMapper>> outputMappings;

  for (size_t i = 0; i < context->output_entry().size(); ++i) {
    const auto mapping = any_cast<std::shared_ptr<OutputMapper>>(visit(context->output_entry(i)));

    outputMappings[mapping->getOutputDomain()] = mapping;
  }

  return outputMappings;
}

antlrcpp::Any TranslationVisitor::visitOutput_entry(FuzzyLanguageParser::Output_entryContext *context) {
  const auto pattern = context->STRING()->getText();

  std::vector<std::pair<double, std::vector<ConfigurationPattern>>> mappings;

  for (size_t i = 0; i < context->pattern_mapping().size(); ++i) {
    const auto mapping =
        any_cast<std::pair<double, std::vector<ConfigurationPattern>>>(visit(context->pattern_mapping(i)));
    mappings.push_back(mapping);
  }

  return std::make_shared<OutputMapper>(pattern, mappings);
}
antlrcpp::Any TranslationVisitor::visitPattern_mapping(FuzzyLanguageParser::Pattern_mappingContext *context) {
  const double value = std::stod(context->NUMBER()->getText());
  std::vector<ConfigurationPattern> configurationPatterns;

  for (size_t i = 0; i < context->configuration_pattern().size(); ++i) {
    auto *val = context->configuration_pattern(i);
    const auto configurationPattern = any_cast<ConfigurationPattern>(visit(val));
    configurationPatterns.push_back(configurationPattern);
  }

  return std::make_pair(value, configurationPatterns);
}

auto makeLiteral(const std::string &optionType, const std::string &optionValue) {
  RuleVM::MemoryCell literal;
  if (optionType == "container") {
    literal = ContainerOption::parseOptionExact(optionValue);
  } else if (optionType == "traversal") {
    literal = TraversalOption::parseOptionExact(optionValue);
  } else if (optionType == "dataLayout") {
    literal = DataLayoutOption::parseOptionExact(optionValue);
  } else if (optionType == "loadEstimator") {
    literal = LoadEstimatorOption::parseOptionExact(optionValue);
  } else if (optionType == "newton3") {
    literal = Newton3Option::parseOptionExact(optionValue);
  } else if (optionType == "cellSizeFactor") {
    literal = RuleVM::MemoryCell{std::stod(optionValue)};
  } else {
    throw std::runtime_error("Unknown option type: " + optionType + " with value: " + optionValue);
  }

  return std::make_shared<Literal>(literal);
}

antlrcpp::Any TranslationVisitor::visitConfiguration_pattern(
    FuzzyLanguageParser::Configuration_patternContext *context) {
  autopas::ConfigurationPattern pattern;

  for (size_t i = 0; i < context->IDENTIFIER().size(); ++i) {
    auto property = context->IDENTIFIER(i)->getText();
    auto value = context->STRING(i)->getText();

    auto literal = makeLiteral(property, value);

    pattern.add(literal->value);
    // sanity check
    const auto knownProperties = {
        "container", "traversal", "dataLayout", "newton3", "loadEstimator", "cellSizeFactor",
    };
    if (std::find(knownProperties.begin(), knownProperties.end(), property) == knownProperties.end()) {
      autopas::utils::ExceptionHandler::exception("RuleBasedProgramParser: Encountered unknown property! (" + property +
                                                  ")");
    }
  }

  return pattern;
}

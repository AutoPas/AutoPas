
// Generated from
// autopas/tuning/tuningStrategy/fuzzyTuning\FuzzyLanguage.g4 by
// ANTLR 4.9.1

#pragma once

#include "FuzzyLanguageParser.h"
#include "antlr4-runtime.h"

namespace autopas_generated_fuzzy_rule_syntax {

/**
 * This class defines an abstract visitor for a parse tree
 * produced by FuzzyLanguageParser.
 */
class FuzzyLanguageVisitor : public antlr4::tree::AbstractParseTreeVisitor {
 public:
  /**
   * Visit parse trees produced by FuzzyLanguageParser.
   */
  virtual antlrcpp::Any visitRule_file(FuzzyLanguageParser::Rule_fileContext *context) = 0;

  virtual antlrcpp::Any visitFuzzy_variable(FuzzyLanguageParser::Fuzzy_variableContext *context) = 0;

  virtual antlrcpp::Any visitMembership_function(FuzzyLanguageParser::Membership_functionContext *context) = 0;

  virtual antlrcpp::Any visitFunction(FuzzyLanguageParser::FunctionContext *context) = 0;

  virtual antlrcpp::Any visitName(FuzzyLanguageParser::NameContext *context) = 0;

  virtual antlrcpp::Any visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *context) = 0;

  virtual antlrcpp::Any visitFuzzy_set(FuzzyLanguageParser::Fuzzy_setContext *context) = 0;

  virtual antlrcpp::Any visitSelection(FuzzyLanguageParser::SelectionContext *context) = 0;
};

}  // namespace autopas_generated_fuzzy_rule_syntax


// Generated from
// autopas/tuning/tuningStrategy/fuzzyTuning\FuzzyLanguage.g4 by
// ANTLR 4.9.1

#pragma once

#include "FuzzyLanguageParser.h"
#include "antlr4-runtime.h"

namespace AutopasGeneratedFuzzyRuleSyntax {

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

  virtual antlrcpp::Any visitSettings(FuzzyLanguageParser::SettingsContext *context) = 0;

  virtual antlrcpp::Any visitLinguistic_variable(FuzzyLanguageParser::Linguistic_variableContext *context) = 0;

  virtual antlrcpp::Any visitFuzzy_term(FuzzyLanguageParser::Fuzzy_termContext *context) = 0;

  virtual antlrcpp::Any visitFunction(FuzzyLanguageParser::FunctionContext *context) = 0;

  virtual antlrcpp::Any visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *context) = 0;

  virtual antlrcpp::Any visitOr(FuzzyLanguageParser::OrContext *context) = 0;

  virtual antlrcpp::Any visitBrackets(FuzzyLanguageParser::BracketsContext *context) = 0;

  virtual antlrcpp::Any visitAnd(FuzzyLanguageParser::AndContext *context) = 0;

  virtual antlrcpp::Any visitSelect(FuzzyLanguageParser::SelectContext *context) = 0;

  virtual antlrcpp::Any visitNegate(FuzzyLanguageParser::NegateContext *context) = 0;

  virtual antlrcpp::Any visitOutput_mapping(FuzzyLanguageParser::Output_mappingContext *context) = 0;

  virtual antlrcpp::Any visitOutput_entry(FuzzyLanguageParser::Output_entryContext *context) = 0;

  virtual antlrcpp::Any visitPattern_mapping(FuzzyLanguageParser::Pattern_mappingContext *context) = 0;

  virtual antlrcpp::Any visitConfiguration_pattern(FuzzyLanguageParser::Configuration_patternContext *context) = 0;
};

}  // namespace AutopasGeneratedFuzzyRuleSyntax

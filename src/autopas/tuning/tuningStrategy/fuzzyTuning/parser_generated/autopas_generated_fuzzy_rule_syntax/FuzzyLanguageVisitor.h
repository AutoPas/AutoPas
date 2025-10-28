
// Generated from FuzzyLanguage.g4 by ANTLR 4.13.2

#pragma once

#include <any>

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
  virtual std::any visitRule_file(FuzzyLanguageParser::Rule_fileContext *context) = 0;

  virtual std::any visitSettings(FuzzyLanguageParser::SettingsContext *context) = 0;

  virtual std::any visitLinguistic_variable(FuzzyLanguageParser::Linguistic_variableContext *context) = 0;

  virtual std::any visitFuzzy_term(FuzzyLanguageParser::Fuzzy_termContext *context) = 0;

  virtual std::any visitFunction(FuzzyLanguageParser::FunctionContext *context) = 0;

  virtual std::any visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *context) = 0;

  virtual std::any visitOr(FuzzyLanguageParser::OrContext *context) = 0;

  virtual std::any visitBrackets(FuzzyLanguageParser::BracketsContext *context) = 0;

  virtual std::any visitAnd(FuzzyLanguageParser::AndContext *context) = 0;

  virtual std::any visitSelect(FuzzyLanguageParser::SelectContext *context) = 0;

  virtual std::any visitNegate(FuzzyLanguageParser::NegateContext *context) = 0;

  virtual std::any visitOutput_mapping(FuzzyLanguageParser::Output_mappingContext *context) = 0;

  virtual std::any visitOutput_entry(FuzzyLanguageParser::Output_entryContext *context) = 0;

  virtual std::any visitPattern_mapping(FuzzyLanguageParser::Pattern_mappingContext *context) = 0;

  virtual std::any visitConfiguration_pattern(FuzzyLanguageParser::Configuration_patternContext *context) = 0;
};

}  // namespace AutopasGeneratedFuzzyRuleSyntax

/**
 * @file TranslationVisitor.h
 * @author Manuel Lerchner
 * @date 09.05.24
 */

#pragma once

#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyControlSystem.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyRule.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzySet.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzySetFactory.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/LinguisticVariable.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageBaseVisitor.h"
#include "autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedProgramParser.h"

using namespace AutopasGeneratedFuzzyRuleSyntax;
using namespace autopas::FuzzyLogic;
using namespace autopas::RuleSyntax;

/**
 * This class implements a visitor for the fuzzy rule parser.
 * The visitRule_file method is the entry point for the visitor and translates the parse tree into linguistic
 * variables and FuzzyControlSystems.
 */
class TranslationVisitor : public FuzzyLanguageVisitor {
 public:
  /**
   * Visit parse trees produced by FuzzyLanguageParser.
   */
  antlrcpp::Any visitRule_file(FuzzyLanguageParser::Rule_fileContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::System_definitionsContext and create a FuzzyControlSystem.
   */
  antlrcpp::Any visitSettings(FuzzyLanguageParser::SettingsContext *ctx) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::Linguistic_variableContext and create a LinguisticVariable.
   */
  antlrcpp::Any visitLinguistic_variable(FuzzyLanguageParser::Linguistic_variableContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::Fuzzy_termContext and create a FuzzySet.
   */
  antlrcpp::Any visitFuzzy_term(FuzzyLanguageParser::Fuzzy_termContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::FunctionContext and collect the function parameters.
   */
  antlrcpp::Any visitFunction(FuzzyLanguageParser::FunctionContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::Fuzzy_ruleContext and create a FuzzyRule.
   */
  antlrcpp::Any visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::Fuzzy_setContext and create a FuzzySet based on the
   * OR operator between the two fuzzy sets.
   */
  antlrcpp::Any visitOr(FuzzyLanguageParser::OrContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::Fuzzy_setContext and create a FuzzySet from the
   * body of the bracket.
   */
  antlrcpp::Any visitBrackets(FuzzyLanguageParser::BracketsContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::BracketsContext and create a FuzzySet based on the
   * AND operator between the two fuzzy sets.
   */
  antlrcpp::Any visitAnd(FuzzyLanguageParser::AndContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::BracketsContext and create a FuzzySet based on the
   * SELECT operator of the LinguisticVariable and a term.
   */
  antlrcpp::Any visitSelect(FuzzyLanguageParser::SelectContext *context) override;

  /**
   * Visit a parse tree produced by FuzzyLanguageParser::BracketsContext and create a FuzzySet based on the
   * NOT operator of the fuzzy set.
   */
  antlrcpp::Any visitNegate(FuzzyLanguageParser::NegateContext *context) override;

  /**
   * Visit a configuration mapping produced by FuzzyLanguageParser::Configuration_mappingContext and create a
   * configuration mapping.
   */
  antlrcpp::Any visitOutput_mapping(FuzzyLanguageParser::Output_mappingContext *context) override;

  /**
   * Visit a configuration entry produced by FuzzyLanguageParser::Configuration_entryContext and create a
   * configuration entry.
   */
  antlrcpp::Any visitOutput_entry(FuzzyLanguageParser::Output_entryContext *context) override;
  /**
   * Visit a pattern mapping produced by FuzzyLanguageParser::Pattern_mappingContext and create a
   * pattern mapping.
   */
  antlrcpp::Any visitPattern_mapping(FuzzyLanguageParser::Pattern_mappingContext *context) override;

  /**
   * Visit an assignment produced by FuzzyLanguageParser::AssignmentContext and create an assignment.
   */
  antlrcpp::Any visitConfiguration_pattern(FuzzyLanguageParser::Configuration_patternContext *context) override;

 private:
  /**
   * The internal state of the visitor.
   */
  struct VisitorState {
    std::map<std::string, std::shared_ptr<LinguisticVariable>> _linguisticVariables;
  } _state;
};

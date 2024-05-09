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
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageVisitor.h"

using namespace autopas_generated_fuzzy_rule_syntax;
using namespace autopas::fuzzy_logic;

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

 private:
  /**
   * The internal state of the visitor.
   */
  struct VisitorState {
    std::map<std::string, std::shared_ptr<LinguisticVariable>> _linguisticVariables;
  } _state;
};

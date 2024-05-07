
// Generated from
// AutoPas/src/autopas/tuning/tuningStrategy/fuzzyTuning\FuzzyLanguage.g4 by
// ANTLR 4.9.1

#pragma once

#include "FuzzyLanguageVisitor.h"
#include "antlr4-runtime.h"

namespace autopas_generated_fuzzy_rule_syntax {

/**
 * This class provides an empty implementation of FuzzyLanguageVisitor, which can be
 * extended to create a visitor which only needs to handle a subset of the available methods.
 */
class FuzzyLanguageBaseVisitor : public FuzzyLanguageVisitor {
 public:
  virtual antlrcpp::Any visitRule_file(FuzzyLanguageParser::Rule_fileContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitFuzzy_variable(FuzzyLanguageParser::Fuzzy_variableContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitMembership_function(FuzzyLanguageParser::Membership_functionContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitFunction(FuzzyLanguageParser::FunctionContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitName(FuzzyLanguageParser::NameContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitFuzzy_set(FuzzyLanguageParser::Fuzzy_setContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitSelection(FuzzyLanguageParser::SelectionContext *ctx) override {
    return visitChildren(ctx);
  }
};

}  // namespace autopas_generated_fuzzy_rule_syntax

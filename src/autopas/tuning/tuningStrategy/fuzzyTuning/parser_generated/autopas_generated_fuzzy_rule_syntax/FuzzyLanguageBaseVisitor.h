
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

  virtual antlrcpp::Any visitSettings(FuzzyLanguageParser::SettingsContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitLinguistic_variable(FuzzyLanguageParser::Linguistic_variableContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitFuzzy_term(FuzzyLanguageParser::Fuzzy_termContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitFunction(FuzzyLanguageParser::FunctionContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitOr(FuzzyLanguageParser::OrContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitBrackets(FuzzyLanguageParser::BracketsContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitAnd(FuzzyLanguageParser::AndContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitSelect(FuzzyLanguageParser::SelectContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitNegate(FuzzyLanguageParser::NegateContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitOutput_mapping(FuzzyLanguageParser::Output_mappingContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitOutput_entry(FuzzyLanguageParser::Output_entryContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitPattern_mapping(FuzzyLanguageParser::Pattern_mappingContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitConfiguration_pattern(FuzzyLanguageParser::Configuration_patternContext *ctx) override {
    return visitChildren(ctx);
  }
};

}  // namespace autopas_generated_fuzzy_rule_syntax

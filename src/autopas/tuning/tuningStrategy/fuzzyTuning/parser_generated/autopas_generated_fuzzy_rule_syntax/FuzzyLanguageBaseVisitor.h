
// Generated from FuzzyLanguage.g4 by ANTLR 4.13.2

#pragma once

#include <any>

#include "FuzzyLanguageVisitor.h"
#include "antlr4-runtime.h"

namespace AutopasGeneratedFuzzyRuleSyntax {

/**
 * This class provides an empty implementation of FuzzyLanguageVisitor, which can be
 * extended to create a visitor which only needs to handle a subset of the available methods.
 */
class FuzzyLanguageBaseVisitor : public FuzzyLanguageVisitor {
 public:
  virtual std::any visitRule_file(FuzzyLanguageParser::Rule_fileContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitSettings(FuzzyLanguageParser::SettingsContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitLinguistic_variable(FuzzyLanguageParser::Linguistic_variableContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitFuzzy_term(FuzzyLanguageParser::Fuzzy_termContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitFunction(FuzzyLanguageParser::FunctionContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitOr(FuzzyLanguageParser::OrContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitBrackets(FuzzyLanguageParser::BracketsContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitAnd(FuzzyLanguageParser::AndContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitSelect(FuzzyLanguageParser::SelectContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitNegate(FuzzyLanguageParser::NegateContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitOutput_mapping(FuzzyLanguageParser::Output_mappingContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitOutput_entry(FuzzyLanguageParser::Output_entryContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitPattern_mapping(FuzzyLanguageParser::Pattern_mappingContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitConfiguration_pattern(FuzzyLanguageParser::Configuration_patternContext *ctx) override {
    return visitChildren(ctx);
  }
};

}  // namespace AutopasGeneratedFuzzyRuleSyntax

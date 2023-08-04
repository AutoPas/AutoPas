
// Generated from /home/tobias/AutoPas2/src/autopas/selectors/tuningStrategy/ruleBasedTuning/RuleLanguage.g4 by
// ANTLR 4.9.1

#pragma once

#include "RuleLanguageVisitor.h"
#include "antlr4-runtime.h"

namespace autopas_generated_rule_syntax {

/**
 * This class provides an empty implementation of RuleLanguageVisitor, which can be
 * extended to create a visitor which only needs to handle a subset of the available methods.
 */
class RuleLanguageBaseVisitor : public RuleLanguageVisitor {
 public:
  virtual antlrcpp::Any visitProgram(RuleLanguageParser::ProgramContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitUnsigned_val(RuleLanguageParser::Unsigned_valContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitLiteral(RuleLanguageParser::LiteralContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitDefine_list(RuleLanguageParser::Define_listContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitVariable(RuleLanguageParser::VariableContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitExpression(RuleLanguageParser::ExpressionContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitDefine(RuleLanguageParser::DefineContext *ctx) override { return visitChildren(ctx); }

  virtual antlrcpp::Any visitProperty_value(RuleLanguageParser::Property_valueContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitConfiguration_order(RuleLanguageParser::Configuration_orderContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitStatement(RuleLanguageParser::StatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual antlrcpp::Any visitIf_statement(RuleLanguageParser::If_statementContext *ctx) override {
    return visitChildren(ctx);
  }
};

}  // namespace autopas_generated_rule_syntax

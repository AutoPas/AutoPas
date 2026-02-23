
// Generated from RuleLanguage.g4 by ANTLR 4.13.2

#pragma once

#include <any>

#include "RuleLanguageVisitor.h"
#include "antlr4-runtime.h"

namespace AutopasGeneratedRuleSyntax {

/**
 * This class provides an empty implementation of RuleLanguageVisitor, which can be
 * extended to create a visitor which only needs to handle a subset of the available methods.
 */
class RuleLanguageBaseVisitor : public RuleLanguageVisitor {
 public:
  virtual std::any visitProgram(RuleLanguageParser::ProgramContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitUnsigned_val(RuleLanguageParser::Unsigned_valContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitLiteral(RuleLanguageParser::LiteralContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitDefine_list(RuleLanguageParser::Define_listContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitVariable(RuleLanguageParser::VariableContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitExpression(RuleLanguageParser::ExpressionContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitDefine(RuleLanguageParser::DefineContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitProperty_value(RuleLanguageParser::Property_valueContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitConfiguration_order(RuleLanguageParser::Configuration_orderContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitStatement(RuleLanguageParser::StatementContext *ctx) override { return visitChildren(ctx); }

  virtual std::any visitIf_statement(RuleLanguageParser::If_statementContext *ctx) override {
    return visitChildren(ctx);
  }
};

}  // namespace AutopasGeneratedRuleSyntax

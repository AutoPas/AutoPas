
// Generated from RuleLanguage.g4 by ANTLR 4.13.2

#pragma once

#include <any>

#include "RuleLanguageParser.h"
#include "antlr4-runtime.h"

namespace AutopasGeneratedRuleSyntax {

/**
 * This class defines an abstract visitor for a parse tree
 * produced by RuleLanguageParser.
 */
class RuleLanguageVisitor : public antlr4::tree::AbstractParseTreeVisitor {
 public:
  /**
   * Visit parse trees produced by RuleLanguageParser.
   */
  virtual std::any visitProgram(RuleLanguageParser::ProgramContext *context) = 0;

  virtual std::any visitUnsigned_val(RuleLanguageParser::Unsigned_valContext *context) = 0;

  virtual std::any visitLiteral(RuleLanguageParser::LiteralContext *context) = 0;

  virtual std::any visitDefine_list(RuleLanguageParser::Define_listContext *context) = 0;

  virtual std::any visitVariable(RuleLanguageParser::VariableContext *context) = 0;

  virtual std::any visitExpression(RuleLanguageParser::ExpressionContext *context) = 0;

  virtual std::any visitDefine(RuleLanguageParser::DefineContext *context) = 0;

  virtual std::any visitProperty_value(RuleLanguageParser::Property_valueContext *context) = 0;

  virtual std::any visitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *context) = 0;

  virtual std::any visitConfiguration_order(RuleLanguageParser::Configuration_orderContext *context) = 0;

  virtual std::any visitStatement(RuleLanguageParser::StatementContext *context) = 0;

  virtual std::any visitIf_statement(RuleLanguageParser::If_statementContext *context) = 0;
};

}  // namespace AutopasGeneratedRuleSyntax

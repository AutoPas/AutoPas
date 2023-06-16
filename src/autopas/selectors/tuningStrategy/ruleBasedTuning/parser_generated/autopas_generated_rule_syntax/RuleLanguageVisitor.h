
// Generated from /home/tobias/AutoPas2/src/autopas/selectors/tuningStrategy/ruleBasedTuning/RuleLanguage.g4 by
// ANTLR 4.9.1

#pragma once

#include "RuleLanguageParser.h"
#include "antlr4-runtime.h"

namespace autopas_generated_rule_syntax {

/**
 * This class defines an abstract visitor for a parse tree
 * produced by RuleLanguageParser.
 */
class RuleLanguageVisitor : public antlr4::tree::AbstractParseTreeVisitor {
 public:
  /**
   * Visit parse trees produced by RuleLanguageParser.
   */
  virtual antlrcpp::Any visitProgram(RuleLanguageParser::ProgramContext *context) = 0;

  virtual antlrcpp::Any visitUnsigned_val(RuleLanguageParser::Unsigned_valContext *context) = 0;

  virtual antlrcpp::Any visitLiteral(RuleLanguageParser::LiteralContext *context) = 0;

  virtual antlrcpp::Any visitDefine_list(RuleLanguageParser::Define_listContext *context) = 0;

  virtual antlrcpp::Any visitVariable(RuleLanguageParser::VariableContext *context) = 0;

  virtual antlrcpp::Any visitExpression(RuleLanguageParser::ExpressionContext *context) = 0;

  virtual antlrcpp::Any visitDefine(RuleLanguageParser::DefineContext *context) = 0;

  virtual antlrcpp::Any visitProperty_value(RuleLanguageParser::Property_valueContext *context) = 0;

  virtual antlrcpp::Any visitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *context) = 0;

  virtual antlrcpp::Any visitConfiguration_order(RuleLanguageParser::Configuration_orderContext *context) = 0;

  virtual antlrcpp::Any visitStatement(RuleLanguageParser::StatementContext *context) = 0;

  virtual antlrcpp::Any visitIf_statement(RuleLanguageParser::If_statementContext *context) = 0;
};

}  // namespace autopas_generated_rule_syntax


// Generated from /home/tobias/AutoPas2/src/autopas/selectors/tuningStrategy/ruleBasedTuning/RuleLanguage.g4 by ANTLR 4.9.1

#pragma once


#include "antlr4-runtime.h"
#include "RuleLanguageParser.h"


namespace autopas::rule_syntax {

/**
 * This interface defines an abstract listener for a parse tree produced by RuleLanguageParser.
 */
class  RuleLanguageListener : public antlr4::tree::ParseTreeListener {
public:

  virtual void enterProgram(RuleLanguageParser::ProgramContext *ctx) = 0;
  virtual void exitProgram(RuleLanguageParser::ProgramContext *ctx) = 0;

  virtual void enterLiteral(RuleLanguageParser::LiteralContext *ctx) = 0;
  virtual void exitLiteral(RuleLanguageParser::LiteralContext *ctx) = 0;

  virtual void enterDefine_list(RuleLanguageParser::Define_listContext *ctx) = 0;
  virtual void exitDefine_list(RuleLanguageParser::Define_listContext *ctx) = 0;

  virtual void enterDefine(RuleLanguageParser::DefineContext *ctx) = 0;
  virtual void exitDefine(RuleLanguageParser::DefineContext *ctx) = 0;

  virtual void enterVariable(RuleLanguageParser::VariableContext *ctx) = 0;
  virtual void exitVariable(RuleLanguageParser::VariableContext *ctx) = 0;

  virtual void enterExpression(RuleLanguageParser::ExpressionContext *ctx) = 0;
  virtual void exitExpression(RuleLanguageParser::ExpressionContext *ctx) = 0;

  virtual void enterProperty_value(RuleLanguageParser::Property_valueContext *ctx) = 0;
  virtual void exitProperty_value(RuleLanguageParser::Property_valueContext *ctx) = 0;

  virtual void enterConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *ctx) = 0;
  virtual void exitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *ctx) = 0;

  virtual void enterConfiguration_order(RuleLanguageParser::Configuration_orderContext *ctx) = 0;
  virtual void exitConfiguration_order(RuleLanguageParser::Configuration_orderContext *ctx) = 0;

  virtual void enterStatement(RuleLanguageParser::StatementContext *ctx) = 0;
  virtual void exitStatement(RuleLanguageParser::StatementContext *ctx) = 0;

  virtual void enterIf_statement(RuleLanguageParser::If_statementContext *ctx) = 0;
  virtual void exitIf_statement(RuleLanguageParser::If_statementContext *ctx) = 0;


};

}  // namespace autopas_rule_syntax

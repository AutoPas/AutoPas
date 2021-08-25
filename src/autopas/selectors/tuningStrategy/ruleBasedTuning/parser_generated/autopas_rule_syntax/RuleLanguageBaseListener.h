
// Generated from /home/tobias/AutoPas2/src/autopas/selectors/tuningStrategy/ruleBasedTuning/RuleLanguage.g4 by ANTLR 4.9.1

#pragma once


#include "antlr4-runtime.h"
#include "RuleLanguageListener.h"


namespace autopas::rule_syntax {

/**
 * This class provides an empty implementation of RuleLanguageListener,
 * which can be extended to create a listener which only needs to handle a subset
 * of the available methods.
 */
class  RuleLanguageBaseListener : public RuleLanguageListener {
public:

  virtual void enterProgram(RuleLanguageParser::ProgramContext * /*ctx*/) override { }
  virtual void exitProgram(RuleLanguageParser::ProgramContext * /*ctx*/) override { }

  virtual void enterLiteral(RuleLanguageParser::LiteralContext * /*ctx*/) override { }
  virtual void exitLiteral(RuleLanguageParser::LiteralContext * /*ctx*/) override { }

  virtual void enterDefine_list(RuleLanguageParser::Define_listContext * /*ctx*/) override { }
  virtual void exitDefine_list(RuleLanguageParser::Define_listContext * /*ctx*/) override { }

  virtual void enterDefine(RuleLanguageParser::DefineContext * /*ctx*/) override { }
  virtual void exitDefine(RuleLanguageParser::DefineContext * /*ctx*/) override { }

  virtual void enterVariable(RuleLanguageParser::VariableContext * /*ctx*/) override { }
  virtual void exitVariable(RuleLanguageParser::VariableContext * /*ctx*/) override { }

  virtual void enterExpression(RuleLanguageParser::ExpressionContext * /*ctx*/) override { }
  virtual void exitExpression(RuleLanguageParser::ExpressionContext * /*ctx*/) override { }

  virtual void enterProperty_value(RuleLanguageParser::Property_valueContext * /*ctx*/) override { }
  virtual void exitProperty_value(RuleLanguageParser::Property_valueContext * /*ctx*/) override { }

  virtual void enterConfiguration_pattern(RuleLanguageParser::Configuration_patternContext * /*ctx*/) override { }
  virtual void exitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext * /*ctx*/) override { }

  virtual void enterConfiguration_order(RuleLanguageParser::Configuration_orderContext * /*ctx*/) override { }
  virtual void exitConfiguration_order(RuleLanguageParser::Configuration_orderContext * /*ctx*/) override { }

  virtual void enterStatement(RuleLanguageParser::StatementContext * /*ctx*/) override { }
  virtual void exitStatement(RuleLanguageParser::StatementContext * /*ctx*/) override { }

  virtual void enterIf_statement(RuleLanguageParser::If_statementContext * /*ctx*/) override { }
  virtual void exitIf_statement(RuleLanguageParser::If_statementContext * /*ctx*/) override { }


  virtual void enterEveryRule(antlr4::ParserRuleContext * /*ctx*/) override { }
  virtual void exitEveryRule(antlr4::ParserRuleContext * /*ctx*/) override { }
  virtual void visitTerminal(antlr4::tree::TerminalNode * /*node*/) override { }
  virtual void visitErrorNode(antlr4::tree::ErrorNode * /*node*/) override { }

};

}  // namespace autopas_rule_syntax

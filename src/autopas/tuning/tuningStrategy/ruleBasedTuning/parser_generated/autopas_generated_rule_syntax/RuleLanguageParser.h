
// Generated from RuleLanguage.g4 by ANTLR 4.13.2

#pragma once

#include <any>

#include "antlr4-runtime.h"

namespace AutopasGeneratedRuleSyntax {

class RuleLanguageParser : public antlr4::Parser {
 public:
  enum {
    T__0 = 1,
    T__1 = 2,
    T__2 = 3,
    T__3 = 4,
    T__4 = 5,
    T__5 = 6,
    T__6 = 7,
    T__7 = 8,
    T__8 = 9,
    T__9 = 10,
    T__10 = 11,
    T__11 = 12,
    T__12 = 13,
    T__13 = 14,
    T__14 = 15,
    T__15 = 16,
    T__16 = 17,
    T__17 = 18,
    T__18 = 19,
    T__19 = 20,
    T__20 = 21,
    T__21 = 22,
    T__22 = 23,
    T__23 = 24,
    COMMENT = 25,
    WS = 26,
    Container_opt = 27,
    Traversal_opt = 28,
    Load_estimator_opt = 29,
    Data_layout_opt = 30,
    Newton3_opt = 31,
    Bool_val = 32,
    Configuration_property = 33,
    DIGIT = 34,
    Unsigned_val = 35,
    Double_val = 36,
    Variable_name = 37
  };

  enum {
    RuleProgram = 0,
    RuleUnsigned_val = 1,
    RuleLiteral = 2,
    RuleDefine_list = 3,
    RuleVariable = 4,
    RuleExpression = 5,
    RuleDefine = 6,
    RuleProperty_value = 7,
    RuleConfiguration_pattern = 8,
    RuleConfiguration_order = 9,
    RuleStatement = 10,
    RuleIf_statement = 11
  };

  explicit RuleLanguageParser(antlr4::TokenStream *input);

  RuleLanguageParser(antlr4::TokenStream *input, const antlr4::atn::ParserATNSimulatorOptions &options);

  ~RuleLanguageParser() override;

  std::string getGrammarFileName() const override;

  const antlr4::atn::ATN &getATN() const override;

  const std::vector<std::string> &getRuleNames() const override;

  const antlr4::dfa::Vocabulary &getVocabulary() const override;

  antlr4::atn::SerializedATNView getSerializedATN() const override;

  class ProgramContext;
  class Unsigned_valContext;
  class LiteralContext;
  class Define_listContext;
  class VariableContext;
  class ExpressionContext;
  class DefineContext;
  class Property_valueContext;
  class Configuration_patternContext;
  class Configuration_orderContext;
  class StatementContext;
  class If_statementContext;

  class ProgramContext : public antlr4::ParserRuleContext {
   public:
    ProgramContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<StatementContext *> statement();
    StatementContext *statement(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  ProgramContext *program();

  class Unsigned_valContext : public antlr4::ParserRuleContext {
   public:
    Unsigned_valContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Unsigned_val();
    antlr4::tree::TerminalNode *DIGIT();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Unsigned_valContext *unsigned_val();

  class LiteralContext : public antlr4::ParserRuleContext {
   public:
    LiteralContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Traversal_opt();
    antlr4::tree::TerminalNode *Container_opt();
    antlr4::tree::TerminalNode *Load_estimator_opt();
    antlr4::tree::TerminalNode *Data_layout_opt();
    antlr4::tree::TerminalNode *Newton3_opt();
    Unsigned_valContext *unsigned_val();
    antlr4::tree::TerminalNode *Double_val();
    antlr4::tree::TerminalNode *Bool_val();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  LiteralContext *literal();

  class Define_listContext : public antlr4::ParserRuleContext {
   public:
    Define_listContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Variable_name();
    std::vector<LiteralContext *> literal();
    LiteralContext *literal(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Define_listContext *define_list();

  class VariableContext : public antlr4::ParserRuleContext {
   public:
    VariableContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Variable_name();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  VariableContext *variable();

  class ExpressionContext : public antlr4::ParserRuleContext {
   public:
    antlr4::Token *op = nullptr;
    ExpressionContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<ExpressionContext *> expression();
    ExpressionContext *expression(size_t i);
    LiteralContext *literal();
    VariableContext *variable();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  ExpressionContext *expression();
  ExpressionContext *expression(int precedence);
  class DefineContext : public antlr4::ParserRuleContext {
   public:
    DefineContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Variable_name();
    ExpressionContext *expression();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  DefineContext *define();

  class Property_valueContext : public antlr4::ParserRuleContext {
   public:
    Property_valueContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Variable_name();
    LiteralContext *literal();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Property_valueContext *property_value();

  class Configuration_patternContext : public antlr4::ParserRuleContext {
   public:
    Configuration_patternContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<antlr4::tree::TerminalNode *> Configuration_property();
    antlr4::tree::TerminalNode *Configuration_property(size_t i);
    std::vector<Property_valueContext *> property_value();
    Property_valueContext *property_value(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Configuration_patternContext *configuration_pattern();

  class Configuration_orderContext : public antlr4::ParserRuleContext {
   public:
    Configuration_orderContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<Configuration_patternContext *> configuration_pattern();
    Configuration_patternContext *configuration_pattern(size_t i);
    std::vector<antlr4::tree::TerminalNode *> Configuration_property();
    antlr4::tree::TerminalNode *Configuration_property(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Configuration_orderContext *configuration_order();

  class StatementContext : public antlr4::ParserRuleContext {
   public:
    StatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    Define_listContext *define_list();
    DefineContext *define();
    If_statementContext *if_statement();
    Configuration_orderContext *configuration_order();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  StatementContext *statement();

  class If_statementContext : public antlr4::ParserRuleContext {
   public:
    If_statementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    ExpressionContext *expression();
    std::vector<StatementContext *> statement();
    StatementContext *statement(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  If_statementContext *if_statement();

  bool sempred(antlr4::RuleContext *_localctx, size_t ruleIndex, size_t predicateIndex) override;

  bool expressionSempred(ExpressionContext *_localctx, size_t predicateIndex);

  // By default the static state used to implement the parser is lazily initialized during the first
  // call to the constructor. You can call this function if you wish to initialize the static state
  // ahead of time.
  static void initialize();

 private:
};

}  // namespace AutopasGeneratedRuleSyntax

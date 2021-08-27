
// Generated from /home/tobias/AutoPas2/src/autopas/selectors/tuningStrategy/ruleBasedTuning/RuleLanguage.g4 by ANTLR 4.9.1

#pragma once


#include "antlr4-runtime.h"


namespace autopas_generated_rule_syntax {


class  RuleLanguageParser : public antlr4::Parser {
public:
  enum {
    T__0 = 1, T__1 = 2, T__2 = 3, T__3 = 4, T__4 = 5, T__5 = 6, T__6 = 7, 
    T__7 = 8, T__8 = 9, T__9 = 10, T__10 = 11, T__11 = 12, T__12 = 13, T__13 = 14, 
    T__14 = 15, WS = 16, Container_opt = 17, Traversal_opt = 18, Load_estimator_opt = 19, 
    Data_layout_opt = 20, Newton3_opt = 21, Bool_val = 22, Configuration_property = 23, 
    DIGIT = 24, Unsigned_val = 25, Double_val = 26, Variable_name = 27
  };

  enum {
    RuleProgram = 0, RuleLiteral = 1, RuleDefine_list = 2, RuleDefine = 3, 
    RuleVariable = 4, RuleAtom_expr = 5, RuleComp_expr = 6, RuleExpression = 7, 
    RuleProperty_value = 8, RuleConfiguration_pattern = 9, RuleConfiguration_order = 10, 
    RuleStatement = 11, RuleIf_statement = 12
  };

  explicit RuleLanguageParser(antlr4::TokenStream *input);
  ~RuleLanguageParser();

  virtual std::string getGrammarFileName() const override;
  virtual const antlr4::atn::ATN& getATN() const override { return _atn; };
  virtual const std::vector<std::string>& getTokenNames() const override { return _tokenNames; }; // deprecated: use vocabulary instead.
  virtual const std::vector<std::string>& getRuleNames() const override;
  virtual antlr4::dfa::Vocabulary& getVocabulary() const override;


  class ProgramContext;
  class LiteralContext;
  class Define_listContext;
  class DefineContext;
  class VariableContext;
  class Atom_exprContext;
  class Comp_exprContext;
  class ExpressionContext;
  class Property_valueContext;
  class Configuration_patternContext;
  class Configuration_orderContext;
  class StatementContext;
  class If_statementContext; 

  class  ProgramContext : public antlr4::ParserRuleContext {
  public:
    ProgramContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<StatementContext *> statement();
    StatementContext* statement(size_t i);


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ProgramContext* program();

  class  LiteralContext : public antlr4::ParserRuleContext {
  public:
    LiteralContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Traversal_opt();
    antlr4::tree::TerminalNode *Container_opt();
    antlr4::tree::TerminalNode *Load_estimator_opt();
    antlr4::tree::TerminalNode *Data_layout_opt();
    antlr4::tree::TerminalNode *Newton3_opt();
    antlr4::tree::TerminalNode *Unsigned_val();
    antlr4::tree::TerminalNode *Double_val();
    antlr4::tree::TerminalNode *Bool_val();


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  LiteralContext* literal();

  class  Define_listContext : public antlr4::ParserRuleContext {
  public:
    Define_listContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Variable_name();
    std::vector<LiteralContext *> literal();
    LiteralContext* literal(size_t i);


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  Define_listContext* define_list();

  class  DefineContext : public antlr4::ParserRuleContext {
  public:
    DefineContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Variable_name();
    LiteralContext *literal();


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  DefineContext* define();

  class  VariableContext : public antlr4::ParserRuleContext {
  public:
    VariableContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Variable_name();


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  VariableContext* variable();

  class  Atom_exprContext : public antlr4::ParserRuleContext {
  public:
    Atom_exprContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    VariableContext *variable();
    LiteralContext *literal();


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  Atom_exprContext* atom_expr();

  class  Comp_exprContext : public antlr4::ParserRuleContext {
  public:
    antlr4::Token *op = nullptr;
    Comp_exprContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<Atom_exprContext *> atom_expr();
    Atom_exprContext* atom_expr(size_t i);


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  Comp_exprContext* comp_expr();

  class  ExpressionContext : public antlr4::ParserRuleContext {
  public:
    ExpressionContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<Comp_exprContext *> comp_expr();
    Comp_exprContext* comp_expr(size_t i);


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ExpressionContext* expression();

  class  Property_valueContext : public antlr4::ParserRuleContext {
  public:
    Property_valueContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *Variable_name();
    LiteralContext *literal();


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  Property_valueContext* property_value();

  class  Configuration_patternContext : public antlr4::ParserRuleContext {
  public:
    Configuration_patternContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<antlr4::tree::TerminalNode *> Configuration_property();
    antlr4::tree::TerminalNode* Configuration_property(size_t i);
    std::vector<Property_valueContext *> property_value();
    Property_valueContext* property_value(size_t i);


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  Configuration_patternContext* configuration_pattern();

  class  Configuration_orderContext : public antlr4::ParserRuleContext {
  public:
    Configuration_orderContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<Configuration_patternContext *> configuration_pattern();
    Configuration_patternContext* configuration_pattern(size_t i);


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  Configuration_orderContext* configuration_order();

  class  StatementContext : public antlr4::ParserRuleContext {
  public:
    StatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    Define_listContext *define_list();
    DefineContext *define();
    If_statementContext *if_statement();
    Configuration_orderContext *configuration_order();


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  StatementContext* statement();

  class  If_statementContext : public antlr4::ParserRuleContext {
  public:
    If_statementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    ExpressionContext *expression();
    std::vector<StatementContext *> statement();
    StatementContext* statement(size_t i);


    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  If_statementContext* if_statement();


private:
  static std::vector<antlr4::dfa::DFA> _decisionToDFA;
  static antlr4::atn::PredictionContextCache _sharedContextCache;
  static std::vector<std::string> _ruleNames;
  static std::vector<std::string> _tokenNames;

  static std::vector<std::string> _literalNames;
  static std::vector<std::string> _symbolicNames;
  static antlr4::dfa::Vocabulary _vocabulary;
  static antlr4::atn::ATN _atn;
  static std::vector<uint16_t> _serializedATN;


  struct Initializer {
    Initializer();
  };
  static Initializer _init;
};

}  // namespace autopas_generated_rule_syntax

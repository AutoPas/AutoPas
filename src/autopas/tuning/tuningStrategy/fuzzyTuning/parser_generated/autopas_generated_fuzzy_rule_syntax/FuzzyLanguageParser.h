
// Generated from
// autopas/tuning/tuningStrategy/fuzzyTuning\FuzzyLanguage.g4 by
// ANTLR 4.9.1

#pragma once

#include "antlr4-runtime.h"

namespace autopas_generated_fuzzy_rule_syntax {

class FuzzyLanguageParser : public antlr4::Parser {
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
    WS = 14,
    COMMENT = 15,
    STRING = 16,
    NUMBER = 17,
    IDENTIFIER = 18
  };

  enum {
    RuleRule_file = 0,
    RuleLinguistic_variable = 1,
    RuleFuzzy_term = 2,
    RuleFunction = 3,
    RuleFuzzy_rule = 4,
    RuleFuzzy_set = 5
  };

  explicit FuzzyLanguageParser(antlr4::TokenStream *input);
  ~FuzzyLanguageParser();

  virtual std::string getGrammarFileName() const override;
  virtual const antlr4::atn::ATN &getATN() const override { return _atn; };
  virtual const std::vector<std::string> &getTokenNames() const override {
    return _tokenNames;
  };  // deprecated: use vocabulary instead.
  virtual const std::vector<std::string> &getRuleNames() const override;
  virtual antlr4::dfa::Vocabulary &getVocabulary() const override;

  class Rule_fileContext;
  class Linguistic_variableContext;
  class Fuzzy_termContext;
  class FunctionContext;
  class Fuzzy_ruleContext;
  class Fuzzy_setContext;

  class Rule_fileContext : public antlr4::ParserRuleContext {
   public:
    Rule_fileContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *EOF();
    std::vector<Linguistic_variableContext *> linguistic_variable();
    Linguistic_variableContext *linguistic_variable(size_t i);
    std::vector<Fuzzy_ruleContext *> fuzzy_rule();
    Fuzzy_ruleContext *fuzzy_rule(size_t i);

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Rule_fileContext *rule_file();

  class Linguistic_variableContext : public antlr4::ParserRuleContext {
   public:
    Linguistic_variableContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *STRING();
    std::vector<antlr4::tree::TerminalNode *> NUMBER();
    antlr4::tree::TerminalNode *NUMBER(size_t i);
    std::vector<Fuzzy_termContext *> fuzzy_term();
    Fuzzy_termContext *fuzzy_term(size_t i);

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Linguistic_variableContext *linguistic_variable();

  class Fuzzy_termContext : public antlr4::ParserRuleContext {
   public:
    Fuzzy_termContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *STRING();
    FunctionContext *function();

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Fuzzy_termContext *fuzzy_term();

  class FunctionContext : public antlr4::ParserRuleContext {
   public:
    FunctionContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *IDENTIFIER();
    std::vector<antlr4::tree::TerminalNode *> NUMBER();
    antlr4::tree::TerminalNode *NUMBER(size_t i);

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  FunctionContext *function();

  class Fuzzy_ruleContext : public antlr4::ParserRuleContext {
   public:
    Fuzzy_ruleContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<Fuzzy_setContext *> fuzzy_set();
    Fuzzy_setContext *fuzzy_set(size_t i);

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Fuzzy_ruleContext *fuzzy_rule();

  class Fuzzy_setContext : public antlr4::ParserRuleContext {
   public:
    Fuzzy_setContext(antlr4::ParserRuleContext *parent, size_t invokingState);

    Fuzzy_setContext() = default;
    void copyFrom(Fuzzy_setContext *context);
    using antlr4::ParserRuleContext::copyFrom;

    virtual size_t getRuleIndex() const override;
  };

  class OrContext : public Fuzzy_setContext {
   public:
    OrContext(Fuzzy_setContext *ctx);

    std::vector<Fuzzy_setContext *> fuzzy_set();
    Fuzzy_setContext *fuzzy_set(size_t i);

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  class BracketsContext : public Fuzzy_setContext {
   public:
    BracketsContext(Fuzzy_setContext *ctx);

    Fuzzy_setContext *fuzzy_set();

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  class AndContext : public Fuzzy_setContext {
   public:
    AndContext(Fuzzy_setContext *ctx);

    std::vector<Fuzzy_setContext *> fuzzy_set();
    Fuzzy_setContext *fuzzy_set(size_t i);

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  class SelectContext : public Fuzzy_setContext {
   public:
    SelectContext(Fuzzy_setContext *ctx);

    std::vector<antlr4::tree::TerminalNode *> STRING();
    antlr4::tree::TerminalNode *STRING(size_t i);

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  class NegateContext : public Fuzzy_setContext {
   public:
    NegateContext(Fuzzy_setContext *ctx);

    Fuzzy_setContext *fuzzy_set();

    virtual antlrcpp::Any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Fuzzy_setContext *fuzzy_set();
  Fuzzy_setContext *fuzzy_set(int precedence);

  virtual bool sempred(antlr4::RuleContext *_localctx, size_t ruleIndex, size_t predicateIndex) override;
  bool fuzzy_setSempred(Fuzzy_setContext *_localctx, size_t predicateIndex);

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

}  // namespace autopas_generated_fuzzy_rule_syntax

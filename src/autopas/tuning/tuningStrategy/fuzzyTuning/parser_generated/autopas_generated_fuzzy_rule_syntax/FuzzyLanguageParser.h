
// Generated from FuzzyLanguage.g4 by ANTLR 4.13.2

#pragma once

#include <any>

#include "antlr4-runtime.h"

namespace AutopasGeneratedFuzzyRuleSyntax {

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
    T__13 = 14,
    T__14 = 15,
    T__15 = 16,
    T__16 = 17,
    T__17 = 18,
    T__18 = 19,
    WS = 20,
    COMMENT = 21,
    STRING = 22,
    NUMBER = 23,
    IDENTIFIER = 24
  };

  enum {
    RuleRule_file = 0,
    RuleSettings = 1,
    RuleLinguistic_variable = 2,
    RuleFuzzy_term = 3,
    RuleFunction = 4,
    RuleFuzzy_rule = 5,
    RuleFuzzy_set = 6,
    RuleOutput_mapping = 7,
    RuleOutput_entry = 8,
    RulePattern_mapping = 9,
    RuleConfiguration_pattern = 10
  };

  explicit FuzzyLanguageParser(antlr4::TokenStream *input);

  FuzzyLanguageParser(antlr4::TokenStream *input, const antlr4::atn::ParserATNSimulatorOptions &options);

  ~FuzzyLanguageParser() override;

  std::string getGrammarFileName() const override;

  const antlr4::atn::ATN &getATN() const override;

  const std::vector<std::string> &getRuleNames() const override;

  const antlr4::dfa::Vocabulary &getVocabulary() const override;

  antlr4::atn::SerializedATNView getSerializedATN() const override;

  class Rule_fileContext;
  class SettingsContext;
  class Linguistic_variableContext;
  class Fuzzy_termContext;
  class FunctionContext;
  class Fuzzy_ruleContext;
  class Fuzzy_setContext;
  class Output_mappingContext;
  class Output_entryContext;
  class Pattern_mappingContext;
  class Configuration_patternContext;

  class Rule_fileContext : public antlr4::ParserRuleContext {
   public:
    Rule_fileContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    SettingsContext *settings();
    Output_mappingContext *output_mapping();
    antlr4::tree::TerminalNode *EOF();
    std::vector<Linguistic_variableContext *> linguistic_variable();
    Linguistic_variableContext *linguistic_variable(size_t i);
    std::vector<Fuzzy_ruleContext *> fuzzy_rule();
    Fuzzy_ruleContext *fuzzy_rule(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Rule_fileContext *rule_file();

  class SettingsContext : public antlr4::ParserRuleContext {
   public:
    SettingsContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<antlr4::tree::TerminalNode *> IDENTIFIER();
    antlr4::tree::TerminalNode *IDENTIFIER(size_t i);
    std::vector<antlr4::tree::TerminalNode *> STRING();
    antlr4::tree::TerminalNode *STRING(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  SettingsContext *settings();

  class Linguistic_variableContext : public antlr4::ParserRuleContext {
   public:
    Linguistic_variableContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *STRING();
    std::vector<antlr4::tree::TerminalNode *> NUMBER();
    antlr4::tree::TerminalNode *NUMBER(size_t i);
    std::vector<Fuzzy_termContext *> fuzzy_term();
    Fuzzy_termContext *fuzzy_term(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Linguistic_variableContext *linguistic_variable();

  class Fuzzy_termContext : public antlr4::ParserRuleContext {
   public:
    Fuzzy_termContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *STRING();
    FunctionContext *function();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Fuzzy_termContext *fuzzy_term();

  class FunctionContext : public antlr4::ParserRuleContext {
   public:
    FunctionContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *IDENTIFIER();
    std::vector<antlr4::tree::TerminalNode *> NUMBER();
    antlr4::tree::TerminalNode *NUMBER(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  FunctionContext *function();

  class Fuzzy_ruleContext : public antlr4::ParserRuleContext {
   public:
    Fuzzy_ruleContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<Fuzzy_setContext *> fuzzy_set();
    Fuzzy_setContext *fuzzy_set(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
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

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  class BracketsContext : public Fuzzy_setContext {
   public:
    BracketsContext(Fuzzy_setContext *ctx);

    Fuzzy_setContext *fuzzy_set();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  class AndContext : public Fuzzy_setContext {
   public:
    AndContext(Fuzzy_setContext *ctx);

    std::vector<Fuzzy_setContext *> fuzzy_set();
    Fuzzy_setContext *fuzzy_set(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  class SelectContext : public Fuzzy_setContext {
   public:
    SelectContext(Fuzzy_setContext *ctx);

    std::vector<antlr4::tree::TerminalNode *> STRING();
    antlr4::tree::TerminalNode *STRING(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  class NegateContext : public Fuzzy_setContext {
   public:
    NegateContext(Fuzzy_setContext *ctx);

    Fuzzy_setContext *fuzzy_set();

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Fuzzy_setContext *fuzzy_set();
  Fuzzy_setContext *fuzzy_set(int precedence);
  class Output_mappingContext : public antlr4::ParserRuleContext {
   public:
    Output_mappingContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<Output_entryContext *> output_entry();
    Output_entryContext *output_entry(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Output_mappingContext *output_mapping();

  class Output_entryContext : public antlr4::ParserRuleContext {
   public:
    Output_entryContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *STRING();
    std::vector<Pattern_mappingContext *> pattern_mapping();
    Pattern_mappingContext *pattern_mapping(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Output_entryContext *output_entry();

  class Pattern_mappingContext : public antlr4::ParserRuleContext {
   public:
    Pattern_mappingContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *NUMBER();
    std::vector<Configuration_patternContext *> configuration_pattern();
    Configuration_patternContext *configuration_pattern(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Pattern_mappingContext *pattern_mapping();

  class Configuration_patternContext : public antlr4::ParserRuleContext {
   public:
    Configuration_patternContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<antlr4::tree::TerminalNode *> IDENTIFIER();
    antlr4::tree::TerminalNode *IDENTIFIER(size_t i);
    std::vector<antlr4::tree::TerminalNode *> STRING();
    antlr4::tree::TerminalNode *STRING(size_t i);

    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
  };

  Configuration_patternContext *configuration_pattern();

  bool sempred(antlr4::RuleContext *_localctx, size_t ruleIndex, size_t predicateIndex) override;

  bool fuzzy_setSempred(Fuzzy_setContext *_localctx, size_t predicateIndex);

  // By default the static state used to implement the parser is lazily initialized during the first
  // call to the constructor. You can call this function if you wish to initialize the static state
  // ahead of time.
  static void initialize();

 private:
};

}  // namespace AutopasGeneratedFuzzyRuleSyntax

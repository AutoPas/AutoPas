
// Generated from RuleLanguage.g4 by ANTLR 4.13.2

#pragma once

#include <any>

#include "antlr4-runtime.h"

namespace AutopasGeneratedRuleSyntax {

class RuleLanguageLexer : public antlr4::Lexer {
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

  explicit RuleLanguageLexer(antlr4::CharStream *input);

  ~RuleLanguageLexer() override;

  std::string getGrammarFileName() const override;

  const std::vector<std::string> &getRuleNames() const override;

  const std::vector<std::string> &getChannelNames() const override;

  const std::vector<std::string> &getModeNames() const override;

  const antlr4::dfa::Vocabulary &getVocabulary() const override;

  antlr4::atn::SerializedATNView getSerializedATN() const override;

  const antlr4::atn::ATN &getATN() const override;

  // By default the static state used to implement the lexer is lazily initialized during the first
  // call to the constructor. You can call this function if you wish to initialize the static state
  // ahead of time.
  static void initialize();

 private:
  // Individual action functions triggered by action() above.

  // Individual semantic predicate functions triggered by sempred() above.
};

}  // namespace AutopasGeneratedRuleSyntax

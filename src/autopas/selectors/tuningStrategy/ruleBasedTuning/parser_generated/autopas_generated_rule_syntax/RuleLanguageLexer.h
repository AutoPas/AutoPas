
// Generated from /home/tobias/AutoPas2/src/autopas/selectors/tuningStrategy/ruleBasedTuning/RuleLanguage.g4 by
// ANTLR 4.9.1

#pragma once

#include "antlr4-runtime.h"

namespace autopas_generated_rule_syntax {

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
    WS = 17,
    Container_opt = 18,
    Traversal_opt = 19,
    Load_estimator_opt = 20,
    Data_layout_opt = 21,
    Newton3_opt = 22,
    Bool_val = 23,
    Configuration_property = 24,
    DIGIT = 25,
    Unsigned_val = 26,
    Double_val = 27,
    Variable_name = 28
  };

  explicit RuleLanguageLexer(antlr4::CharStream *input);
  ~RuleLanguageLexer();

  virtual std::string getGrammarFileName() const override;
  virtual const std::vector<std::string> &getRuleNames() const override;

  virtual const std::vector<std::string> &getChannelNames() const override;
  virtual const std::vector<std::string> &getModeNames() const override;
  virtual const std::vector<std::string> &getTokenNames() const override;  // deprecated, use vocabulary instead
  virtual antlr4::dfa::Vocabulary &getVocabulary() const override;

  virtual const std::vector<uint16_t> getSerializedATN() const override;
  virtual const antlr4::atn::ATN &getATN() const override;

 private:
  static std::vector<antlr4::dfa::DFA> _decisionToDFA;
  static antlr4::atn::PredictionContextCache _sharedContextCache;
  static std::vector<std::string> _ruleNames;
  static std::vector<std::string> _tokenNames;
  static std::vector<std::string> _channelNames;
  static std::vector<std::string> _modeNames;

  static std::vector<std::string> _literalNames;
  static std::vector<std::string> _symbolicNames;
  static antlr4::dfa::Vocabulary _vocabulary;
  static antlr4::atn::ATN _atn;
  static std::vector<uint16_t> _serializedATN;

  // Individual action functions triggered by action() above.

  // Individual semantic predicate functions triggered by sempred() above.

  struct Initializer {
    Initializer();
  };
  static Initializer _init;
};

}  // namespace autopas_generated_rule_syntax

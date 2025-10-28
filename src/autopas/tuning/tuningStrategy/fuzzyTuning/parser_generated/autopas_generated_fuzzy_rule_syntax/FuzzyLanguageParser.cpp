
// Generated from FuzzyLanguage.g4 by ANTLR 4.13.2

#include "FuzzyLanguageParser.h"

#include "FuzzyLanguageVisitor.h"

using namespace antlrcpp;
using namespace AutopasGeneratedFuzzyRuleSyntax;

using namespace antlr4;

namespace {

struct FuzzyLanguageParserStaticData final {
  FuzzyLanguageParserStaticData(std::vector<std::string> ruleNames, std::vector<std::string> literalNames,
                                std::vector<std::string> symbolicNames)
      : ruleNames(std::move(ruleNames)),
        literalNames(std::move(literalNames)),
        symbolicNames(std::move(symbolicNames)),
        vocabulary(this->literalNames, this->symbolicNames) {}

  FuzzyLanguageParserStaticData(const FuzzyLanguageParserStaticData &) = delete;
  FuzzyLanguageParserStaticData(FuzzyLanguageParserStaticData &&) = delete;
  FuzzyLanguageParserStaticData &operator=(const FuzzyLanguageParserStaticData &) = delete;
  FuzzyLanguageParserStaticData &operator=(FuzzyLanguageParserStaticData &&) = delete;

  std::vector<antlr4::dfa::DFA> decisionToDFA;
  antlr4::atn::PredictionContextCache sharedContextCache;
  const std::vector<std::string> ruleNames;
  const std::vector<std::string> literalNames;
  const std::vector<std::string> symbolicNames;
  const antlr4::dfa::Vocabulary vocabulary;
  antlr4::atn::SerializedATNView serializedATN;
  std::unique_ptr<antlr4::atn::ATN> atn;
};

::antlr4::internal::OnceFlag fuzzylanguageParserOnceFlag;
#if ANTLR4_USE_THREAD_LOCAL_CACHE
static thread_local
#endif
    std::unique_ptr<FuzzyLanguageParserStaticData>
        fuzzylanguageParserStaticData = nullptr;

void fuzzylanguageParserInitialize() {
#if ANTLR4_USE_THREAD_LOCAL_CACHE
  if (fuzzylanguageParserStaticData != nullptr) {
    return;
  }
#else
  assert(fuzzylanguageParserStaticData == nullptr);
#endif
  auto staticData = std::make_unique<FuzzyLanguageParserStaticData>(
      std::vector<std::string>{"rule_file", "settings", "linguistic_variable", "fuzzy_term", "function", "fuzzy_rule",
                               "fuzzy_set", "output_mapping", "output_entry", "pattern_mapping",
                               "configuration_pattern"},
      std::vector<std::string>{"",         "'FuzzySystemSettings'",
                               "':'",      "'FuzzyVariable'",
                               "'domain'", "'range'",
                               "'('",      "','",
                               "')'",      "'if'",
                               "'then'",   "'&&'",
                               "'||'",     "'!'",
                               "'=='",     "'OutputMapping'",
                               "'=>'",     "'['",
                               "'='",      "']'"},
      std::vector<std::string>{"", "", "", "", "", "", "", "",   "",        "",       "",       "",          "",
                               "", "", "", "", "", "", "", "WS", "COMMENT", "STRING", "NUMBER", "IDENTIFIER"});
  static const int32_t serializedATNSegment[] = {
      4,  1,   24,  150, 2,   0,  7,   0,   2,   1,   7,   1,   2,   2,   7,   2,   2,   3,  7,  3,   2,   4,   7,   4,
      2,  5,   7,   5,   2,   6,  7,   6,   2,   7,   7,   7,   2,   8,   7,   8,   2,   9,  7,  9,   2,   10,  7,   10,
      1,  0,   1,   0,   5,   0,  25,  8,   0,   10,  0,   12,  0,   28,  9,   0,   1,   0,  1,  0,   5,   0,   32,  8,
      0,  10,  0,   12,  0,   35, 9,   0,   1,   0,   1,   0,   1,   1,   1,   1,   1,   1,  1,  1,   1,   1,   5,   1,
      44, 8,   1,   10,  1,   12, 1,   47,  9,   1,   1,   2,   1,   2,   1,   2,   1,   2,  1,  2,   1,   2,   1,   2,
      1,  2,   1,   2,   1,   2,  1,   2,   1,   2,   1,   2,   4,   2,   62,  8,   2,   11, 2,  12,  2,   63,  1,   3,
      1,  3,   1,   3,   1,   3,  1,   4,   1,   4,   1,   4,   1,   4,   1,   4,   5,   4,  75, 8,   4,   10,  4,   12,
      4,  78,  9,   4,   1,   4,  1,   4,   1,   5,   1,   5,   1,   5,   1,   5,   1,   5,  1,  6,   1,   6,   1,   6,
      1,  6,   1,   6,   1,   6,  1,   6,   1,   6,   1,   6,   1,   6,   3,   6,   97,  8,  6,  1,   6,   1,   6,   1,
      6,  1,   6,   1,   6,   1,  6,   5,   6,   105, 8,   6,   10,  6,   12,  6,   108, 9,  6,  1,   7,   1,   7,   1,
      7,  4,   7,   113, 8,   7,  11,  7,   12,  7,   114, 1,   8,   1,   8,   1,   8,   4,  8,  120, 8,   8,   11,  8,
      12, 8,   121, 1,   9,   1,  9,   1,   9,   1,   9,   1,   9,   5,   9,   129, 8,   9,  10, 9,   12,  9,   132, 9,
      9,  1,   10,  1,   10,  1,  10,  1,   10,  1,   10,  1,   10,  1,   10,  1,   10,  1,  10, 5,   10,  143, 8,   10,
      10, 10,  12,  10,  146, 9,  10,  1,   10,  1,   10,  1,   10,  0,   1,   12,  11,  0,  2,  4,   6,   8,   10,  12,
      14, 16,  18,  20,  0,   0,  151, 0,   22,  1,   0,   0,   0,   2,   38,  1,   0,   0,  0,  4,   48,  1,   0,   0,
      0,  6,   65,  1,   0,   0,  0,   8,   69,  1,   0,   0,   0,   10,  81,  1,   0,   0,  0,  12,  96,  1,   0,   0,
      0,  14,  109, 1,   0,   0,  0,   16,  116, 1,   0,   0,   0,   18,  123, 1,   0,   0,  0,  20,  133, 1,   0,   0,
      0,  22,  26,  3,   2,   1,  0,   23,  25,  3,   4,   2,   0,   24,  23,  1,   0,   0,  0,  25,  28,  1,   0,   0,
      0,  26,  24,  1,   0,   0,  0,   26,  27,  1,   0,   0,   0,   27,  29,  1,   0,   0,  0,  28,  26,  1,   0,   0,
      0,  29,  33,  3,   14,  7,  0,   30,  32,  3,   10,  5,   0,   31,  30,  1,   0,   0,  0,  32,  35,  1,   0,   0,
      0,  33,  31,  1,   0,   0,  0,   33,  34,  1,   0,   0,   0,   34,  36,  1,   0,   0,  0,  35,  33,  1,   0,   0,
      0,  36,  37,  5,   0,   0,  1,   37,  1,   1,   0,   0,   0,   38,  39,  5,   1,   0,  0,  39,  45,  5,   2,   0,
      0,  40,  41,  5,   24,  0,  0,   41,  42,  5,   2,   0,   0,   42,  44,  5,   22,  0,  0,  43,  40,  1,   0,   0,
      0,  44,  47,  1,   0,   0,  0,   45,  43,  1,   0,   0,   0,   45,  46,  1,   0,   0,  0,  46,  3,   1,   0,   0,
      0,  47,  45,  1,   0,   0,  0,   48,  49,  5,   3,   0,   0,   49,  50,  5,   2,   0,  0,  50,  51,  5,   4,   0,
      0,  51,  52,  5,   2,   0,  0,   52,  53,  5,   22,  0,   0,   53,  54,  5,   5,   0,  0,  54,  55,  5,   2,   0,
      0,  55,  56,  5,   6,   0,  0,   56,  57,  5,   23,  0,   0,   57,  58,  5,   7,   0,  0,  58,  59,  5,   23,  0,
      0,  59,  61,  5,   8,   0,  0,   60,  62,  3,   6,   3,   0,   61,  60,  1,   0,   0,  0,  62,  63,  1,   0,   0,
      0,  63,  61,  1,   0,   0,  0,   63,  64,  1,   0,   0,   0,   64,  5,   1,   0,   0,  0,  65,  66,  5,   22,  0,
      0,  66,  67,  5,   2,   0,  0,   67,  68,  3,   8,   4,   0,   68,  7,   1,   0,   0,  0,  69,  70,  5,   24,  0,
      0,  70,  71,  5,   6,   0,  0,   71,  76,  5,   23,  0,   0,   72,  73,  5,   7,   0,  0,  73,  75,  5,   23,  0,
      0,  74,  72,  1,   0,   0,  0,   75,  78,  1,   0,   0,   0,   76,  74,  1,   0,   0,  0,  76,  77,  1,   0,   0,
      0,  77,  79,  1,   0,   0,  0,   78,  76,  1,   0,   0,   0,   79,  80,  5,   8,   0,  0,  80,  9,   1,   0,   0,
      0,  81,  82,  5,   9,   0,  0,   82,  83,  3,   12,  6,   0,   83,  84,  5,   10,  0,  0,  84,  85,  3,   12,  6,
      0,  85,  11,  1,   0,   0,  0,   86,  87,  6,   6,   -1,  0,   87,  88,  5,   6,   0,  0,  88,  89,  3,   12,  6,
      0,  89,  90,  5,   8,   0,  0,   90,  97,  1,   0,   0,   0,   91,  92,  5,   13,  0,  0,  92,  97,  3,   12,  6,
      2,  93,  94,  5,   22,  0,  0,   94,  95,  5,   14,  0,   0,   95,  97,  5,   22,  0,  0,  96,  86,  1,   0,   0,
      0,  96,  91,  1,   0,   0,  0,   96,  93,  1,   0,   0,   0,   97,  106, 1,   0,   0,  0,  98,  99,  10,  4,   0,
      0,  99,  100, 5,   11,  0,  0,   100, 105, 3,   12,  6,   5,   101, 102, 10,  3,   0,  0,  102, 103, 5,   12,  0,
      0,  103, 105, 3,   12,  6,  4,   104, 98,  1,   0,   0,   0,   104, 101, 1,   0,   0,  0,  105, 108, 1,   0,   0,
      0,  106, 104, 1,   0,   0,  0,   106, 107, 1,   0,   0,   0,   107, 13,  1,   0,   0,  0,  108, 106, 1,   0,   0,
      0,  109, 110, 5,   15,  0,  0,   110, 112, 5,   2,   0,   0,   111, 113, 3,   16,  8,  0,  112, 111, 1,   0,   0,
      0,  113, 114, 1,   0,   0,  0,   114, 112, 1,   0,   0,   0,   114, 115, 1,   0,   0,  0,  115, 15,  1,   0,   0,
      0,  116, 117, 5,   22,  0,  0,   117, 119, 5,   2,   0,   0,   118, 120, 3,   18,  9,  0,  119, 118, 1,   0,   0,
      0,  120, 121, 1,   0,   0,  0,   121, 119, 1,   0,   0,   0,   121, 122, 1,   0,   0,  0,  122, 17,  1,   0,   0,
      0,  123, 124, 5,   23,  0,  0,   124, 125, 5,   16,  0,   0,   125, 130, 3,   20,  10, 0,  126, 127, 5,   7,   0,
      0,  127, 129, 3,   20,  10, 0,   128, 126, 1,   0,   0,   0,   129, 132, 1,   0,   0,  0,  130, 128, 1,   0,   0,
      0,  130, 131, 1,   0,   0,  0,   131, 19,  1,   0,   0,   0,   132, 130, 1,   0,   0,  0,  133, 134, 5,   17,  0,
      0,  134, 135, 5,   24,  0,  0,   135, 136, 5,   18,  0,   0,   136, 137, 5,   22,  0,  0,  137, 144, 1,   0,   0,
      0,  138, 139, 5,   7,   0,  0,   139, 140, 5,   24,  0,   0,   140, 141, 5,   18,  0,  0,  141, 143, 5,   22,  0,
      0,  142, 138, 1,   0,   0,  0,   143, 146, 1,   0,   0,   0,   144, 142, 1,   0,   0,  0,  144, 145, 1,   0,   0,
      0,  145, 147, 1,   0,   0,  0,   146, 144, 1,   0,   0,   0,   147, 148, 5,   19,  0,  0,  148, 21,  1,   0,   0,
      0,  12,  26,  33,  45,  63, 76,  96,  104, 106, 114, 121, 130, 144};
  staticData->serializedATN = antlr4::atn::SerializedATNView(
      serializedATNSegment, sizeof(serializedATNSegment) / sizeof(serializedATNSegment[0]));

  antlr4::atn::ATNDeserializer deserializer;
  staticData->atn = deserializer.deserialize(staticData->serializedATN);

  const size_t count = staticData->atn->getNumberOfDecisions();
  staticData->decisionToDFA.reserve(count);
  for (size_t i = 0; i < count; i++) {
    staticData->decisionToDFA.emplace_back(staticData->atn->getDecisionState(i), i);
  }
  fuzzylanguageParserStaticData = std::move(staticData);
}

}  // namespace

FuzzyLanguageParser::FuzzyLanguageParser(TokenStream *input)
    : FuzzyLanguageParser(input, antlr4::atn::ParserATNSimulatorOptions()) {}

FuzzyLanguageParser::FuzzyLanguageParser(TokenStream *input, const antlr4::atn::ParserATNSimulatorOptions &options)
    : Parser(input) {
  FuzzyLanguageParser::initialize();
  _interpreter = new atn::ParserATNSimulator(this, *fuzzylanguageParserStaticData->atn,
                                             fuzzylanguageParserStaticData->decisionToDFA,
                                             fuzzylanguageParserStaticData->sharedContextCache, options);
}

FuzzyLanguageParser::~FuzzyLanguageParser() { delete _interpreter; }

const atn::ATN &FuzzyLanguageParser::getATN() const { return *fuzzylanguageParserStaticData->atn; }

std::string FuzzyLanguageParser::getGrammarFileName() const { return "FuzzyLanguage.g4"; }

const std::vector<std::string> &FuzzyLanguageParser::getRuleNames() const {
  return fuzzylanguageParserStaticData->ruleNames;
}

const dfa::Vocabulary &FuzzyLanguageParser::getVocabulary() const { return fuzzylanguageParserStaticData->vocabulary; }

antlr4::atn::SerializedATNView FuzzyLanguageParser::getSerializedATN() const {
  return fuzzylanguageParserStaticData->serializedATN;
}

//----------------- Rule_fileContext ------------------------------------------------------------------

FuzzyLanguageParser::Rule_fileContext::Rule_fileContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

FuzzyLanguageParser::SettingsContext *FuzzyLanguageParser::Rule_fileContext::settings() {
  return getRuleContext<FuzzyLanguageParser::SettingsContext>(0);
}

FuzzyLanguageParser::Output_mappingContext *FuzzyLanguageParser::Rule_fileContext::output_mapping() {
  return getRuleContext<FuzzyLanguageParser::Output_mappingContext>(0);
}

tree::TerminalNode *FuzzyLanguageParser::Rule_fileContext::EOF() { return getToken(FuzzyLanguageParser::EOF, 0); }

std::vector<FuzzyLanguageParser::Linguistic_variableContext *>
FuzzyLanguageParser::Rule_fileContext::linguistic_variable() {
  return getRuleContexts<FuzzyLanguageParser::Linguistic_variableContext>();
}

FuzzyLanguageParser::Linguistic_variableContext *FuzzyLanguageParser::Rule_fileContext::linguistic_variable(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Linguistic_variableContext>(i);
}

std::vector<FuzzyLanguageParser::Fuzzy_ruleContext *> FuzzyLanguageParser::Rule_fileContext::fuzzy_rule() {
  return getRuleContexts<FuzzyLanguageParser::Fuzzy_ruleContext>();
}

FuzzyLanguageParser::Fuzzy_ruleContext *FuzzyLanguageParser::Rule_fileContext::fuzzy_rule(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_ruleContext>(i);
}

size_t FuzzyLanguageParser::Rule_fileContext::getRuleIndex() const { return FuzzyLanguageParser::RuleRule_file; }

std::any FuzzyLanguageParser::Rule_fileContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitRule_file(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Rule_fileContext *FuzzyLanguageParser::rule_file() {
  Rule_fileContext *_localctx = _tracker.createInstance<Rule_fileContext>(_ctx, getState());
  enterRule(_localctx, 0, FuzzyLanguageParser::RuleRule_file);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(22);
    settings();
    setState(26);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__2) {
      setState(23);
      linguistic_variable();
      setState(28);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(29);
    output_mapping();
    setState(33);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__8) {
      setState(30);
      fuzzy_rule();
      setState(35);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(36);
    match(FuzzyLanguageParser::EOF);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- SettingsContext ------------------------------------------------------------------

FuzzyLanguageParser::SettingsContext::SettingsContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<tree::TerminalNode *> FuzzyLanguageParser::SettingsContext::IDENTIFIER() {
  return getTokens(FuzzyLanguageParser::IDENTIFIER);
}

tree::TerminalNode *FuzzyLanguageParser::SettingsContext::IDENTIFIER(size_t i) {
  return getToken(FuzzyLanguageParser::IDENTIFIER, i);
}

std::vector<tree::TerminalNode *> FuzzyLanguageParser::SettingsContext::STRING() {
  return getTokens(FuzzyLanguageParser::STRING);
}

tree::TerminalNode *FuzzyLanguageParser::SettingsContext::STRING(size_t i) {
  return getToken(FuzzyLanguageParser::STRING, i);
}

size_t FuzzyLanguageParser::SettingsContext::getRuleIndex() const { return FuzzyLanguageParser::RuleSettings; }

std::any FuzzyLanguageParser::SettingsContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitSettings(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::SettingsContext *FuzzyLanguageParser::settings() {
  SettingsContext *_localctx = _tracker.createInstance<SettingsContext>(_ctx, getState());
  enterRule(_localctx, 2, FuzzyLanguageParser::RuleSettings);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(38);
    match(FuzzyLanguageParser::T__0);
    setState(39);
    match(FuzzyLanguageParser::T__1);
    setState(45);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::IDENTIFIER) {
      setState(40);
      match(FuzzyLanguageParser::IDENTIFIER);
      setState(41);
      match(FuzzyLanguageParser::T__1);
      setState(42);
      match(FuzzyLanguageParser::STRING);
      setState(47);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Linguistic_variableContext ------------------------------------------------------------------

FuzzyLanguageParser::Linguistic_variableContext::Linguistic_variableContext(ParserRuleContext *parent,
                                                                            size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *FuzzyLanguageParser::Linguistic_variableContext::STRING() {
  return getToken(FuzzyLanguageParser::STRING, 0);
}

std::vector<tree::TerminalNode *> FuzzyLanguageParser::Linguistic_variableContext::NUMBER() {
  return getTokens(FuzzyLanguageParser::NUMBER);
}

tree::TerminalNode *FuzzyLanguageParser::Linguistic_variableContext::NUMBER(size_t i) {
  return getToken(FuzzyLanguageParser::NUMBER, i);
}

std::vector<FuzzyLanguageParser::Fuzzy_termContext *> FuzzyLanguageParser::Linguistic_variableContext::fuzzy_term() {
  return getRuleContexts<FuzzyLanguageParser::Fuzzy_termContext>();
}

FuzzyLanguageParser::Fuzzy_termContext *FuzzyLanguageParser::Linguistic_variableContext::fuzzy_term(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_termContext>(i);
}

size_t FuzzyLanguageParser::Linguistic_variableContext::getRuleIndex() const {
  return FuzzyLanguageParser::RuleLinguistic_variable;
}

std::any FuzzyLanguageParser::Linguistic_variableContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitLinguistic_variable(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Linguistic_variableContext *FuzzyLanguageParser::linguistic_variable() {
  Linguistic_variableContext *_localctx = _tracker.createInstance<Linguistic_variableContext>(_ctx, getState());
  enterRule(_localctx, 4, FuzzyLanguageParser::RuleLinguistic_variable);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(48);
    match(FuzzyLanguageParser::T__2);
    setState(49);
    match(FuzzyLanguageParser::T__1);
    setState(50);
    match(FuzzyLanguageParser::T__3);
    setState(51);
    match(FuzzyLanguageParser::T__1);
    setState(52);
    match(FuzzyLanguageParser::STRING);
    setState(53);
    match(FuzzyLanguageParser::T__4);
    setState(54);
    match(FuzzyLanguageParser::T__1);
    setState(55);
    match(FuzzyLanguageParser::T__5);
    setState(56);
    match(FuzzyLanguageParser::NUMBER);
    setState(57);
    match(FuzzyLanguageParser::T__6);
    setState(58);
    match(FuzzyLanguageParser::NUMBER);
    setState(59);
    match(FuzzyLanguageParser::T__7);
    setState(61);
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(60);
      fuzzy_term();
      setState(63);
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while (_la == FuzzyLanguageParser::STRING);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Fuzzy_termContext ------------------------------------------------------------------

FuzzyLanguageParser::Fuzzy_termContext::Fuzzy_termContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *FuzzyLanguageParser::Fuzzy_termContext::STRING() {
  return getToken(FuzzyLanguageParser::STRING, 0);
}

FuzzyLanguageParser::FunctionContext *FuzzyLanguageParser::Fuzzy_termContext::function() {
  return getRuleContext<FuzzyLanguageParser::FunctionContext>(0);
}

size_t FuzzyLanguageParser::Fuzzy_termContext::getRuleIndex() const { return FuzzyLanguageParser::RuleFuzzy_term; }

std::any FuzzyLanguageParser::Fuzzy_termContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitFuzzy_term(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Fuzzy_termContext *FuzzyLanguageParser::fuzzy_term() {
  Fuzzy_termContext *_localctx = _tracker.createInstance<Fuzzy_termContext>(_ctx, getState());
  enterRule(_localctx, 6, FuzzyLanguageParser::RuleFuzzy_term);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(65);
    match(FuzzyLanguageParser::STRING);
    setState(66);
    match(FuzzyLanguageParser::T__1);
    setState(67);
    function();

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- FunctionContext ------------------------------------------------------------------

FuzzyLanguageParser::FunctionContext::FunctionContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *FuzzyLanguageParser::FunctionContext::IDENTIFIER() {
  return getToken(FuzzyLanguageParser::IDENTIFIER, 0);
}

std::vector<tree::TerminalNode *> FuzzyLanguageParser::FunctionContext::NUMBER() {
  return getTokens(FuzzyLanguageParser::NUMBER);
}

tree::TerminalNode *FuzzyLanguageParser::FunctionContext::NUMBER(size_t i) {
  return getToken(FuzzyLanguageParser::NUMBER, i);
}

size_t FuzzyLanguageParser::FunctionContext::getRuleIndex() const { return FuzzyLanguageParser::RuleFunction; }

std::any FuzzyLanguageParser::FunctionContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitFunction(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::FunctionContext *FuzzyLanguageParser::function() {
  FunctionContext *_localctx = _tracker.createInstance<FunctionContext>(_ctx, getState());
  enterRule(_localctx, 8, FuzzyLanguageParser::RuleFunction);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(69);
    match(FuzzyLanguageParser::IDENTIFIER);
    setState(70);
    match(FuzzyLanguageParser::T__5);
    setState(71);
    match(FuzzyLanguageParser::NUMBER);
    setState(76);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__6) {
      setState(72);
      match(FuzzyLanguageParser::T__6);
      setState(73);
      match(FuzzyLanguageParser::NUMBER);
      setState(78);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(79);
    match(FuzzyLanguageParser::T__7);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Fuzzy_ruleContext ------------------------------------------------------------------

FuzzyLanguageParser::Fuzzy_ruleContext::Fuzzy_ruleContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<FuzzyLanguageParser::Fuzzy_setContext *> FuzzyLanguageParser::Fuzzy_ruleContext::fuzzy_set() {
  return getRuleContexts<FuzzyLanguageParser::Fuzzy_setContext>();
}

FuzzyLanguageParser::Fuzzy_setContext *FuzzyLanguageParser::Fuzzy_ruleContext::fuzzy_set(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_setContext>(i);
}

size_t FuzzyLanguageParser::Fuzzy_ruleContext::getRuleIndex() const { return FuzzyLanguageParser::RuleFuzzy_rule; }

std::any FuzzyLanguageParser::Fuzzy_ruleContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitFuzzy_rule(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Fuzzy_ruleContext *FuzzyLanguageParser::fuzzy_rule() {
  Fuzzy_ruleContext *_localctx = _tracker.createInstance<Fuzzy_ruleContext>(_ctx, getState());
  enterRule(_localctx, 10, FuzzyLanguageParser::RuleFuzzy_rule);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(81);
    match(FuzzyLanguageParser::T__8);
    setState(82);
    fuzzy_set(0);
    setState(83);
    match(FuzzyLanguageParser::T__9);
    setState(84);
    fuzzy_set(0);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Fuzzy_setContext ------------------------------------------------------------------

FuzzyLanguageParser::Fuzzy_setContext::Fuzzy_setContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

size_t FuzzyLanguageParser::Fuzzy_setContext::getRuleIndex() const { return FuzzyLanguageParser::RuleFuzzy_set; }

void FuzzyLanguageParser::Fuzzy_setContext::copyFrom(Fuzzy_setContext *ctx) { ParserRuleContext::copyFrom(ctx); }

//----------------- OrContext ------------------------------------------------------------------

std::vector<FuzzyLanguageParser::Fuzzy_setContext *> FuzzyLanguageParser::OrContext::fuzzy_set() {
  return getRuleContexts<FuzzyLanguageParser::Fuzzy_setContext>();
}

FuzzyLanguageParser::Fuzzy_setContext *FuzzyLanguageParser::OrContext::fuzzy_set(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_setContext>(i);
}

FuzzyLanguageParser::OrContext::OrContext(Fuzzy_setContext *ctx) { copyFrom(ctx); }

std::any FuzzyLanguageParser::OrContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitOr(this);
  else
    return visitor->visitChildren(this);
}
//----------------- BracketsContext ------------------------------------------------------------------

FuzzyLanguageParser::Fuzzy_setContext *FuzzyLanguageParser::BracketsContext::fuzzy_set() {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_setContext>(0);
}

FuzzyLanguageParser::BracketsContext::BracketsContext(Fuzzy_setContext *ctx) { copyFrom(ctx); }

std::any FuzzyLanguageParser::BracketsContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitBrackets(this);
  else
    return visitor->visitChildren(this);
}
//----------------- AndContext ------------------------------------------------------------------

std::vector<FuzzyLanguageParser::Fuzzy_setContext *> FuzzyLanguageParser::AndContext::fuzzy_set() {
  return getRuleContexts<FuzzyLanguageParser::Fuzzy_setContext>();
}

FuzzyLanguageParser::Fuzzy_setContext *FuzzyLanguageParser::AndContext::fuzzy_set(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_setContext>(i);
}

FuzzyLanguageParser::AndContext::AndContext(Fuzzy_setContext *ctx) { copyFrom(ctx); }

std::any FuzzyLanguageParser::AndContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitAnd(this);
  else
    return visitor->visitChildren(this);
}
//----------------- SelectContext ------------------------------------------------------------------

std::vector<tree::TerminalNode *> FuzzyLanguageParser::SelectContext::STRING() {
  return getTokens(FuzzyLanguageParser::STRING);
}

tree::TerminalNode *FuzzyLanguageParser::SelectContext::STRING(size_t i) {
  return getToken(FuzzyLanguageParser::STRING, i);
}

FuzzyLanguageParser::SelectContext::SelectContext(Fuzzy_setContext *ctx) { copyFrom(ctx); }

std::any FuzzyLanguageParser::SelectContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitSelect(this);
  else
    return visitor->visitChildren(this);
}
//----------------- NegateContext ------------------------------------------------------------------

FuzzyLanguageParser::Fuzzy_setContext *FuzzyLanguageParser::NegateContext::fuzzy_set() {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_setContext>(0);
}

FuzzyLanguageParser::NegateContext::NegateContext(Fuzzy_setContext *ctx) { copyFrom(ctx); }

std::any FuzzyLanguageParser::NegateContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitNegate(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Fuzzy_setContext *FuzzyLanguageParser::fuzzy_set() { return fuzzy_set(0); }

FuzzyLanguageParser::Fuzzy_setContext *FuzzyLanguageParser::fuzzy_set(int precedence) {
  ParserRuleContext *parentContext = _ctx;
  size_t parentState = getState();
  FuzzyLanguageParser::Fuzzy_setContext *_localctx = _tracker.createInstance<Fuzzy_setContext>(_ctx, parentState);
  FuzzyLanguageParser::Fuzzy_setContext *previousContext = _localctx;
  (void)previousContext;  // Silence compiler, in case the context is not used by generated code.
  size_t startState = 12;
  enterRecursionRule(_localctx, 12, FuzzyLanguageParser::RuleFuzzy_set, precedence);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    unrollRecursionContexts(parentContext);
  });
  try {
    size_t alt;
    enterOuterAlt(_localctx, 1);
    setState(96);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case FuzzyLanguageParser::T__5: {
        _localctx = _tracker.createInstance<BracketsContext>(_localctx);
        _ctx = _localctx;
        previousContext = _localctx;

        setState(87);
        match(FuzzyLanguageParser::T__5);
        setState(88);
        fuzzy_set(0);
        setState(89);
        match(FuzzyLanguageParser::T__7);
        break;
      }

      case FuzzyLanguageParser::T__12: {
        _localctx = _tracker.createInstance<NegateContext>(_localctx);
        _ctx = _localctx;
        previousContext = _localctx;
        setState(91);
        match(FuzzyLanguageParser::T__12);
        setState(92);
        fuzzy_set(2);
        break;
      }

      case FuzzyLanguageParser::STRING: {
        _localctx = _tracker.createInstance<SelectContext>(_localctx);
        _ctx = _localctx;
        previousContext = _localctx;
        setState(93);
        match(FuzzyLanguageParser::STRING);
        setState(94);
        match(FuzzyLanguageParser::T__13);
        setState(95);
        match(FuzzyLanguageParser::STRING);
        break;
      }

      default:
        throw NoViableAltException(this);
    }
    _ctx->stop = _input->LT(-1);
    setState(106);
    _errHandler->sync(this);
    alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 7, _ctx);
    while (alt != 2 && alt != atn::ATN::INVALID_ALT_NUMBER) {
      if (alt == 1) {
        if (!_parseListeners.empty()) triggerExitRuleEvent();
        previousContext = _localctx;
        setState(104);
        _errHandler->sync(this);
        switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 6, _ctx)) {
          case 1: {
            auto newContext = _tracker.createInstance<AndContext>(
                _tracker.createInstance<Fuzzy_setContext>(parentContext, parentState));
            _localctx = newContext;
            pushNewRecursionContext(newContext, startState, RuleFuzzy_set);
            setState(98);

            if (!(precpred(_ctx, 4))) throw FailedPredicateException(this, "precpred(_ctx, 4)");
            setState(99);
            match(FuzzyLanguageParser::T__10);
            setState(100);
            fuzzy_set(5);
            break;
          }

          case 2: {
            auto newContext = _tracker.createInstance<OrContext>(
                _tracker.createInstance<Fuzzy_setContext>(parentContext, parentState));
            _localctx = newContext;
            pushNewRecursionContext(newContext, startState, RuleFuzzy_set);
            setState(101);

            if (!(precpred(_ctx, 3))) throw FailedPredicateException(this, "precpred(_ctx, 3)");
            setState(102);
            match(FuzzyLanguageParser::T__11);
            setState(103);
            fuzzy_set(4);
            break;
          }

          default:
            break;
        }
      }
      setState(108);
      _errHandler->sync(this);
      alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 7, _ctx);
    }
  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }
  return _localctx;
}

//----------------- Output_mappingContext ------------------------------------------------------------------

FuzzyLanguageParser::Output_mappingContext::Output_mappingContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<FuzzyLanguageParser::Output_entryContext *> FuzzyLanguageParser::Output_mappingContext::output_entry() {
  return getRuleContexts<FuzzyLanguageParser::Output_entryContext>();
}

FuzzyLanguageParser::Output_entryContext *FuzzyLanguageParser::Output_mappingContext::output_entry(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Output_entryContext>(i);
}

size_t FuzzyLanguageParser::Output_mappingContext::getRuleIndex() const {
  return FuzzyLanguageParser::RuleOutput_mapping;
}

std::any FuzzyLanguageParser::Output_mappingContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitOutput_mapping(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Output_mappingContext *FuzzyLanguageParser::output_mapping() {
  Output_mappingContext *_localctx = _tracker.createInstance<Output_mappingContext>(_ctx, getState());
  enterRule(_localctx, 14, FuzzyLanguageParser::RuleOutput_mapping);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(109);
    match(FuzzyLanguageParser::T__14);
    setState(110);
    match(FuzzyLanguageParser::T__1);
    setState(112);
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(111);
      output_entry();
      setState(114);
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while (_la == FuzzyLanguageParser::STRING);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Output_entryContext ------------------------------------------------------------------

FuzzyLanguageParser::Output_entryContext::Output_entryContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *FuzzyLanguageParser::Output_entryContext::STRING() {
  return getToken(FuzzyLanguageParser::STRING, 0);
}

std::vector<FuzzyLanguageParser::Pattern_mappingContext *> FuzzyLanguageParser::Output_entryContext::pattern_mapping() {
  return getRuleContexts<FuzzyLanguageParser::Pattern_mappingContext>();
}

FuzzyLanguageParser::Pattern_mappingContext *FuzzyLanguageParser::Output_entryContext::pattern_mapping(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Pattern_mappingContext>(i);
}

size_t FuzzyLanguageParser::Output_entryContext::getRuleIndex() const { return FuzzyLanguageParser::RuleOutput_entry; }

std::any FuzzyLanguageParser::Output_entryContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitOutput_entry(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Output_entryContext *FuzzyLanguageParser::output_entry() {
  Output_entryContext *_localctx = _tracker.createInstance<Output_entryContext>(_ctx, getState());
  enterRule(_localctx, 16, FuzzyLanguageParser::RuleOutput_entry);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(116);
    match(FuzzyLanguageParser::STRING);
    setState(117);
    match(FuzzyLanguageParser::T__1);
    setState(119);
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(118);
      pattern_mapping();
      setState(121);
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while (_la == FuzzyLanguageParser::NUMBER);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Pattern_mappingContext ------------------------------------------------------------------

FuzzyLanguageParser::Pattern_mappingContext::Pattern_mappingContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *FuzzyLanguageParser::Pattern_mappingContext::NUMBER() {
  return getToken(FuzzyLanguageParser::NUMBER, 0);
}

std::vector<FuzzyLanguageParser::Configuration_patternContext *>
FuzzyLanguageParser::Pattern_mappingContext::configuration_pattern() {
  return getRuleContexts<FuzzyLanguageParser::Configuration_patternContext>();
}

FuzzyLanguageParser::Configuration_patternContext *FuzzyLanguageParser::Pattern_mappingContext::configuration_pattern(
    size_t i) {
  return getRuleContext<FuzzyLanguageParser::Configuration_patternContext>(i);
}

size_t FuzzyLanguageParser::Pattern_mappingContext::getRuleIndex() const {
  return FuzzyLanguageParser::RulePattern_mapping;
}

std::any FuzzyLanguageParser::Pattern_mappingContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitPattern_mapping(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Pattern_mappingContext *FuzzyLanguageParser::pattern_mapping() {
  Pattern_mappingContext *_localctx = _tracker.createInstance<Pattern_mappingContext>(_ctx, getState());
  enterRule(_localctx, 18, FuzzyLanguageParser::RulePattern_mapping);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(123);
    match(FuzzyLanguageParser::NUMBER);
    setState(124);
    match(FuzzyLanguageParser::T__15);
    setState(125);
    configuration_pattern();
    setState(130);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__6) {
      setState(126);
      match(FuzzyLanguageParser::T__6);
      setState(127);
      configuration_pattern();
      setState(132);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Configuration_patternContext ------------------------------------------------------------------

FuzzyLanguageParser::Configuration_patternContext::Configuration_patternContext(ParserRuleContext *parent,
                                                                                size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<tree::TerminalNode *> FuzzyLanguageParser::Configuration_patternContext::IDENTIFIER() {
  return getTokens(FuzzyLanguageParser::IDENTIFIER);
}

tree::TerminalNode *FuzzyLanguageParser::Configuration_patternContext::IDENTIFIER(size_t i) {
  return getToken(FuzzyLanguageParser::IDENTIFIER, i);
}

std::vector<tree::TerminalNode *> FuzzyLanguageParser::Configuration_patternContext::STRING() {
  return getTokens(FuzzyLanguageParser::STRING);
}

tree::TerminalNode *FuzzyLanguageParser::Configuration_patternContext::STRING(size_t i) {
  return getToken(FuzzyLanguageParser::STRING, i);
}

size_t FuzzyLanguageParser::Configuration_patternContext::getRuleIndex() const {
  return FuzzyLanguageParser::RuleConfiguration_pattern;
}

std::any FuzzyLanguageParser::Configuration_patternContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitConfiguration_pattern(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Configuration_patternContext *FuzzyLanguageParser::configuration_pattern() {
  Configuration_patternContext *_localctx = _tracker.createInstance<Configuration_patternContext>(_ctx, getState());
  enterRule(_localctx, 20, FuzzyLanguageParser::RuleConfiguration_pattern);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(133);
    match(FuzzyLanguageParser::T__16);

    setState(134);
    match(FuzzyLanguageParser::IDENTIFIER);
    setState(135);
    match(FuzzyLanguageParser::T__17);
    setState(136);
    match(FuzzyLanguageParser::STRING);
    setState(144);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__6) {
      setState(138);
      match(FuzzyLanguageParser::T__6);
      setState(139);
      match(FuzzyLanguageParser::IDENTIFIER);
      setState(140);
      match(FuzzyLanguageParser::T__17);
      setState(141);
      match(FuzzyLanguageParser::STRING);
      setState(146);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(147);
    match(FuzzyLanguageParser::T__18);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

bool FuzzyLanguageParser::sempred(RuleContext *context, size_t ruleIndex, size_t predicateIndex) {
  switch (ruleIndex) {
    case 6:
      return fuzzy_setSempred(antlrcpp::downCast<Fuzzy_setContext *>(context), predicateIndex);

    default:
      break;
  }
  return true;
}

bool FuzzyLanguageParser::fuzzy_setSempred(Fuzzy_setContext *_localctx, size_t predicateIndex) {
  switch (predicateIndex) {
    case 0:
      return precpred(_ctx, 4);
    case 1:
      return precpred(_ctx, 3);

    default:
      break;
  }
  return true;
}

void FuzzyLanguageParser::initialize() {
#if ANTLR4_USE_THREAD_LOCAL_CACHE
  fuzzylanguageParserInitialize();
#else
  ::antlr4::internal::call_once(fuzzylanguageParserOnceFlag, fuzzylanguageParserInitialize);
#endif
}

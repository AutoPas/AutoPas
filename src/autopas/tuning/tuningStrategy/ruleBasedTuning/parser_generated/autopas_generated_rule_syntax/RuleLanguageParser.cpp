
// Generated from RuleLanguage.g4 by ANTLR 4.13.2

#include "RuleLanguageParser.h"

#include "RuleLanguageVisitor.h"

using namespace antlrcpp;
using namespace AutopasGeneratedRuleSyntax;

using namespace antlr4;

namespace {

struct RuleLanguageParserStaticData final {
  RuleLanguageParserStaticData(std::vector<std::string> ruleNames, std::vector<std::string> literalNames,
                               std::vector<std::string> symbolicNames)
      : ruleNames(std::move(ruleNames)),
        literalNames(std::move(literalNames)),
        symbolicNames(std::move(symbolicNames)),
        vocabulary(this->literalNames, this->symbolicNames) {}

  RuleLanguageParserStaticData(const RuleLanguageParserStaticData &) = delete;
  RuleLanguageParserStaticData(RuleLanguageParserStaticData &&) = delete;
  RuleLanguageParserStaticData &operator=(const RuleLanguageParserStaticData &) = delete;
  RuleLanguageParserStaticData &operator=(RuleLanguageParserStaticData &&) = delete;

  std::vector<antlr4::dfa::DFA> decisionToDFA;
  antlr4::atn::PredictionContextCache sharedContextCache;
  const std::vector<std::string> ruleNames;
  const std::vector<std::string> literalNames;
  const std::vector<std::string> symbolicNames;
  const antlr4::dfa::Vocabulary vocabulary;
  antlr4::atn::SerializedATNView serializedATN;
  std::unique_ptr<antlr4::atn::ATN> atn;
};

::antlr4::internal::OnceFlag rulelanguageParserOnceFlag;
#if ANTLR4_USE_THREAD_LOCAL_CACHE
static thread_local
#endif
    std::unique_ptr<RuleLanguageParserStaticData>
        rulelanguageParserStaticData = nullptr;

void rulelanguageParserInitialize() {
#if ANTLR4_USE_THREAD_LOCAL_CACHE
  if (rulelanguageParserStaticData != nullptr) {
    return;
  }
#else
  assert(rulelanguageParserStaticData == nullptr);
#endif
  auto staticData = std::make_unique<RuleLanguageParserStaticData>(
      std::vector<std::string>{"program", "unsigned_val", "literal", "define_list", "variable", "expression", "define",
                               "property_value", "configuration_pattern", "configuration_order", "statement",
                               "if_statement"},
      std::vector<std::string>{"",       "'\"'",     "'define_list'", "'='", "','",   "';'",         "'*'",  "'/'",
                               "'+'",    "'-'",      "'>'",           "'<'", "'not'", "'and'",       "'or'", "'('",
                               "')'",    "'define'", "'['",           "']'", "'>='",  "'with same'", "'if'", "':'",
                               "'endif'"},
      std::vector<std::string>{"",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "",
                               "COMMENT",
                               "WS",
                               "Container_opt",
                               "Traversal_opt",
                               "Load_estimator_opt",
                               "Data_layout_opt",
                               "Newton3_opt",
                               "Bool_val",
                               "Configuration_property",
                               "DIGIT",
                               "Unsigned_val",
                               "Double_val",
                               "Variable_name"});
  static const int32_t serializedATNSegment[] = {
      4,  1,  37, 155, 2,   0,   7,   0,   2,   1,   7,   1,   2,   2,  7,  2,  2,   3,   7,   3,  2,  4,  7,   4,
      2,  5,  7,  5,   2,   6,   7,   6,   2,   7,   7,   7,   2,   8,  7,  8,  2,   9,   7,   9,  2,  10, 7,   10,
      2,  11, 7,  11,  1,   0,   4,   0,   26,  8,   0,   11,  0,   12, 0,  27, 1,   1,   1,   1,  1,  2,  1,   2,
      1,  2,  1,  2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,  1,  2,  1,   2,   1,   2,  1,  2,  1,   2,
      1,  2,  1,  2,   1,   2,   1,   2,   1,   2,   1,   2,   3,   2,  52, 8,  2,   1,   3,   1,  3,  1,  3,   1,
      3,  1,  3,  1,   3,   5,   3,   60,  8,   3,   10,  3,   12,  3,  63, 9,  3,   1,   3,   1,  3,  1,  4,   1,
      4,  1,  5,  1,   5,   1,   5,   1,   5,   1,   5,   1,   5,   1,  5,  1,  5,   1,   5,   3,  5,  78, 8,   5,
      1,  5,  1,  5,   1,   5,   1,   5,   1,   5,   1,   5,   1,   5,  1,  5,  1,   5,   1,   5,  1,  5,  1,   5,
      5,  5,  92, 8,   5,   10,  5,   12,  5,   95,  9,   5,   1,   6,  1,  6,  1,   6,   1,   6,  1,  6,  1,   6,
      1,  7,  1,  7,   3,   7,   105, 8,   7,   1,   8,   1,   8,   1,  8,  1,  8,   1,   8,   1,  8,  1,  8,   1,
      8,  1,  8,  5,   8,   116, 8,   8,   10,  8,   12,  8,   119, 9,  8,  1,  8,   1,   8,   1,  9,  1,  9,   1,
      9,  1,  9,  1,   9,   1,   9,   1,   9,   5,   9,   130, 8,   9,  10, 9,  12,  9,   133, 9,  9,  3,  9,   135,
      8,  9,  1,  9,   1,   9,   1,   10,  1,   10,  1,   10,  1,   10, 3,  10, 143, 8,   10,  1,  11, 1,  11,  1,
      11, 1,  11, 4,   11,  149, 8,   11,  11,  11,  12,  11,  150, 1,  11, 1,  11,  1,   11,  0,  1,  10, 12,  0,
      2,  4,  6,  8,   10,  12,  14,  16,  18,  20,  22,  0,   5,   1,  0,  34, 35,  1,   0,   6,  7,  1,  0,   8,
      9,  1,  0,  10,  11,  1,   0,   13,  14,  166, 0,   25,  1,   0,  0,  0,  2,   29,  1,   0,  0,  0,  4,   51,
      1,  0,  0,  0,   6,   53,  1,   0,   0,   0,   8,   66,  1,   0,  0,  0,  10,  77,  1,   0,  0,  0,  12,  96,
      1,  0,  0,  0,   14,  104, 1,   0,   0,   0,   16,  106, 1,   0,  0,  0,  18,  122, 1,   0,  0,  0,  20,  142,
      1,  0,  0,  0,   22,  144, 1,   0,   0,   0,   24,  26,  3,   20, 10, 0,  25,  24,  1,   0,  0,  0,  26,  27,
      1,  0,  0,  0,   27,  25,  1,   0,   0,   0,   27,  28,  1,   0,  0,  0,  28,  1,   1,   0,  0,  0,  29,  30,
      7,  0,  0,  0,   30,  3,   1,   0,   0,   0,   31,  32,  5,   1,  0,  0,  32,  33,  5,   28, 0,  0,  33,  52,
      5,  1,  0,  0,   34,  35,  5,   1,   0,   0,   35,  36,  5,   27, 0,  0,  36,  52,  5,   1,  0,  0,  37,  38,
      5,  1,  0,  0,   38,  39,  5,   29,  0,   0,   39,  52,  5,   1,  0,  0,  40,  41,  5,   1,  0,  0,  41,  42,
      5,  30, 0,  0,   42,  52,  5,   1,   0,   0,   43,  44,  5,   1,  0,  0,  44,  45,  5,   31, 0,  0,  45,  52,
      5,  1,  0,  0,   46,  52,  3,   2,   1,   0,   47,  52,  5,   36, 0,  0,  48,  49,  5,   1,  0,  0,  49,  50,
      5,  32, 0,  0,   50,  52,  5,   1,   0,   0,   51,  31,  1,   0,  0,  0,  51,  34,  1,   0,  0,  0,  51,  37,
      1,  0,  0,  0,   51,  40,  1,   0,   0,   0,   51,  43,  1,   0,  0,  0,  51,  46,  1,   0,  0,  0,  51,  47,
      1,  0,  0,  0,   51,  48,  1,   0,   0,   0,   52,  5,   1,   0,  0,  0,  53,  54,  5,   2,  0,  0,  54,  55,
      5,  37, 0,  0,   55,  56,  5,   3,   0,   0,   56,  61,  3,   4,  2,  0,  57,  58,  5,   4,  0,  0,  58,  60,
      3,  4,  2,  0,   59,  57,  1,   0,   0,   0,   60,  63,  1,   0,  0,  0,  61,  59,  1,   0,  0,  0,  61,  62,
      1,  0,  0,  0,   62,  64,  1,   0,   0,   0,   63,  61,  1,   0,  0,  0,  64,  65,  5,   5,  0,  0,  65,  7,
      1,  0,  0,  0,   66,  67,  5,   37,  0,   0,   67,  9,   1,   0,  0,  0,  68,  69,  6,   5,  -1, 0,  69,  70,
      5,  12, 0,  0,   70,  78,  3,   10,  5,   5,   71,  78,  3,   4,  2,  0,  72,  78,  3,   8,  4,  0,  73,  74,
      5,  15, 0,  0,   74,  75,  3,   10,  5,   0,   75,  76,  5,   16, 0,  0,  76,  78,  1,   0,  0,  0,  77,  68,
      1,  0,  0,  0,   77,  71,  1,   0,   0,   0,   77,  72,  1,   0,  0,  0,  77,  73,  1,   0,  0,  0,  78,  93,
      1,  0,  0,  0,   79,  80,  10,  8,   0,   0,   80,  81,  7,   1,  0,  0,  81,  92,  3,   10, 5,  9,  82,  83,
      10, 7,  0,  0,   83,  84,  7,   2,   0,   0,   84,  92,  3,   10, 5,  8,  85,  86,  10,  6,  0,  0,  86,  87,
      7,  3,  0,  0,   87,  92,  3,   10,  5,   7,   88,  89,  10,  4,  0,  0,  89,  90,  7,   4,  0,  0,  90,  92,
      3,  10, 5,  5,   91,  79,  1,   0,   0,   0,   91,  82,  1,   0,  0,  0,  91,  85,  1,   0,  0,  0,  91,  88,
      1,  0,  0,  0,   92,  95,  1,   0,   0,   0,   93,  91,  1,   0,  0,  0,  93,  94,  1,   0,  0,  0,  94,  11,
      1,  0,  0,  0,   95,  93,  1,   0,   0,   0,   96,  97,  5,   17, 0,  0,  97,  98,  5,   37, 0,  0,  98,  99,
      5,  3,  0,  0,   99,  100, 3,   10,  5,   0,   100, 101, 5,   5,  0,  0,  101, 13,  1,   0,  0,  0,  102, 105,
      5,  37, 0,  0,   103, 105, 3,   4,   2,   0,   104, 102, 1,   0,  0,  0,  104, 103, 1,   0,  0,  0,  105, 15,
      1,  0,  0,  0,   106, 107, 5,   18,  0,   0,   107, 108, 5,   33, 0,  0,  108, 109, 5,   3,  0,  0,  109, 110,
      3,  14, 7,  0,   110, 117, 1,   0,   0,   0,   111, 112, 5,   4,  0,  0,  112, 113, 5,   33, 0,  0,  113, 114,
      5,  3,  0,  0,   114, 116, 3,   14,  7,   0,   115, 111, 1,   0,  0,  0,  116, 119, 1,   0,  0,  0,  117, 115,
      1,  0,  0,  0,   117, 118, 1,   0,   0,   0,   118, 120, 1,   0,  0,  0,  119, 117, 1,   0,  0,  0,  120, 121,
      5,  19, 0,  0,   121, 17,  1,   0,   0,   0,   122, 123, 3,   16, 8,  0,  123, 124, 5,   20, 0,  0,  124, 134,
      3,  16, 8,  0,   125, 126, 5,   21,  0,   0,   126, 131, 5,   33, 0,  0,  127, 128, 5,   4,  0,  0,  128, 130,
      5,  33, 0,  0,   129, 127, 1,   0,   0,   0,   130, 133, 1,   0,  0,  0,  131, 129, 1,   0,  0,  0,  131, 132,
      1,  0,  0,  0,   132, 135, 1,   0,   0,   0,   133, 131, 1,   0,  0,  0,  134, 125, 1,   0,  0,  0,  134, 135,
      1,  0,  0,  0,   135, 136, 1,   0,   0,   0,   136, 137, 5,   5,  0,  0,  137, 19,  1,   0,  0,  0,  138, 143,
      3,  6,  3,  0,   139, 143, 3,   12,  6,   0,   140, 143, 3,   22, 11, 0,  141, 143, 3,   18, 9,  0,  142, 138,
      1,  0,  0,  0,   142, 139, 1,   0,   0,   0,   142, 140, 1,   0,  0,  0,  142, 141, 1,   0,  0,  0,  143, 21,
      1,  0,  0,  0,   144, 145, 5,   22,  0,   0,   145, 146, 3,   10, 5,  0,  146, 148, 5,   23, 0,  0,  147, 149,
      3,  20, 10, 0,   148, 147, 1,   0,   0,   0,   149, 150, 1,   0,  0,  0,  150, 148, 1,   0,  0,  0,  150, 151,
      1,  0,  0,  0,   151, 152, 1,   0,   0,   0,   152, 153, 5,   24, 0,  0,  153, 23,  1,   0,  0,  0,  12,  27,
      51, 61, 77, 91,  93,  104, 117, 131, 134, 142, 150};
  staticData->serializedATN = antlr4::atn::SerializedATNView(
      serializedATNSegment, sizeof(serializedATNSegment) / sizeof(serializedATNSegment[0]));

  antlr4::atn::ATNDeserializer deserializer;
  staticData->atn = deserializer.deserialize(staticData->serializedATN);

  const size_t count = staticData->atn->getNumberOfDecisions();
  staticData->decisionToDFA.reserve(count);
  for (size_t i = 0; i < count; i++) {
    staticData->decisionToDFA.emplace_back(staticData->atn->getDecisionState(i), i);
  }
  rulelanguageParserStaticData = std::move(staticData);
}

}  // namespace

RuleLanguageParser::RuleLanguageParser(TokenStream *input)
    : RuleLanguageParser(input, antlr4::atn::ParserATNSimulatorOptions()) {}

RuleLanguageParser::RuleLanguageParser(TokenStream *input, const antlr4::atn::ParserATNSimulatorOptions &options)
    : Parser(input) {
  RuleLanguageParser::initialize();
  _interpreter =
      new atn::ParserATNSimulator(this, *rulelanguageParserStaticData->atn, rulelanguageParserStaticData->decisionToDFA,
                                  rulelanguageParserStaticData->sharedContextCache, options);
}

RuleLanguageParser::~RuleLanguageParser() { delete _interpreter; }

const atn::ATN &RuleLanguageParser::getATN() const { return *rulelanguageParserStaticData->atn; }

std::string RuleLanguageParser::getGrammarFileName() const { return "RuleLanguage.g4"; }

const std::vector<std::string> &RuleLanguageParser::getRuleNames() const {
  return rulelanguageParserStaticData->ruleNames;
}

const dfa::Vocabulary &RuleLanguageParser::getVocabulary() const { return rulelanguageParserStaticData->vocabulary; }

antlr4::atn::SerializedATNView RuleLanguageParser::getSerializedATN() const {
  return rulelanguageParserStaticData->serializedATN;
}

//----------------- ProgramContext ------------------------------------------------------------------

RuleLanguageParser::ProgramContext::ProgramContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<RuleLanguageParser::StatementContext *> RuleLanguageParser::ProgramContext::statement() {
  return getRuleContexts<RuleLanguageParser::StatementContext>();
}

RuleLanguageParser::StatementContext *RuleLanguageParser::ProgramContext::statement(size_t i) {
  return getRuleContext<RuleLanguageParser::StatementContext>(i);
}

size_t RuleLanguageParser::ProgramContext::getRuleIndex() const { return RuleLanguageParser::RuleProgram; }

std::any RuleLanguageParser::ProgramContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitProgram(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::ProgramContext *RuleLanguageParser::program() {
  ProgramContext *_localctx = _tracker.createInstance<ProgramContext>(_ctx, getState());
  enterRule(_localctx, 0, RuleLanguageParser::RuleProgram);
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
    setState(25);
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(24);
      statement();
      setState(27);
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while ((((_la & ~0x3fULL) == 0) && ((1ULL << _la) & 4587524) != 0));

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Unsigned_valContext ------------------------------------------------------------------

RuleLanguageParser::Unsigned_valContext::Unsigned_valContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *RuleLanguageParser::Unsigned_valContext::Unsigned_val() {
  return getToken(RuleLanguageParser::Unsigned_val, 0);
}

tree::TerminalNode *RuleLanguageParser::Unsigned_valContext::DIGIT() { return getToken(RuleLanguageParser::DIGIT, 0); }

size_t RuleLanguageParser::Unsigned_valContext::getRuleIndex() const { return RuleLanguageParser::RuleUnsigned_val; }

std::any RuleLanguageParser::Unsigned_valContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitUnsigned_val(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Unsigned_valContext *RuleLanguageParser::unsigned_val() {
  Unsigned_valContext *_localctx = _tracker.createInstance<Unsigned_valContext>(_ctx, getState());
  enterRule(_localctx, 2, RuleLanguageParser::RuleUnsigned_val);
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
    setState(29);
    _la = _input->LA(1);
    if (!(_la == RuleLanguageParser::DIGIT

          || _la == RuleLanguageParser::Unsigned_val)) {
      _errHandler->recoverInline(this);
    } else {
      _errHandler->reportMatch(this);
      consume();
    }

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- LiteralContext ------------------------------------------------------------------

RuleLanguageParser::LiteralContext::LiteralContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *RuleLanguageParser::LiteralContext::Traversal_opt() {
  return getToken(RuleLanguageParser::Traversal_opt, 0);
}

tree::TerminalNode *RuleLanguageParser::LiteralContext::Container_opt() {
  return getToken(RuleLanguageParser::Container_opt, 0);
}

tree::TerminalNode *RuleLanguageParser::LiteralContext::Load_estimator_opt() {
  return getToken(RuleLanguageParser::Load_estimator_opt, 0);
}

tree::TerminalNode *RuleLanguageParser::LiteralContext::Data_layout_opt() {
  return getToken(RuleLanguageParser::Data_layout_opt, 0);
}

tree::TerminalNode *RuleLanguageParser::LiteralContext::Newton3_opt() {
  return getToken(RuleLanguageParser::Newton3_opt, 0);
}

RuleLanguageParser::Unsigned_valContext *RuleLanguageParser::LiteralContext::unsigned_val() {
  return getRuleContext<RuleLanguageParser::Unsigned_valContext>(0);
}

tree::TerminalNode *RuleLanguageParser::LiteralContext::Double_val() {
  return getToken(RuleLanguageParser::Double_val, 0);
}

tree::TerminalNode *RuleLanguageParser::LiteralContext::Bool_val() { return getToken(RuleLanguageParser::Bool_val, 0); }

size_t RuleLanguageParser::LiteralContext::getRuleIndex() const { return RuleLanguageParser::RuleLiteral; }

std::any RuleLanguageParser::LiteralContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitLiteral(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::LiteralContext *RuleLanguageParser::literal() {
  LiteralContext *_localctx = _tracker.createInstance<LiteralContext>(_ctx, getState());
  enterRule(_localctx, 4, RuleLanguageParser::RuleLiteral);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(51);
    _errHandler->sync(this);
    switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 1, _ctx)) {
      case 1: {
        enterOuterAlt(_localctx, 1);
        setState(31);
        match(RuleLanguageParser::T__0);
        setState(32);
        match(RuleLanguageParser::Traversal_opt);
        setState(33);
        match(RuleLanguageParser::T__0);
        break;
      }

      case 2: {
        enterOuterAlt(_localctx, 2);
        setState(34);
        match(RuleLanguageParser::T__0);
        setState(35);
        match(RuleLanguageParser::Container_opt);
        setState(36);
        match(RuleLanguageParser::T__0);
        break;
      }

      case 3: {
        enterOuterAlt(_localctx, 3);
        setState(37);
        match(RuleLanguageParser::T__0);
        setState(38);
        match(RuleLanguageParser::Load_estimator_opt);
        setState(39);
        match(RuleLanguageParser::T__0);
        break;
      }

      case 4: {
        enterOuterAlt(_localctx, 4);
        setState(40);
        match(RuleLanguageParser::T__0);
        setState(41);
        match(RuleLanguageParser::Data_layout_opt);
        setState(42);
        match(RuleLanguageParser::T__0);
        break;
      }

      case 5: {
        enterOuterAlt(_localctx, 5);
        setState(43);
        match(RuleLanguageParser::T__0);
        setState(44);
        match(RuleLanguageParser::Newton3_opt);
        setState(45);
        match(RuleLanguageParser::T__0);
        break;
      }

      case 6: {
        enterOuterAlt(_localctx, 6);
        setState(46);
        unsigned_val();
        break;
      }

      case 7: {
        enterOuterAlt(_localctx, 7);
        setState(47);
        match(RuleLanguageParser::Double_val);
        break;
      }

      case 8: {
        enterOuterAlt(_localctx, 8);
        setState(48);
        match(RuleLanguageParser::T__0);
        setState(49);
        match(RuleLanguageParser::Bool_val);
        setState(50);
        match(RuleLanguageParser::T__0);
        break;
      }

      default:
        break;
    }

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Define_listContext ------------------------------------------------------------------

RuleLanguageParser::Define_listContext::Define_listContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *RuleLanguageParser::Define_listContext::Variable_name() {
  return getToken(RuleLanguageParser::Variable_name, 0);
}

std::vector<RuleLanguageParser::LiteralContext *> RuleLanguageParser::Define_listContext::literal() {
  return getRuleContexts<RuleLanguageParser::LiteralContext>();
}

RuleLanguageParser::LiteralContext *RuleLanguageParser::Define_listContext::literal(size_t i) {
  return getRuleContext<RuleLanguageParser::LiteralContext>(i);
}

size_t RuleLanguageParser::Define_listContext::getRuleIndex() const { return RuleLanguageParser::RuleDefine_list; }

std::any RuleLanguageParser::Define_listContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitDefine_list(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Define_listContext *RuleLanguageParser::define_list() {
  Define_listContext *_localctx = _tracker.createInstance<Define_listContext>(_ctx, getState());
  enterRule(_localctx, 6, RuleLanguageParser::RuleDefine_list);
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
    setState(53);
    match(RuleLanguageParser::T__1);
    setState(54);
    match(RuleLanguageParser::Variable_name);
    setState(55);
    match(RuleLanguageParser::T__2);
    setState(56);
    literal();
    setState(61);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == RuleLanguageParser::T__3) {
      setState(57);
      match(RuleLanguageParser::T__3);
      setState(58);
      literal();
      setState(63);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(64);
    match(RuleLanguageParser::T__4);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- VariableContext ------------------------------------------------------------------

RuleLanguageParser::VariableContext::VariableContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *RuleLanguageParser::VariableContext::Variable_name() {
  return getToken(RuleLanguageParser::Variable_name, 0);
}

size_t RuleLanguageParser::VariableContext::getRuleIndex() const { return RuleLanguageParser::RuleVariable; }

std::any RuleLanguageParser::VariableContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitVariable(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::VariableContext *RuleLanguageParser::variable() {
  VariableContext *_localctx = _tracker.createInstance<VariableContext>(_ctx, getState());
  enterRule(_localctx, 8, RuleLanguageParser::RuleVariable);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(66);
    match(RuleLanguageParser::Variable_name);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ExpressionContext ------------------------------------------------------------------

RuleLanguageParser::ExpressionContext::ExpressionContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<RuleLanguageParser::ExpressionContext *> RuleLanguageParser::ExpressionContext::expression() {
  return getRuleContexts<RuleLanguageParser::ExpressionContext>();
}

RuleLanguageParser::ExpressionContext *RuleLanguageParser::ExpressionContext::expression(size_t i) {
  return getRuleContext<RuleLanguageParser::ExpressionContext>(i);
}

RuleLanguageParser::LiteralContext *RuleLanguageParser::ExpressionContext::literal() {
  return getRuleContext<RuleLanguageParser::LiteralContext>(0);
}

RuleLanguageParser::VariableContext *RuleLanguageParser::ExpressionContext::variable() {
  return getRuleContext<RuleLanguageParser::VariableContext>(0);
}

size_t RuleLanguageParser::ExpressionContext::getRuleIndex() const { return RuleLanguageParser::RuleExpression; }

std::any RuleLanguageParser::ExpressionContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitExpression(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::ExpressionContext *RuleLanguageParser::expression() { return expression(0); }

RuleLanguageParser::ExpressionContext *RuleLanguageParser::expression(int precedence) {
  ParserRuleContext *parentContext = _ctx;
  size_t parentState = getState();
  RuleLanguageParser::ExpressionContext *_localctx = _tracker.createInstance<ExpressionContext>(_ctx, parentState);
  RuleLanguageParser::ExpressionContext *previousContext = _localctx;
  (void)previousContext;  // Silence compiler, in case the context is not used by generated code.
  size_t startState = 10;
  enterRecursionRule(_localctx, 10, RuleLanguageParser::RuleExpression, precedence);

  size_t _la = 0;

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
    setState(77);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case RuleLanguageParser::T__11: {
        setState(69);
        antlrcpp::downCast<ExpressionContext *>(_localctx)->op = match(RuleLanguageParser::T__11);
        setState(70);
        expression(5);
        break;
      }

      case RuleLanguageParser::T__0:
      case RuleLanguageParser::DIGIT:
      case RuleLanguageParser::Unsigned_val:
      case RuleLanguageParser::Double_val: {
        setState(71);
        literal();
        break;
      }

      case RuleLanguageParser::Variable_name: {
        setState(72);
        variable();
        break;
      }

      case RuleLanguageParser::T__14: {
        setState(73);
        match(RuleLanguageParser::T__14);
        setState(74);
        expression(0);
        setState(75);
        match(RuleLanguageParser::T__15);
        break;
      }

      default:
        throw NoViableAltException(this);
    }
    _ctx->stop = _input->LT(-1);
    setState(93);
    _errHandler->sync(this);
    alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 5, _ctx);
    while (alt != 2 && alt != atn::ATN::INVALID_ALT_NUMBER) {
      if (alt == 1) {
        if (!_parseListeners.empty()) triggerExitRuleEvent();
        previousContext = _localctx;
        setState(91);
        _errHandler->sync(this);
        switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 4, _ctx)) {
          case 1: {
            _localctx = _tracker.createInstance<ExpressionContext>(parentContext, parentState);
            pushNewRecursionContext(_localctx, startState, RuleExpression);
            setState(79);

            if (!(precpred(_ctx, 8))) throw FailedPredicateException(this, "precpred(_ctx, 8)");
            setState(80);
            antlrcpp::downCast<ExpressionContext *>(_localctx)->op = _input->LT(1);
            _la = _input->LA(1);
            if (!(_la == RuleLanguageParser::T__5

                  || _la == RuleLanguageParser::T__6)) {
              antlrcpp::downCast<ExpressionContext *>(_localctx)->op = _errHandler->recoverInline(this);
            } else {
              _errHandler->reportMatch(this);
              consume();
            }
            setState(81);
            expression(9);
            break;
          }

          case 2: {
            _localctx = _tracker.createInstance<ExpressionContext>(parentContext, parentState);
            pushNewRecursionContext(_localctx, startState, RuleExpression);
            setState(82);

            if (!(precpred(_ctx, 7))) throw FailedPredicateException(this, "precpred(_ctx, 7)");
            setState(83);
            antlrcpp::downCast<ExpressionContext *>(_localctx)->op = _input->LT(1);
            _la = _input->LA(1);
            if (!(_la == RuleLanguageParser::T__7

                  || _la == RuleLanguageParser::T__8)) {
              antlrcpp::downCast<ExpressionContext *>(_localctx)->op = _errHandler->recoverInline(this);
            } else {
              _errHandler->reportMatch(this);
              consume();
            }
            setState(84);
            expression(8);
            break;
          }

          case 3: {
            _localctx = _tracker.createInstance<ExpressionContext>(parentContext, parentState);
            pushNewRecursionContext(_localctx, startState, RuleExpression);
            setState(85);

            if (!(precpred(_ctx, 6))) throw FailedPredicateException(this, "precpred(_ctx, 6)");
            setState(86);
            antlrcpp::downCast<ExpressionContext *>(_localctx)->op = _input->LT(1);
            _la = _input->LA(1);
            if (!(_la == RuleLanguageParser::T__9

                  || _la == RuleLanguageParser::T__10)) {
              antlrcpp::downCast<ExpressionContext *>(_localctx)->op = _errHandler->recoverInline(this);
            } else {
              _errHandler->reportMatch(this);
              consume();
            }
            setState(87);
            expression(7);
            break;
          }

          case 4: {
            _localctx = _tracker.createInstance<ExpressionContext>(parentContext, parentState);
            pushNewRecursionContext(_localctx, startState, RuleExpression);
            setState(88);

            if (!(precpred(_ctx, 4))) throw FailedPredicateException(this, "precpred(_ctx, 4)");
            setState(89);
            antlrcpp::downCast<ExpressionContext *>(_localctx)->op = _input->LT(1);
            _la = _input->LA(1);
            if (!(_la == RuleLanguageParser::T__12

                  || _la == RuleLanguageParser::T__13)) {
              antlrcpp::downCast<ExpressionContext *>(_localctx)->op = _errHandler->recoverInline(this);
            } else {
              _errHandler->reportMatch(this);
              consume();
            }
            setState(90);
            expression(5);
            break;
          }

          default:
            break;
        }
      }
      setState(95);
      _errHandler->sync(this);
      alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 5, _ctx);
    }
  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }
  return _localctx;
}

//----------------- DefineContext ------------------------------------------------------------------

RuleLanguageParser::DefineContext::DefineContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *RuleLanguageParser::DefineContext::Variable_name() {
  return getToken(RuleLanguageParser::Variable_name, 0);
}

RuleLanguageParser::ExpressionContext *RuleLanguageParser::DefineContext::expression() {
  return getRuleContext<RuleLanguageParser::ExpressionContext>(0);
}

size_t RuleLanguageParser::DefineContext::getRuleIndex() const { return RuleLanguageParser::RuleDefine; }

std::any RuleLanguageParser::DefineContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitDefine(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::DefineContext *RuleLanguageParser::define() {
  DefineContext *_localctx = _tracker.createInstance<DefineContext>(_ctx, getState());
  enterRule(_localctx, 12, RuleLanguageParser::RuleDefine);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(96);
    match(RuleLanguageParser::T__16);
    setState(97);
    match(RuleLanguageParser::Variable_name);
    setState(98);
    match(RuleLanguageParser::T__2);
    setState(99);
    expression(0);
    setState(100);
    match(RuleLanguageParser::T__4);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Property_valueContext ------------------------------------------------------------------

RuleLanguageParser::Property_valueContext::Property_valueContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *RuleLanguageParser::Property_valueContext::Variable_name() {
  return getToken(RuleLanguageParser::Variable_name, 0);
}

RuleLanguageParser::LiteralContext *RuleLanguageParser::Property_valueContext::literal() {
  return getRuleContext<RuleLanguageParser::LiteralContext>(0);
}

size_t RuleLanguageParser::Property_valueContext::getRuleIndex() const {
  return RuleLanguageParser::RuleProperty_value;
}

std::any RuleLanguageParser::Property_valueContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitProperty_value(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Property_valueContext *RuleLanguageParser::property_value() {
  Property_valueContext *_localctx = _tracker.createInstance<Property_valueContext>(_ctx, getState());
  enterRule(_localctx, 14, RuleLanguageParser::RuleProperty_value);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(104);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case RuleLanguageParser::Variable_name: {
        enterOuterAlt(_localctx, 1);
        setState(102);
        match(RuleLanguageParser::Variable_name);
        break;
      }

      case RuleLanguageParser::T__0:
      case RuleLanguageParser::DIGIT:
      case RuleLanguageParser::Unsigned_val:
      case RuleLanguageParser::Double_val: {
        enterOuterAlt(_localctx, 2);
        setState(103);
        literal();
        break;
      }

      default:
        throw NoViableAltException(this);
    }

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Configuration_patternContext ------------------------------------------------------------------

RuleLanguageParser::Configuration_patternContext::Configuration_patternContext(ParserRuleContext *parent,
                                                                               size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<tree::TerminalNode *> RuleLanguageParser::Configuration_patternContext::Configuration_property() {
  return getTokens(RuleLanguageParser::Configuration_property);
}

tree::TerminalNode *RuleLanguageParser::Configuration_patternContext::Configuration_property(size_t i) {
  return getToken(RuleLanguageParser::Configuration_property, i);
}

std::vector<RuleLanguageParser::Property_valueContext *>
RuleLanguageParser::Configuration_patternContext::property_value() {
  return getRuleContexts<RuleLanguageParser::Property_valueContext>();
}

RuleLanguageParser::Property_valueContext *RuleLanguageParser::Configuration_patternContext::property_value(size_t i) {
  return getRuleContext<RuleLanguageParser::Property_valueContext>(i);
}

size_t RuleLanguageParser::Configuration_patternContext::getRuleIndex() const {
  return RuleLanguageParser::RuleConfiguration_pattern;
}

std::any RuleLanguageParser::Configuration_patternContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitConfiguration_pattern(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Configuration_patternContext *RuleLanguageParser::configuration_pattern() {
  Configuration_patternContext *_localctx = _tracker.createInstance<Configuration_patternContext>(_ctx, getState());
  enterRule(_localctx, 16, RuleLanguageParser::RuleConfiguration_pattern);
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
    setState(106);
    match(RuleLanguageParser::T__17);

    setState(107);
    match(RuleLanguageParser::Configuration_property);
    setState(108);
    match(RuleLanguageParser::T__2);
    setState(109);
    property_value();
    setState(117);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == RuleLanguageParser::T__3) {
      setState(111);
      match(RuleLanguageParser::T__3);
      setState(112);
      match(RuleLanguageParser::Configuration_property);
      setState(113);
      match(RuleLanguageParser::T__2);
      setState(114);
      property_value();
      setState(119);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(120);
    match(RuleLanguageParser::T__18);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Configuration_orderContext ------------------------------------------------------------------

RuleLanguageParser::Configuration_orderContext::Configuration_orderContext(ParserRuleContext *parent,
                                                                           size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<RuleLanguageParser::Configuration_patternContext *>
RuleLanguageParser::Configuration_orderContext::configuration_pattern() {
  return getRuleContexts<RuleLanguageParser::Configuration_patternContext>();
}

RuleLanguageParser::Configuration_patternContext *RuleLanguageParser::Configuration_orderContext::configuration_pattern(
    size_t i) {
  return getRuleContext<RuleLanguageParser::Configuration_patternContext>(i);
}

std::vector<tree::TerminalNode *> RuleLanguageParser::Configuration_orderContext::Configuration_property() {
  return getTokens(RuleLanguageParser::Configuration_property);
}

tree::TerminalNode *RuleLanguageParser::Configuration_orderContext::Configuration_property(size_t i) {
  return getToken(RuleLanguageParser::Configuration_property, i);
}

size_t RuleLanguageParser::Configuration_orderContext::getRuleIndex() const {
  return RuleLanguageParser::RuleConfiguration_order;
}

std::any RuleLanguageParser::Configuration_orderContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitConfiguration_order(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Configuration_orderContext *RuleLanguageParser::configuration_order() {
  Configuration_orderContext *_localctx = _tracker.createInstance<Configuration_orderContext>(_ctx, getState());
  enterRule(_localctx, 18, RuleLanguageParser::RuleConfiguration_order);
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
    setState(122);
    configuration_pattern();
    setState(123);
    match(RuleLanguageParser::T__19);
    setState(124);
    configuration_pattern();
    setState(134);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == RuleLanguageParser::T__20) {
      setState(125);
      match(RuleLanguageParser::T__20);
      setState(126);
      match(RuleLanguageParser::Configuration_property);
      setState(131);
      _errHandler->sync(this);
      _la = _input->LA(1);
      while (_la == RuleLanguageParser::T__3) {
        setState(127);
        match(RuleLanguageParser::T__3);
        setState(128);
        match(RuleLanguageParser::Configuration_property);
        setState(133);
        _errHandler->sync(this);
        _la = _input->LA(1);
      }
    }
    setState(136);
    match(RuleLanguageParser::T__4);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- StatementContext ------------------------------------------------------------------

RuleLanguageParser::StatementContext::StatementContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

RuleLanguageParser::Define_listContext *RuleLanguageParser::StatementContext::define_list() {
  return getRuleContext<RuleLanguageParser::Define_listContext>(0);
}

RuleLanguageParser::DefineContext *RuleLanguageParser::StatementContext::define() {
  return getRuleContext<RuleLanguageParser::DefineContext>(0);
}

RuleLanguageParser::If_statementContext *RuleLanguageParser::StatementContext::if_statement() {
  return getRuleContext<RuleLanguageParser::If_statementContext>(0);
}

RuleLanguageParser::Configuration_orderContext *RuleLanguageParser::StatementContext::configuration_order() {
  return getRuleContext<RuleLanguageParser::Configuration_orderContext>(0);
}

size_t RuleLanguageParser::StatementContext::getRuleIndex() const { return RuleLanguageParser::RuleStatement; }

std::any RuleLanguageParser::StatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitStatement(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::StatementContext *RuleLanguageParser::statement() {
  StatementContext *_localctx = _tracker.createInstance<StatementContext>(_ctx, getState());
  enterRule(_localctx, 20, RuleLanguageParser::RuleStatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(142);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case RuleLanguageParser::T__1: {
        enterOuterAlt(_localctx, 1);
        setState(138);
        define_list();
        break;
      }

      case RuleLanguageParser::T__16: {
        enterOuterAlt(_localctx, 2);
        setState(139);
        define();
        break;
      }

      case RuleLanguageParser::T__21: {
        enterOuterAlt(_localctx, 3);
        setState(140);
        if_statement();
        break;
      }

      case RuleLanguageParser::T__17: {
        enterOuterAlt(_localctx, 4);
        setState(141);
        configuration_order();
        break;
      }

      default:
        throw NoViableAltException(this);
    }

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- If_statementContext ------------------------------------------------------------------

RuleLanguageParser::If_statementContext::If_statementContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

RuleLanguageParser::ExpressionContext *RuleLanguageParser::If_statementContext::expression() {
  return getRuleContext<RuleLanguageParser::ExpressionContext>(0);
}

std::vector<RuleLanguageParser::StatementContext *> RuleLanguageParser::If_statementContext::statement() {
  return getRuleContexts<RuleLanguageParser::StatementContext>();
}

RuleLanguageParser::StatementContext *RuleLanguageParser::If_statementContext::statement(size_t i) {
  return getRuleContext<RuleLanguageParser::StatementContext>(i);
}

size_t RuleLanguageParser::If_statementContext::getRuleIndex() const { return RuleLanguageParser::RuleIf_statement; }

std::any RuleLanguageParser::If_statementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor *>(visitor))
    return parserVisitor->visitIf_statement(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::If_statementContext *RuleLanguageParser::if_statement() {
  If_statementContext *_localctx = _tracker.createInstance<If_statementContext>(_ctx, getState());
  enterRule(_localctx, 22, RuleLanguageParser::RuleIf_statement);
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
    setState(144);
    match(RuleLanguageParser::T__21);
    setState(145);
    expression(0);
    setState(146);
    match(RuleLanguageParser::T__22);
    setState(148);
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(147);
      statement();
      setState(150);
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while ((((_la & ~0x3fULL) == 0) && ((1ULL << _la) & 4587524) != 0));
    setState(152);
    match(RuleLanguageParser::T__23);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

bool RuleLanguageParser::sempred(RuleContext *context, size_t ruleIndex, size_t predicateIndex) {
  switch (ruleIndex) {
    case 5:
      return expressionSempred(antlrcpp::downCast<ExpressionContext *>(context), predicateIndex);

    default:
      break;
  }
  return true;
}

bool RuleLanguageParser::expressionSempred(ExpressionContext *_localctx, size_t predicateIndex) {
  switch (predicateIndex) {
    case 0:
      return precpred(_ctx, 8);
    case 1:
      return precpred(_ctx, 7);
    case 2:
      return precpred(_ctx, 6);
    case 3:
      return precpred(_ctx, 4);

    default:
      break;
  }
  return true;
}

void RuleLanguageParser::initialize() {
#if ANTLR4_USE_THREAD_LOCAL_CACHE
  rulelanguageParserInitialize();
#else
  ::antlr4::internal::call_once(rulelanguageParserOnceFlag, rulelanguageParserInitialize);
#endif
}

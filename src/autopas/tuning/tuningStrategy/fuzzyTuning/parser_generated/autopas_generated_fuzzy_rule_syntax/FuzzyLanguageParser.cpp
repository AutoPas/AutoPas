
// Generated from
// AutoPas/src/autopas/tuning/tuningStrategy/fuzzyTuning\FuzzyLanguage.g4 by
// ANTLR 4.9.1

#include "FuzzyLanguageParser.h"

#include "FuzzyLanguageVisitor.h"

using namespace antlrcpp;
using namespace autopas_generated_fuzzy_rule_syntax;
using namespace antlr4;

FuzzyLanguageParser::FuzzyLanguageParser(TokenStream *input) : Parser(input) {
  _interpreter = new atn::ParserATNSimulator(this, _atn, _decisionToDFA, _sharedContextCache);
}

FuzzyLanguageParser::~FuzzyLanguageParser() { delete _interpreter; }

std::string FuzzyLanguageParser::getGrammarFileName() const { return "FuzzyLanguage.g4"; }

const std::vector<std::string> &FuzzyLanguageParser::getRuleNames() const { return _ruleNames; }

dfa::Vocabulary &FuzzyLanguageParser::getVocabulary() const { return _vocabulary; }

//----------------- Rule_fileContext ------------------------------------------------------------------

FuzzyLanguageParser::Rule_fileContext::Rule_fileContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

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

antlrcpp::Any FuzzyLanguageParser::Rule_fileContext::accept(tree::ParseTreeVisitor *visitor) {
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
    setState(15);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__0) {
      setState(12);
      linguistic_variable();
      setState(17);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(21);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__7) {
      setState(18);
      fuzzy_rule();
      setState(23);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(24);
    match(FuzzyLanguageParser::EOF);

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

antlrcpp::Any FuzzyLanguageParser::Linguistic_variableContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitLinguistic_variable(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Linguistic_variableContext *FuzzyLanguageParser::linguistic_variable() {
  Linguistic_variableContext *_localctx = _tracker.createInstance<Linguistic_variableContext>(_ctx, getState());
  enterRule(_localctx, 2, FuzzyLanguageParser::RuleLinguistic_variable);
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
    setState(26);
    match(FuzzyLanguageParser::T__0);
    setState(27);
    match(FuzzyLanguageParser::T__1);
    setState(28);
    match(FuzzyLanguageParser::T__2);
    setState(29);
    match(FuzzyLanguageParser::T__1);
    setState(30);
    match(FuzzyLanguageParser::STRING);
    setState(31);
    match(FuzzyLanguageParser::T__3);
    setState(32);
    match(FuzzyLanguageParser::T__1);
    setState(33);
    match(FuzzyLanguageParser::T__4);
    setState(34);
    match(FuzzyLanguageParser::NUMBER);
    setState(35);
    match(FuzzyLanguageParser::T__5);
    setState(36);
    match(FuzzyLanguageParser::NUMBER);
    setState(37);
    match(FuzzyLanguageParser::T__6);
    setState(39);
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(38);
      fuzzy_term();
      setState(41);
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

antlrcpp::Any FuzzyLanguageParser::Fuzzy_termContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitFuzzy_term(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Fuzzy_termContext *FuzzyLanguageParser::fuzzy_term() {
  Fuzzy_termContext *_localctx = _tracker.createInstance<Fuzzy_termContext>(_ctx, getState());
  enterRule(_localctx, 4, FuzzyLanguageParser::RuleFuzzy_term);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(43);
    match(FuzzyLanguageParser::STRING);
    setState(44);
    match(FuzzyLanguageParser::T__1);
    setState(45);
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

antlrcpp::Any FuzzyLanguageParser::FunctionContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitFunction(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::FunctionContext *FuzzyLanguageParser::function() {
  FunctionContext *_localctx = _tracker.createInstance<FunctionContext>(_ctx, getState());
  enterRule(_localctx, 6, FuzzyLanguageParser::RuleFunction);
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
    setState(47);
    match(FuzzyLanguageParser::IDENTIFIER);
    setState(48);
    match(FuzzyLanguageParser::T__4);
    setState(49);
    match(FuzzyLanguageParser::NUMBER);
    setState(54);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__5) {
      setState(50);
      match(FuzzyLanguageParser::T__5);
      setState(51);
      match(FuzzyLanguageParser::NUMBER);
      setState(56);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(57);
    match(FuzzyLanguageParser::T__6);

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

antlrcpp::Any FuzzyLanguageParser::Fuzzy_ruleContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitFuzzy_rule(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Fuzzy_ruleContext *FuzzyLanguageParser::fuzzy_rule() {
  Fuzzy_ruleContext *_localctx = _tracker.createInstance<Fuzzy_ruleContext>(_ctx, getState());
  enterRule(_localctx, 8, FuzzyLanguageParser::RuleFuzzy_rule);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(59);
    match(FuzzyLanguageParser::T__7);
    setState(60);
    fuzzy_set(0);
    setState(61);
    match(FuzzyLanguageParser::T__8);
    setState(62);
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

antlrcpp::Any FuzzyLanguageParser::OrContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::BracketsContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::AndContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::SelectContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::NegateContext::accept(tree::ParseTreeVisitor *visitor) {
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
  size_t startState = 10;
  enterRecursionRule(_localctx, 10, FuzzyLanguageParser::RuleFuzzy_set, precedence);

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
    setState(74);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case FuzzyLanguageParser::T__4: {
        _localctx = _tracker.createInstance<BracketsContext>(_localctx);
        _ctx = _localctx;
        previousContext = _localctx;

        setState(65);
        match(FuzzyLanguageParser::T__4);
        setState(66);
        fuzzy_set(0);
        setState(67);
        match(FuzzyLanguageParser::T__6);
        break;
      }

      case FuzzyLanguageParser::T__11: {
        _localctx = _tracker.createInstance<NegateContext>(_localctx);
        _ctx = _localctx;
        previousContext = _localctx;
        setState(69);
        match(FuzzyLanguageParser::T__11);
        setState(70);
        fuzzy_set(2);
        break;
      }

      case FuzzyLanguageParser::STRING: {
        _localctx = _tracker.createInstance<SelectContext>(_localctx);
        _ctx = _localctx;
        previousContext = _localctx;
        setState(71);
        match(FuzzyLanguageParser::STRING);
        setState(72);
        match(FuzzyLanguageParser::T__12);
        setState(73);
        match(FuzzyLanguageParser::STRING);
        break;
      }

      default:
        throw NoViableAltException(this);
    }
    _ctx->stop = _input->LT(-1);
    setState(84);
    _errHandler->sync(this);
    alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 6, _ctx);
    while (alt != 2 && alt != atn::ATN::INVALID_ALT_NUMBER) {
      if (alt == 1) {
        if (!_parseListeners.empty()) triggerExitRuleEvent();
        previousContext = _localctx;
        setState(82);
        _errHandler->sync(this);
        switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 5, _ctx)) {
          case 1: {
            auto newContext = _tracker.createInstance<AndContext>(
                _tracker.createInstance<Fuzzy_setContext>(parentContext, parentState));
            _localctx = newContext;
            pushNewRecursionContext(newContext, startState, RuleFuzzy_set);
            setState(76);

            if (!(precpred(_ctx, 4))) throw FailedPredicateException(this, "precpred(_ctx, 4)");
            setState(77);
            match(FuzzyLanguageParser::T__9);
            setState(78);
            fuzzy_set(5);
            break;
          }

          case 2: {
            auto newContext = _tracker.createInstance<OrContext>(
                _tracker.createInstance<Fuzzy_setContext>(parentContext, parentState));
            _localctx = newContext;
            pushNewRecursionContext(newContext, startState, RuleFuzzy_set);
            setState(79);

            if (!(precpred(_ctx, 3))) throw FailedPredicateException(this, "precpred(_ctx, 3)");
            setState(80);
            match(FuzzyLanguageParser::T__10);
            setState(81);
            fuzzy_set(4);
            break;
          }

          default:
            break;
        }
      }
      setState(86);
      _errHandler->sync(this);
      alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 6, _ctx);
    }
  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }
  return _localctx;
}

bool FuzzyLanguageParser::sempred(RuleContext *context, size_t ruleIndex, size_t predicateIndex) {
  switch (ruleIndex) {
    case 5:
      return fuzzy_setSempred(dynamic_cast<Fuzzy_setContext *>(context), predicateIndex);

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

// Static vars and initialization.
std::vector<dfa::DFA> FuzzyLanguageParser::_decisionToDFA;
atn::PredictionContextCache FuzzyLanguageParser::_sharedContextCache;

// We own the ATN which in turn owns the ATN states.
atn::ATN FuzzyLanguageParser::_atn;
std::vector<uint16_t> FuzzyLanguageParser::_serializedATN;

std::vector<std::string> FuzzyLanguageParser::_ruleNames = {"rule_file", "linguistic_variable", "fuzzy_term",
                                                            "function",  "fuzzy_rule",          "fuzzy_set"};

std::vector<std::string> FuzzyLanguageParser::_literalNames = {
    "",     "'FuzzyVariable'", "':'",  "'domain'", "'range'", "'('", "','", "')'",
    "'if'", "'then'",          "'&&'", "'||'",     "'!'",     "'=='"};

std::vector<std::string> FuzzyLanguageParser::_symbolicNames = {
    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "WS", "COMMENT", "STRING", "NUMBER", "IDENTIFIER"};

dfa::Vocabulary FuzzyLanguageParser::_vocabulary(_literalNames, _symbolicNames);

std::vector<std::string> FuzzyLanguageParser::_tokenNames;

FuzzyLanguageParser::Initializer::Initializer() {
  for (size_t i = 0; i < _symbolicNames.size(); ++i) {
    std::string name = _vocabulary.getLiteralName(i);
    if (name.empty()) {
      name = _vocabulary.getSymbolicName(i);
    }

    if (name.empty()) {
      _tokenNames.push_back("<INVALID>");
    } else {
      _tokenNames.push_back(name);
    }
  }

  _serializedATN = {
      0x3,  0x608b, 0xa72a, 0x8133, 0xb9ed, 0x417c, 0x3be7, 0x7786, 0x5964, 0x3,  0x14, 0x5a, 0x4,  0x2,  0x9,  0x2,
      0x4,  0x3,    0x9,    0x3,    0x4,    0x4,    0x9,    0x4,    0x4,    0x5,  0x9,  0x5,  0x4,  0x6,  0x9,  0x6,
      0x4,  0x7,    0x9,    0x7,    0x3,    0x2,    0x7,    0x2,    0x10,   0xa,  0x2,  0xc,  0x2,  0xe,  0x2,  0x13,
      0xb,  0x2,    0x3,    0x2,    0x7,    0x2,    0x16,   0xa,    0x2,    0xc,  0x2,  0xe,  0x2,  0x19, 0xb,  0x2,
      0x3,  0x2,    0x3,    0x2,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,  0x3,  0x3,  0x3,  0x3,  0x3,  0x3,
      0x3,  0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,  0x3,  0x3,  0x3,  0x3,  0x6,  0x3,
      0x2a, 0xa,    0x3,    0xd,    0x3,    0xe,    0x3,    0x2b,   0x3,    0x4,  0x3,  0x4,  0x3,  0x4,  0x3,  0x4,
      0x3,  0x5,    0x3,    0x5,    0x3,    0x5,    0x3,    0x5,    0x3,    0x5,  0x7,  0x5,  0x37, 0xa,  0x5,  0xc,
      0x5,  0xe,    0x5,    0x3a,   0xb,    0x5,    0x3,    0x5,    0x3,    0x5,  0x3,  0x6,  0x3,  0x6,  0x3,  0x6,
      0x3,  0x6,    0x3,    0x6,    0x3,    0x7,    0x3,    0x7,    0x3,    0x7,  0x3,  0x7,  0x3,  0x7,  0x3,  0x7,
      0x3,  0x7,    0x3,    0x7,    0x3,    0x7,    0x3,    0x7,    0x5,    0x7,  0x4d, 0xa,  0x7,  0x3,  0x7,  0x3,
      0x7,  0x3,    0x7,    0x3,    0x7,    0x3,    0x7,    0x3,    0x7,    0x7,  0x7,  0x55, 0xa,  0x7,  0xc,  0x7,
      0xe,  0x7,    0x58,   0xb,    0x7,    0x3,    0x7,    0x2,    0x3,    0xc,  0x8,  0x2,  0x4,  0x6,  0x8,  0xa,
      0xc,  0x2,    0x2,    0x2,    0x5b,   0x2,    0x11,   0x3,    0x2,    0x2,  0x2,  0x4,  0x1c, 0x3,  0x2,  0x2,
      0x2,  0x6,    0x2d,   0x3,    0x2,    0x2,    0x2,    0x8,    0x31,   0x3,  0x2,  0x2,  0x2,  0xa,  0x3d, 0x3,
      0x2,  0x2,    0x2,    0xc,    0x4c,   0x3,    0x2,    0x2,    0x2,    0xe,  0x10, 0x5,  0x4,  0x3,  0x2,  0xf,
      0xe,  0x3,    0x2,    0x2,    0x2,    0x10,   0x13,   0x3,    0x2,    0x2,  0x2,  0x11, 0xf,  0x3,  0x2,  0x2,
      0x2,  0x11,   0x12,   0x3,    0x2,    0x2,    0x2,    0x12,   0x17,   0x3,  0x2,  0x2,  0x2,  0x13, 0x11, 0x3,
      0x2,  0x2,    0x2,    0x14,   0x16,   0x5,    0xa,    0x6,    0x2,    0x15, 0x14, 0x3,  0x2,  0x2,  0x2,  0x16,
      0x19, 0x3,    0x2,    0x2,    0x2,    0x17,   0x15,   0x3,    0x2,    0x2,  0x2,  0x17, 0x18, 0x3,  0x2,  0x2,
      0x2,  0x18,   0x1a,   0x3,    0x2,    0x2,    0x2,    0x19,   0x17,   0x3,  0x2,  0x2,  0x2,  0x1a, 0x1b, 0x7,
      0x2,  0x2,    0x3,    0x1b,   0x3,    0x3,    0x2,    0x2,    0x2,    0x1c, 0x1d, 0x7,  0x3,  0x2,  0x2,  0x1d,
      0x1e, 0x7,    0x4,    0x2,    0x2,    0x1e,   0x1f,   0x7,    0x5,    0x2,  0x2,  0x1f, 0x20, 0x7,  0x4,  0x2,
      0x2,  0x20,   0x21,   0x7,    0x12,   0x2,    0x2,    0x21,   0x22,   0x7,  0x6,  0x2,  0x2,  0x22, 0x23, 0x7,
      0x4,  0x2,    0x2,    0x23,   0x24,   0x7,    0x7,    0x2,    0x2,    0x24, 0x25, 0x7,  0x13, 0x2,  0x2,  0x25,
      0x26, 0x7,    0x8,    0x2,    0x2,    0x26,   0x27,   0x7,    0x13,   0x2,  0x2,  0x27, 0x29, 0x7,  0x9,  0x2,
      0x2,  0x28,   0x2a,   0x5,    0x6,    0x4,    0x2,    0x29,   0x28,   0x3,  0x2,  0x2,  0x2,  0x2a, 0x2b, 0x3,
      0x2,  0x2,    0x2,    0x2b,   0x29,   0x3,    0x2,    0x2,    0x2,    0x2b, 0x2c, 0x3,  0x2,  0x2,  0x2,  0x2c,
      0x5,  0x3,    0x2,    0x2,    0x2,    0x2d,   0x2e,   0x7,    0x12,   0x2,  0x2,  0x2e, 0x2f, 0x7,  0x4,  0x2,
      0x2,  0x2f,   0x30,   0x5,    0x8,    0x5,    0x2,    0x30,   0x7,    0x3,  0x2,  0x2,  0x2,  0x31, 0x32, 0x7,
      0x14, 0x2,    0x2,    0x32,   0x33,   0x7,    0x7,    0x2,    0x2,    0x33, 0x38, 0x7,  0x13, 0x2,  0x2,  0x34,
      0x35, 0x7,    0x8,    0x2,    0x2,    0x35,   0x37,   0x7,    0x13,   0x2,  0x2,  0x36, 0x34, 0x3,  0x2,  0x2,
      0x2,  0x37,   0x3a,   0x3,    0x2,    0x2,    0x2,    0x38,   0x36,   0x3,  0x2,  0x2,  0x2,  0x38, 0x39, 0x3,
      0x2,  0x2,    0x2,    0x39,   0x3b,   0x3,    0x2,    0x2,    0x2,    0x3a, 0x38, 0x3,  0x2,  0x2,  0x2,  0x3b,
      0x3c, 0x7,    0x9,    0x2,    0x2,    0x3c,   0x9,    0x3,    0x2,    0x2,  0x2,  0x3d, 0x3e, 0x7,  0xa,  0x2,
      0x2,  0x3e,   0x3f,   0x5,    0xc,    0x7,    0x2,    0x3f,   0x40,   0x7,  0xb,  0x2,  0x2,  0x40, 0x41, 0x5,
      0xc,  0x7,    0x2,    0x41,   0xb,    0x3,    0x2,    0x2,    0x2,    0x42, 0x43, 0x8,  0x7,  0x1,  0x2,  0x43,
      0x44, 0x7,    0x7,    0x2,    0x2,    0x44,   0x45,   0x5,    0xc,    0x7,  0x2,  0x45, 0x46, 0x7,  0x9,  0x2,
      0x2,  0x46,   0x4d,   0x3,    0x2,    0x2,    0x2,    0x47,   0x48,   0x7,  0xe,  0x2,  0x2,  0x48, 0x4d, 0x5,
      0xc,  0x7,    0x4,    0x49,   0x4a,   0x7,    0x12,   0x2,    0x2,    0x4a, 0x4b, 0x7,  0xf,  0x2,  0x2,  0x4b,
      0x4d, 0x7,    0x12,   0x2,    0x2,    0x4c,   0x42,   0x3,    0x2,    0x2,  0x2,  0x4c, 0x47, 0x3,  0x2,  0x2,
      0x2,  0x4c,   0x49,   0x3,    0x2,    0x2,    0x2,    0x4d,   0x56,   0x3,  0x2,  0x2,  0x2,  0x4e, 0x4f, 0xc,
      0x6,  0x2,    0x2,    0x4f,   0x50,   0x7,    0xc,    0x2,    0x2,    0x50, 0x55, 0x5,  0xc,  0x7,  0x7,  0x51,
      0x52, 0xc,    0x5,    0x2,    0x2,    0x52,   0x53,   0x7,    0xd,    0x2,  0x2,  0x53, 0x55, 0x5,  0xc,  0x7,
      0x6,  0x54,   0x4e,   0x3,    0x2,    0x2,    0x2,    0x54,   0x51,   0x3,  0x2,  0x2,  0x2,  0x55, 0x58, 0x3,
      0x2,  0x2,    0x2,    0x56,   0x54,   0x3,    0x2,    0x2,    0x2,    0x56, 0x57, 0x3,  0x2,  0x2,  0x2,  0x57,
      0xd,  0x3,    0x2,    0x2,    0x2,    0x58,   0x56,   0x3,    0x2,    0x2,  0x2,  0x9,  0x11, 0x17, 0x2b, 0x38,
      0x4c, 0x54,   0x56,
  };

  atn::ATNDeserializer deserializer;
  _atn = deserializer.deserialize(_serializedATN);

  size_t count = _atn.getNumberOfDecisions();
  _decisionToDFA.reserve(count);
  for (size_t i = 0; i < count; i++) {
    _decisionToDFA.emplace_back(_atn.getDecisionState(i), i);
  }
}

FuzzyLanguageParser::Initializer FuzzyLanguageParser::_init;

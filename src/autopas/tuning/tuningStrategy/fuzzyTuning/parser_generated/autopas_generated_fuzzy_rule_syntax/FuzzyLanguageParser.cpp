
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

std::vector<FuzzyLanguageParser::Fuzzy_variableContext *> FuzzyLanguageParser::Rule_fileContext::fuzzy_variable() {
  return getRuleContexts<FuzzyLanguageParser::Fuzzy_variableContext>();
}

FuzzyLanguageParser::Fuzzy_variableContext *FuzzyLanguageParser::Rule_fileContext::fuzzy_variable(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_variableContext>(i);
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
    setState(19);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__0) {
      setState(16);
      fuzzy_variable();
      setState(21);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(25);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__7) {
      setState(22);
      fuzzy_rule();
      setState(27);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(28);
    match(FuzzyLanguageParser::EOF);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Fuzzy_variableContext ------------------------------------------------------------------

FuzzyLanguageParser::Fuzzy_variableContext::Fuzzy_variableContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

FuzzyLanguageParser::NameContext *FuzzyLanguageParser::Fuzzy_variableContext::name() {
  return getRuleContext<FuzzyLanguageParser::NameContext>(0);
}

std::vector<tree::TerminalNode *> FuzzyLanguageParser::Fuzzy_variableContext::NUMBER() {
  return getTokens(FuzzyLanguageParser::NUMBER);
}

tree::TerminalNode *FuzzyLanguageParser::Fuzzy_variableContext::NUMBER(size_t i) {
  return getToken(FuzzyLanguageParser::NUMBER, i);
}

std::vector<FuzzyLanguageParser::Membership_functionContext *>
FuzzyLanguageParser::Fuzzy_variableContext::membership_function() {
  return getRuleContexts<FuzzyLanguageParser::Membership_functionContext>();
}

FuzzyLanguageParser::Membership_functionContext *FuzzyLanguageParser::Fuzzy_variableContext::membership_function(
    size_t i) {
  return getRuleContext<FuzzyLanguageParser::Membership_functionContext>(i);
}

size_t FuzzyLanguageParser::Fuzzy_variableContext::getRuleIndex() const {
  return FuzzyLanguageParser::RuleFuzzy_variable;
}

antlrcpp::Any FuzzyLanguageParser::Fuzzy_variableContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitFuzzy_variable(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Fuzzy_variableContext *FuzzyLanguageParser::fuzzy_variable() {
  Fuzzy_variableContext *_localctx = _tracker.createInstance<Fuzzy_variableContext>(_ctx, getState());
  enterRule(_localctx, 2, FuzzyLanguageParser::RuleFuzzy_variable);
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
    setState(30);
    match(FuzzyLanguageParser::T__0);
    setState(31);
    match(FuzzyLanguageParser::T__1);
    setState(32);
    match(FuzzyLanguageParser::T__2);
    setState(33);
    match(FuzzyLanguageParser::T__1);
    setState(34);
    name();
    setState(35);
    match(FuzzyLanguageParser::T__3);
    setState(36);
    match(FuzzyLanguageParser::T__1);
    setState(37);
    match(FuzzyLanguageParser::T__4);
    setState(38);
    match(FuzzyLanguageParser::NUMBER);
    setState(39);
    match(FuzzyLanguageParser::T__5);
    setState(40);
    match(FuzzyLanguageParser::NUMBER);
    setState(41);
    match(FuzzyLanguageParser::T__6);
    setState(43);
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(42);
      membership_function();
      setState(45);
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

//----------------- Membership_functionContext ------------------------------------------------------------------

FuzzyLanguageParser::Membership_functionContext::Membership_functionContext(ParserRuleContext *parent,
                                                                            size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

FuzzyLanguageParser::NameContext *FuzzyLanguageParser::Membership_functionContext::name() {
  return getRuleContext<FuzzyLanguageParser::NameContext>(0);
}

FuzzyLanguageParser::FunctionContext *FuzzyLanguageParser::Membership_functionContext::function() {
  return getRuleContext<FuzzyLanguageParser::FunctionContext>(0);
}

size_t FuzzyLanguageParser::Membership_functionContext::getRuleIndex() const {
  return FuzzyLanguageParser::RuleMembership_function;
}

antlrcpp::Any FuzzyLanguageParser::Membership_functionContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitMembership_function(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::Membership_functionContext *FuzzyLanguageParser::membership_function() {
  Membership_functionContext *_localctx = _tracker.createInstance<Membership_functionContext>(_ctx, getState());
  enterRule(_localctx, 4, FuzzyLanguageParser::RuleMembership_function);

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
    name();
    setState(48);
    match(FuzzyLanguageParser::T__1);
    setState(49);
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
    setState(51);
    match(FuzzyLanguageParser::IDENTIFIER);
    setState(52);
    match(FuzzyLanguageParser::T__4);
    setState(53);
    match(FuzzyLanguageParser::NUMBER);
    setState(58);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == FuzzyLanguageParser::T__5) {
      setState(54);
      match(FuzzyLanguageParser::T__5);
      setState(55);
      match(FuzzyLanguageParser::NUMBER);
      setState(60);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(61);
    match(FuzzyLanguageParser::T__6);

  } catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- NameContext ------------------------------------------------------------------

FuzzyLanguageParser::NameContext::NameContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

tree::TerminalNode *FuzzyLanguageParser::NameContext::STRING() { return getToken(FuzzyLanguageParser::STRING, 0); }

size_t FuzzyLanguageParser::NameContext::getRuleIndex() const { return FuzzyLanguageParser::RuleName; }

antlrcpp::Any FuzzyLanguageParser::NameContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitName(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::NameContext *FuzzyLanguageParser::name() {
  NameContext *_localctx = _tracker.createInstance<NameContext>(_ctx, getState());
  enterRule(_localctx, 8, FuzzyLanguageParser::RuleName);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(63);
    match(FuzzyLanguageParser::STRING);

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
    setState(65);
    match(FuzzyLanguageParser::T__7);
    setState(66);
    fuzzy_set(0);
    setState(67);
    match(FuzzyLanguageParser::T__8);
    setState(68);
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

std::vector<FuzzyLanguageParser::Fuzzy_setContext *> FuzzyLanguageParser::Fuzzy_setContext::fuzzy_set() {
  return getRuleContexts<FuzzyLanguageParser::Fuzzy_setContext>();
}

FuzzyLanguageParser::Fuzzy_setContext *FuzzyLanguageParser::Fuzzy_setContext::fuzzy_set(size_t i) {
  return getRuleContext<FuzzyLanguageParser::Fuzzy_setContext>(i);
}

FuzzyLanguageParser::SelectionContext *FuzzyLanguageParser::Fuzzy_setContext::selection() {
  return getRuleContext<FuzzyLanguageParser::SelectionContext>(0);
}

size_t FuzzyLanguageParser::Fuzzy_setContext::getRuleIndex() const { return FuzzyLanguageParser::RuleFuzzy_set; }

antlrcpp::Any FuzzyLanguageParser::Fuzzy_setContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitFuzzy_set(this);
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
    setState(78);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case FuzzyLanguageParser::T__4: {
        setState(71);
        match(FuzzyLanguageParser::T__4);
        setState(72);
        fuzzy_set(0);
        setState(73);
        match(FuzzyLanguageParser::T__6);
        break;
      }

      case FuzzyLanguageParser::T__11: {
        setState(75);
        match(FuzzyLanguageParser::T__11);
        setState(76);
        fuzzy_set(2);
        break;
      }

      case FuzzyLanguageParser::STRING: {
        setState(77);
        selection();
        break;
      }

      default:
        throw NoViableAltException(this);
    }
    _ctx->stop = _input->LT(-1);
    setState(88);
    _errHandler->sync(this);
    alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 6, _ctx);
    while (alt != 2 && alt != atn::ATN::INVALID_ALT_NUMBER) {
      if (alt == 1) {
        if (!_parseListeners.empty()) triggerExitRuleEvent();
        previousContext = _localctx;
        setState(86);
        _errHandler->sync(this);
        switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 5, _ctx)) {
          case 1: {
            _localctx = _tracker.createInstance<Fuzzy_setContext>(parentContext, parentState);
            pushNewRecursionContext(_localctx, startState, RuleFuzzy_set);
            setState(80);

            if (!(precpred(_ctx, 4))) throw FailedPredicateException(this, "precpred(_ctx, 4)");
            setState(81);
            match(FuzzyLanguageParser::T__9);
            setState(82);
            fuzzy_set(5);
            break;
          }

          case 2: {
            _localctx = _tracker.createInstance<Fuzzy_setContext>(parentContext, parentState);
            pushNewRecursionContext(_localctx, startState, RuleFuzzy_set);
            setState(83);

            if (!(precpred(_ctx, 3))) throw FailedPredicateException(this, "precpred(_ctx, 3)");
            setState(84);
            match(FuzzyLanguageParser::T__10);
            setState(85);
            fuzzy_set(4);
            break;
          }

          default:
            break;
        }
      }
      setState(90);
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

//----------------- SelectionContext ------------------------------------------------------------------

FuzzyLanguageParser::SelectionContext::SelectionContext(ParserRuleContext *parent, size_t invokingState)
    : ParserRuleContext(parent, invokingState) {}

std::vector<FuzzyLanguageParser::NameContext *> FuzzyLanguageParser::SelectionContext::name() {
  return getRuleContexts<FuzzyLanguageParser::NameContext>();
}

FuzzyLanguageParser::NameContext *FuzzyLanguageParser::SelectionContext::name(size_t i) {
  return getRuleContext<FuzzyLanguageParser::NameContext>(i);
}

size_t FuzzyLanguageParser::SelectionContext::getRuleIndex() const { return FuzzyLanguageParser::RuleSelection; }

antlrcpp::Any FuzzyLanguageParser::SelectionContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<FuzzyLanguageVisitor *>(visitor))
    return parserVisitor->visitSelection(this);
  else
    return visitor->visitChildren(this);
}

FuzzyLanguageParser::SelectionContext *FuzzyLanguageParser::selection() {
  SelectionContext *_localctx = _tracker.createInstance<SelectionContext>(_ctx, getState());
  enterRule(_localctx, 14, FuzzyLanguageParser::RuleSelection);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(91);
    name();
    setState(92);
    match(FuzzyLanguageParser::T__12);
    setState(93);
    name();

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

std::vector<std::string> FuzzyLanguageParser::_ruleNames = {
    "rule_file", "fuzzy_variable", "membership_function", "function", "name", "fuzzy_rule", "fuzzy_set", "selection"};

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
      0x3,  0x608b, 0xa72a, 0x8133, 0xb9ed, 0x417c, 0x3be7, 0x7786, 0x5964, 0x3,  0x14, 0x62, 0x4,  0x2,  0x9,  0x2,
      0x4,  0x3,    0x9,    0x3,    0x4,    0x4,    0x9,    0x4,    0x4,    0x5,  0x9,  0x5,  0x4,  0x6,  0x9,  0x6,
      0x4,  0x7,    0x9,    0x7,    0x4,    0x8,    0x9,    0x8,    0x4,    0x9,  0x9,  0x9,  0x3,  0x2,  0x7,  0x2,
      0x14, 0xa,    0x2,    0xc,    0x2,    0xe,    0x2,    0x17,   0xb,    0x2,  0x3,  0x2,  0x7,  0x2,  0x1a, 0xa,
      0x2,  0xc,    0x2,    0xe,    0x2,    0x1d,   0xb,    0x2,    0x3,    0x2,  0x3,  0x2,  0x3,  0x3,  0x3,  0x3,
      0x3,  0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,  0x3,  0x3,  0x3,  0x3,  0x3,  0x3,
      0x3,  0x3,    0x3,    0x3,    0x3,    0x3,    0x6,    0x3,    0x2e,   0xa,  0x3,  0xd,  0x3,  0xe,  0x3,  0x2f,
      0x3,  0x4,    0x3,    0x4,    0x3,    0x4,    0x3,    0x4,    0x3,    0x5,  0x3,  0x5,  0x3,  0x5,  0x3,  0x5,
      0x3,  0x5,    0x7,    0x5,    0x3b,   0xa,    0x5,    0xc,    0x5,    0xe,  0x5,  0x3e, 0xb,  0x5,  0x3,  0x5,
      0x3,  0x5,    0x3,    0x6,    0x3,    0x6,    0x3,    0x7,    0x3,    0x7,  0x3,  0x7,  0x3,  0x7,  0x3,  0x7,
      0x3,  0x8,    0x3,    0x8,    0x3,    0x8,    0x3,    0x8,    0x3,    0x8,  0x3,  0x8,  0x3,  0x8,  0x3,  0x8,
      0x5,  0x8,    0x51,   0xa,    0x8,    0x3,    0x8,    0x3,    0x8,    0x3,  0x8,  0x3,  0x8,  0x3,  0x8,  0x3,
      0x8,  0x7,    0x8,    0x59,   0xa,    0x8,    0xc,    0x8,    0xe,    0x8,  0x5c, 0xb,  0x8,  0x3,  0x9,  0x3,
      0x9,  0x3,    0x9,    0x3,    0x9,    0x3,    0x9,    0x2,    0x3,    0xe,  0xa,  0x2,  0x4,  0x6,  0x8,  0xa,
      0xc,  0xe,    0x10,   0x2,    0x2,    0x2,    0x61,   0x2,    0x15,   0x3,  0x2,  0x2,  0x2,  0x4,  0x20, 0x3,
      0x2,  0x2,    0x2,    0x6,    0x31,   0x3,    0x2,    0x2,    0x2,    0x8,  0x35, 0x3,  0x2,  0x2,  0x2,  0xa,
      0x41, 0x3,    0x2,    0x2,    0x2,    0xc,    0x43,   0x3,    0x2,    0x2,  0x2,  0xe,  0x50, 0x3,  0x2,  0x2,
      0x2,  0x10,   0x5d,   0x3,    0x2,    0x2,    0x2,    0x12,   0x14,   0x5,  0x4,  0x3,  0x2,  0x13, 0x12, 0x3,
      0x2,  0x2,    0x2,    0x14,   0x17,   0x3,    0x2,    0x2,    0x2,    0x15, 0x13, 0x3,  0x2,  0x2,  0x2,  0x15,
      0x16, 0x3,    0x2,    0x2,    0x2,    0x16,   0x1b,   0x3,    0x2,    0x2,  0x2,  0x17, 0x15, 0x3,  0x2,  0x2,
      0x2,  0x18,   0x1a,   0x5,    0xc,    0x7,    0x2,    0x19,   0x18,   0x3,  0x2,  0x2,  0x2,  0x1a, 0x1d, 0x3,
      0x2,  0x2,    0x2,    0x1b,   0x19,   0x3,    0x2,    0x2,    0x2,    0x1b, 0x1c, 0x3,  0x2,  0x2,  0x2,  0x1c,
      0x1e, 0x3,    0x2,    0x2,    0x2,    0x1d,   0x1b,   0x3,    0x2,    0x2,  0x2,  0x1e, 0x1f, 0x7,  0x2,  0x2,
      0x3,  0x1f,   0x3,    0x3,    0x2,    0x2,    0x2,    0x20,   0x21,   0x7,  0x3,  0x2,  0x2,  0x21, 0x22, 0x7,
      0x4,  0x2,    0x2,    0x22,   0x23,   0x7,    0x5,    0x2,    0x2,    0x23, 0x24, 0x7,  0x4,  0x2,  0x2,  0x24,
      0x25, 0x5,    0xa,    0x6,    0x2,    0x25,   0x26,   0x7,    0x6,    0x2,  0x2,  0x26, 0x27, 0x7,  0x4,  0x2,
      0x2,  0x27,   0x28,   0x7,    0x7,    0x2,    0x2,    0x28,   0x29,   0x7,  0x13, 0x2,  0x2,  0x29, 0x2a, 0x7,
      0x8,  0x2,    0x2,    0x2a,   0x2b,   0x7,    0x13,   0x2,    0x2,    0x2b, 0x2d, 0x7,  0x9,  0x2,  0x2,  0x2c,
      0x2e, 0x5,    0x6,    0x4,    0x2,    0x2d,   0x2c,   0x3,    0x2,    0x2,  0x2,  0x2e, 0x2f, 0x3,  0x2,  0x2,
      0x2,  0x2f,   0x2d,   0x3,    0x2,    0x2,    0x2,    0x2f,   0x30,   0x3,  0x2,  0x2,  0x2,  0x30, 0x5,  0x3,
      0x2,  0x2,    0x2,    0x31,   0x32,   0x5,    0xa,    0x6,    0x2,    0x32, 0x33, 0x7,  0x4,  0x2,  0x2,  0x33,
      0x34, 0x5,    0x8,    0x5,    0x2,    0x34,   0x7,    0x3,    0x2,    0x2,  0x2,  0x35, 0x36, 0x7,  0x14, 0x2,
      0x2,  0x36,   0x37,   0x7,    0x7,    0x2,    0x2,    0x37,   0x3c,   0x7,  0x13, 0x2,  0x2,  0x38, 0x39, 0x7,
      0x8,  0x2,    0x2,    0x39,   0x3b,   0x7,    0x13,   0x2,    0x2,    0x3a, 0x38, 0x3,  0x2,  0x2,  0x2,  0x3b,
      0x3e, 0x3,    0x2,    0x2,    0x2,    0x3c,   0x3a,   0x3,    0x2,    0x2,  0x2,  0x3c, 0x3d, 0x3,  0x2,  0x2,
      0x2,  0x3d,   0x3f,   0x3,    0x2,    0x2,    0x2,    0x3e,   0x3c,   0x3,  0x2,  0x2,  0x2,  0x3f, 0x40, 0x7,
      0x9,  0x2,    0x2,    0x40,   0x9,    0x3,    0x2,    0x2,    0x2,    0x41, 0x42, 0x7,  0x12, 0x2,  0x2,  0x42,
      0xb,  0x3,    0x2,    0x2,    0x2,    0x43,   0x44,   0x7,    0xa,    0x2,  0x2,  0x44, 0x45, 0x5,  0xe,  0x8,
      0x2,  0x45,   0x46,   0x7,    0xb,    0x2,    0x2,    0x46,   0x47,   0x5,  0xe,  0x8,  0x2,  0x47, 0xd,  0x3,
      0x2,  0x2,    0x2,    0x48,   0x49,   0x8,    0x8,    0x1,    0x2,    0x49, 0x4a, 0x7,  0x7,  0x2,  0x2,  0x4a,
      0x4b, 0x5,    0xe,    0x8,    0x2,    0x4b,   0x4c,   0x7,    0x9,    0x2,  0x2,  0x4c, 0x51, 0x3,  0x2,  0x2,
      0x2,  0x4d,   0x4e,   0x7,    0xe,    0x2,    0x2,    0x4e,   0x51,   0x5,  0xe,  0x8,  0x4,  0x4f, 0x51, 0x5,
      0x10, 0x9,    0x2,    0x50,   0x48,   0x3,    0x2,    0x2,    0x2,    0x50, 0x4d, 0x3,  0x2,  0x2,  0x2,  0x50,
      0x4f, 0x3,    0x2,    0x2,    0x2,    0x51,   0x5a,   0x3,    0x2,    0x2,  0x2,  0x52, 0x53, 0xc,  0x6,  0x2,
      0x2,  0x53,   0x54,   0x7,    0xc,    0x2,    0x2,    0x54,   0x59,   0x5,  0xe,  0x8,  0x7,  0x55, 0x56, 0xc,
      0x5,  0x2,    0x2,    0x56,   0x57,   0x7,    0xd,    0x2,    0x2,    0x57, 0x59, 0x5,  0xe,  0x8,  0x6,  0x58,
      0x52, 0x3,    0x2,    0x2,    0x2,    0x58,   0x55,   0x3,    0x2,    0x2,  0x2,  0x59, 0x5c, 0x3,  0x2,  0x2,
      0x2,  0x5a,   0x58,   0x3,    0x2,    0x2,    0x2,    0x5a,   0x5b,   0x3,  0x2,  0x2,  0x2,  0x5b, 0xf,  0x3,
      0x2,  0x2,    0x2,    0x5c,   0x5a,   0x3,    0x2,    0x2,    0x2,    0x5d, 0x5e, 0x5,  0xa,  0x6,  0x2,  0x5e,
      0x5f, 0x7,    0xf,    0x2,    0x2,    0x5f,   0x60,   0x5,    0xa,    0x6,  0x2,  0x60, 0x11, 0x3,  0x2,  0x2,
      0x2,  0x9,    0x15,   0x1b,   0x2f,   0x3c,   0x50,   0x58,   0x5a,
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


// Generated from
// AutoPas/src/autopas/tuning/tuningStrategy/fuzzyTuning\FuzzyLanguage.g4 by
// ANTLR 4.9.1

#include "FuzzyLanguageParser.h"

#include "FuzzyLanguageVisitor.h"

using namespace antlrcpp;
using namespace AutopasGeneratedFuzzyRuleSyntax;
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

antlrcpp::Any FuzzyLanguageParser::SettingsContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::Linguistic_variableContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::Fuzzy_termContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::FunctionContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::Output_mappingContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::Output_entryContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::Pattern_mappingContext::accept(tree::ParseTreeVisitor *visitor) {
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

antlrcpp::Any FuzzyLanguageParser::Configuration_patternContext::accept(tree::ParseTreeVisitor *visitor) {
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

std::vector<std::string> FuzzyLanguageParser::_ruleNames = {"rule_file",
                                                            "settings",
                                                            "linguistic_variable",
                                                            "fuzzy_term",
                                                            "function",
                                                            "fuzzy_rule",
                                                            "fuzzy_set",
                                                            "output_mapping",
                                                            "output_entry",
                                                            "pattern_mapping",
                                                            "configuration_pattern"};

std::vector<std::string> FuzzyLanguageParser::_literalNames = {"",         "'FuzzySystemSettings'",
                                                               "':'",      "'FuzzyVariable'",
                                                               "'domain'", "'range'",
                                                               "'('",      "','",
                                                               "')'",      "'if'",
                                                               "'then'",   "'&&'",
                                                               "'||'",     "'!'",
                                                               "'=='",     "'OutputMapping'",
                                                               "'=>'",     "'['",
                                                               "'='",      "']'"};

std::vector<std::string> FuzzyLanguageParser::_symbolicNames = {
    "", "", "", "", "", "", "", "",   "",        "",       "",       "",          "",
    "", "", "", "", "", "", "", "WS", "COMMENT", "STRING", "NUMBER", "IDENTIFIER"};

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
      0x3,  0x608b, 0xa72a, 0x8133, 0xb9ed, 0x417c, 0x3be7, 0x7786, 0x5964, 0x3,  0x1a, 0x98, 0x4,  0x2,  0x9,  0x2,
      0x4,  0x3,    0x9,    0x3,    0x4,    0x4,    0x9,    0x4,    0x4,    0x5,  0x9,  0x5,  0x4,  0x6,  0x9,  0x6,
      0x4,  0x7,    0x9,    0x7,    0x4,    0x8,    0x9,    0x8,    0x4,    0x9,  0x9,  0x9,  0x4,  0xa,  0x9,  0xa,
      0x4,  0xb,    0x9,    0xb,    0x4,    0xc,    0x9,    0xc,    0x3,    0x2,  0x3,  0x2,  0x7,  0x2,  0x1b, 0xa,
      0x2,  0xc,    0x2,    0xe,    0x2,    0x1e,   0xb,    0x2,    0x3,    0x2,  0x3,  0x2,  0x7,  0x2,  0x22, 0xa,
      0x2,  0xc,    0x2,    0xe,    0x2,    0x25,   0xb,    0x2,    0x3,    0x2,  0x3,  0x2,  0x3,  0x3,  0x3,  0x3,
      0x3,  0x3,    0x3,    0x3,    0x3,    0x3,    0x7,    0x3,    0x2e,   0xa,  0x3,  0xc,  0x3,  0xe,  0x3,  0x31,
      0xb,  0x3,    0x3,    0x4,    0x3,    0x4,    0x3,    0x4,    0x3,    0x4,  0x3,  0x4,  0x3,  0x4,  0x3,  0x4,
      0x3,  0x4,    0x3,    0x4,    0x3,    0x4,    0x3,    0x4,    0x3,    0x4,  0x3,  0x4,  0x6,  0x4,  0x40, 0xa,
      0x4,  0xd,    0x4,    0xe,    0x4,    0x41,   0x3,    0x5,    0x3,    0x5,  0x3,  0x5,  0x3,  0x5,  0x3,  0x6,
      0x3,  0x6,    0x3,    0x6,    0x3,    0x6,    0x3,    0x6,    0x7,    0x6,  0x4d, 0xa,  0x6,  0xc,  0x6,  0xe,
      0x6,  0x50,   0xb,    0x6,    0x3,    0x6,    0x3,    0x6,    0x3,    0x7,  0x3,  0x7,  0x3,  0x7,  0x3,  0x7,
      0x3,  0x7,    0x3,    0x8,    0x3,    0x8,    0x3,    0x8,    0x3,    0x8,  0x3,  0x8,  0x3,  0x8,  0x3,  0x8,
      0x3,  0x8,    0x3,    0x8,    0x3,    0x8,    0x5,    0x8,    0x63,   0xa,  0x8,  0x3,  0x8,  0x3,  0x8,  0x3,
      0x8,  0x3,    0x8,    0x3,    0x8,    0x3,    0x8,    0x7,    0x8,    0x6b, 0xa,  0x8,  0xc,  0x8,  0xe,  0x8,
      0x6e, 0xb,    0x8,    0x3,    0x9,    0x3,    0x9,    0x3,    0x9,    0x6,  0x9,  0x73, 0xa,  0x9,  0xd,  0x9,
      0xe,  0x9,    0x74,   0x3,    0xa,    0x3,    0xa,    0x3,    0xa,    0x6,  0xa,  0x7a, 0xa,  0xa,  0xd,  0xa,
      0xe,  0xa,    0x7b,   0x3,    0xb,    0x3,    0xb,    0x3,    0xb,    0x3,  0xb,  0x3,  0xb,  0x7,  0xb,  0x83,
      0xa,  0xb,    0xc,    0xb,    0xe,    0xb,    0x86,   0xb,    0xb,    0x3,  0xc,  0x3,  0xc,  0x3,  0xc,  0x3,
      0xc,  0x3,    0xc,    0x3,    0xc,    0x3,    0xc,    0x3,    0xc,    0x3,  0xc,  0x7,  0xc,  0x91, 0xa,  0xc,
      0xc,  0xc,    0xe,    0xc,    0x94,   0xb,    0xc,    0x3,    0xc,    0x3,  0xc,  0x3,  0xc,  0x2,  0x3,  0xe,
      0xd,  0x2,    0x4,    0x6,    0x8,    0xa,    0xc,    0xe,    0x10,   0x12, 0x14, 0x16, 0x2,  0x2,  0x2,  0x99,
      0x2,  0x18,   0x3,    0x2,    0x2,    0x2,    0x4,    0x28,   0x3,    0x2,  0x2,  0x2,  0x6,  0x32, 0x3,  0x2,
      0x2,  0x2,    0x8,    0x43,   0x3,    0x2,    0x2,    0x2,    0xa,    0x47, 0x3,  0x2,  0x2,  0x2,  0xc,  0x53,
      0x3,  0x2,    0x2,    0x2,    0xe,    0x62,   0x3,    0x2,    0x2,    0x2,  0x10, 0x6f, 0x3,  0x2,  0x2,  0x2,
      0x12, 0x76,   0x3,    0x2,    0x2,    0x2,    0x14,   0x7d,   0x3,    0x2,  0x2,  0x2,  0x16, 0x87, 0x3,  0x2,
      0x2,  0x2,    0x18,   0x1c,   0x5,    0x4,    0x3,    0x2,    0x19,   0x1b, 0x5,  0x6,  0x4,  0x2,  0x1a, 0x19,
      0x3,  0x2,    0x2,    0x2,    0x1b,   0x1e,   0x3,    0x2,    0x2,    0x2,  0x1c, 0x1a, 0x3,  0x2,  0x2,  0x2,
      0x1c, 0x1d,   0x3,    0x2,    0x2,    0x2,    0x1d,   0x1f,   0x3,    0x2,  0x2,  0x2,  0x1e, 0x1c, 0x3,  0x2,
      0x2,  0x2,    0x1f,   0x23,   0x5,    0x10,   0x9,    0x2,    0x20,   0x22, 0x5,  0xc,  0x7,  0x2,  0x21, 0x20,
      0x3,  0x2,    0x2,    0x2,    0x22,   0x25,   0x3,    0x2,    0x2,    0x2,  0x23, 0x21, 0x3,  0x2,  0x2,  0x2,
      0x23, 0x24,   0x3,    0x2,    0x2,    0x2,    0x24,   0x26,   0x3,    0x2,  0x2,  0x2,  0x25, 0x23, 0x3,  0x2,
      0x2,  0x2,    0x26,   0x27,   0x7,    0x2,    0x2,    0x3,    0x27,   0x3,  0x3,  0x2,  0x2,  0x2,  0x28, 0x29,
      0x7,  0x3,    0x2,    0x2,    0x29,   0x2f,   0x7,    0x4,    0x2,    0x2,  0x2a, 0x2b, 0x7,  0x1a, 0x2,  0x2,
      0x2b, 0x2c,   0x7,    0x4,    0x2,    0x2,    0x2c,   0x2e,   0x7,    0x18, 0x2,  0x2,  0x2d, 0x2a, 0x3,  0x2,
      0x2,  0x2,    0x2e,   0x31,   0x3,    0x2,    0x2,    0x2,    0x2f,   0x2d, 0x3,  0x2,  0x2,  0x2,  0x2f, 0x30,
      0x3,  0x2,    0x2,    0x2,    0x30,   0x5,    0x3,    0x2,    0x2,    0x2,  0x31, 0x2f, 0x3,  0x2,  0x2,  0x2,
      0x32, 0x33,   0x7,    0x5,    0x2,    0x2,    0x33,   0x34,   0x7,    0x4,  0x2,  0x2,  0x34, 0x35, 0x7,  0x6,
      0x2,  0x2,    0x35,   0x36,   0x7,    0x4,    0x2,    0x2,    0x36,   0x37, 0x7,  0x18, 0x2,  0x2,  0x37, 0x38,
      0x7,  0x7,    0x2,    0x2,    0x38,   0x39,   0x7,    0x4,    0x2,    0x2,  0x39, 0x3a, 0x7,  0x8,  0x2,  0x2,
      0x3a, 0x3b,   0x7,    0x19,   0x2,    0x2,    0x3b,   0x3c,   0x7,    0x9,  0x2,  0x2,  0x3c, 0x3d, 0x7,  0x19,
      0x2,  0x2,    0x3d,   0x3f,   0x7,    0xa,    0x2,    0x2,    0x3e,   0x40, 0x5,  0x8,  0x5,  0x2,  0x3f, 0x3e,
      0x3,  0x2,    0x2,    0x2,    0x40,   0x41,   0x3,    0x2,    0x2,    0x2,  0x41, 0x3f, 0x3,  0x2,  0x2,  0x2,
      0x41, 0x42,   0x3,    0x2,    0x2,    0x2,    0x42,   0x7,    0x3,    0x2,  0x2,  0x2,  0x43, 0x44, 0x7,  0x18,
      0x2,  0x2,    0x44,   0x45,   0x7,    0x4,    0x2,    0x2,    0x45,   0x46, 0x5,  0xa,  0x6,  0x2,  0x46, 0x9,
      0x3,  0x2,    0x2,    0x2,    0x47,   0x48,   0x7,    0x1a,   0x2,    0x2,  0x48, 0x49, 0x7,  0x8,  0x2,  0x2,
      0x49, 0x4e,   0x7,    0x19,   0x2,    0x2,    0x4a,   0x4b,   0x7,    0x9,  0x2,  0x2,  0x4b, 0x4d, 0x7,  0x19,
      0x2,  0x2,    0x4c,   0x4a,   0x3,    0x2,    0x2,    0x2,    0x4d,   0x50, 0x3,  0x2,  0x2,  0x2,  0x4e, 0x4c,
      0x3,  0x2,    0x2,    0x2,    0x4e,   0x4f,   0x3,    0x2,    0x2,    0x2,  0x4f, 0x51, 0x3,  0x2,  0x2,  0x2,
      0x50, 0x4e,   0x3,    0x2,    0x2,    0x2,    0x51,   0x52,   0x7,    0xa,  0x2,  0x2,  0x52, 0xb,  0x3,  0x2,
      0x2,  0x2,    0x53,   0x54,   0x7,    0xb,    0x2,    0x2,    0x54,   0x55, 0x5,  0xe,  0x8,  0x2,  0x55, 0x56,
      0x7,  0xc,    0x2,    0x2,    0x56,   0x57,   0x5,    0xe,    0x8,    0x2,  0x57, 0xd,  0x3,  0x2,  0x2,  0x2,
      0x58, 0x59,   0x8,    0x8,    0x1,    0x2,    0x59,   0x5a,   0x7,    0x8,  0x2,  0x2,  0x5a, 0x5b, 0x5,  0xe,
      0x8,  0x2,    0x5b,   0x5c,   0x7,    0xa,    0x2,    0x2,    0x5c,   0x63, 0x3,  0x2,  0x2,  0x2,  0x5d, 0x5e,
      0x7,  0xf,    0x2,    0x2,    0x5e,   0x63,   0x5,    0xe,    0x8,    0x4,  0x5f, 0x60, 0x7,  0x18, 0x2,  0x2,
      0x60, 0x61,   0x7,    0x10,   0x2,    0x2,    0x61,   0x63,   0x7,    0x18, 0x2,  0x2,  0x62, 0x58, 0x3,  0x2,
      0x2,  0x2,    0x62,   0x5d,   0x3,    0x2,    0x2,    0x2,    0x62,   0x5f, 0x3,  0x2,  0x2,  0x2,  0x63, 0x6c,
      0x3,  0x2,    0x2,    0x2,    0x64,   0x65,   0xc,    0x6,    0x2,    0x2,  0x65, 0x66, 0x7,  0xd,  0x2,  0x2,
      0x66, 0x6b,   0x5,    0xe,    0x8,    0x7,    0x67,   0x68,   0xc,    0x5,  0x2,  0x2,  0x68, 0x69, 0x7,  0xe,
      0x2,  0x2,    0x69,   0x6b,   0x5,    0xe,    0x8,    0x6,    0x6a,   0x64, 0x3,  0x2,  0x2,  0x2,  0x6a, 0x67,
      0x3,  0x2,    0x2,    0x2,    0x6b,   0x6e,   0x3,    0x2,    0x2,    0x2,  0x6c, 0x6a, 0x3,  0x2,  0x2,  0x2,
      0x6c, 0x6d,   0x3,    0x2,    0x2,    0x2,    0x6d,   0xf,    0x3,    0x2,  0x2,  0x2,  0x6e, 0x6c, 0x3,  0x2,
      0x2,  0x2,    0x6f,   0x70,   0x7,    0x11,   0x2,    0x2,    0x70,   0x72, 0x7,  0x4,  0x2,  0x2,  0x71, 0x73,
      0x5,  0x12,   0xa,    0x2,    0x72,   0x71,   0x3,    0x2,    0x2,    0x2,  0x73, 0x74, 0x3,  0x2,  0x2,  0x2,
      0x74, 0x72,   0x3,    0x2,    0x2,    0x2,    0x74,   0x75,   0x3,    0x2,  0x2,  0x2,  0x75, 0x11, 0x3,  0x2,
      0x2,  0x2,    0x76,   0x77,   0x7,    0x18,   0x2,    0x2,    0x77,   0x79, 0x7,  0x4,  0x2,  0x2,  0x78, 0x7a,
      0x5,  0x14,   0xb,    0x2,    0x79,   0x78,   0x3,    0x2,    0x2,    0x2,  0x7a, 0x7b, 0x3,  0x2,  0x2,  0x2,
      0x7b, 0x79,   0x3,    0x2,    0x2,    0x2,    0x7b,   0x7c,   0x3,    0x2,  0x2,  0x2,  0x7c, 0x13, 0x3,  0x2,
      0x2,  0x2,    0x7d,   0x7e,   0x7,    0x19,   0x2,    0x2,    0x7e,   0x7f, 0x7,  0x12, 0x2,  0x2,  0x7f, 0x84,
      0x5,  0x16,   0xc,    0x2,    0x80,   0x81,   0x7,    0x9,    0x2,    0x2,  0x81, 0x83, 0x5,  0x16, 0xc,  0x2,
      0x82, 0x80,   0x3,    0x2,    0x2,    0x2,    0x83,   0x86,   0x3,    0x2,  0x2,  0x2,  0x84, 0x82, 0x3,  0x2,
      0x2,  0x2,    0x84,   0x85,   0x3,    0x2,    0x2,    0x2,    0x85,   0x15, 0x3,  0x2,  0x2,  0x2,  0x86, 0x84,
      0x3,  0x2,    0x2,    0x2,    0x87,   0x88,   0x7,    0x13,   0x2,    0x2,  0x88, 0x89, 0x7,  0x1a, 0x2,  0x2,
      0x89, 0x8a,   0x7,    0x14,   0x2,    0x2,    0x8a,   0x8b,   0x7,    0x18, 0x2,  0x2,  0x8b, 0x92, 0x3,  0x2,
      0x2,  0x2,    0x8c,   0x8d,   0x7,    0x9,    0x2,    0x2,    0x8d,   0x8e, 0x7,  0x1a, 0x2,  0x2,  0x8e, 0x8f,
      0x7,  0x14,   0x2,    0x2,    0x8f,   0x91,   0x7,    0x18,   0x2,    0x2,  0x90, 0x8c, 0x3,  0x2,  0x2,  0x2,
      0x91, 0x94,   0x3,    0x2,    0x2,    0x2,    0x92,   0x90,   0x3,    0x2,  0x2,  0x2,  0x92, 0x93, 0x3,  0x2,
      0x2,  0x2,    0x93,   0x95,   0x3,    0x2,    0x2,    0x2,    0x94,   0x92, 0x3,  0x2,  0x2,  0x2,  0x95, 0x96,
      0x7,  0x15,   0x2,    0x2,    0x96,   0x17,   0x3,    0x2,    0x2,    0x2,  0xe,  0x1c, 0x23, 0x2f, 0x41, 0x4e,
      0x62, 0x6a,   0x6c,   0x74,   0x7b,   0x84,   0x92,
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

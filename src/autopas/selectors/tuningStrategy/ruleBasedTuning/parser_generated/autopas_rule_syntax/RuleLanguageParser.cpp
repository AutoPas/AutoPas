
// Generated from /home/tobias/AutoPas2/src/autopas/selectors/tuningStrategy/ruleBasedTuning/RuleLanguage.g4 by ANTLR 4.9.1


#include "RuleLanguageListener.h"
#include "RuleLanguageVisitor.h"

#include "RuleLanguageParser.h"


using namespace antlrcpp;
using namespace autopas::rule_syntax;
using namespace antlr4;

RuleLanguageParser::RuleLanguageParser(TokenStream *input) : Parser(input) {
  _interpreter = new atn::ParserATNSimulator(this, _atn, _decisionToDFA, _sharedContextCache);
}

RuleLanguageParser::~RuleLanguageParser() {
  delete _interpreter;
}

std::string RuleLanguageParser::getGrammarFileName() const {
  return "RuleLanguage.g4";
}

const std::vector<std::string>& RuleLanguageParser::getRuleNames() const {
  return _ruleNames;
}

dfa::Vocabulary& RuleLanguageParser::getVocabulary() const {
  return _vocabulary;
}


//----------------- ProgramContext ------------------------------------------------------------------

RuleLanguageParser::ProgramContext::ProgramContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<RuleLanguageParser::StatementContext *> RuleLanguageParser::ProgramContext::statement() {
  return getRuleContexts<RuleLanguageParser::StatementContext>();
}

RuleLanguageParser::StatementContext* RuleLanguageParser::ProgramContext::statement(size_t i) {
  return getRuleContext<RuleLanguageParser::StatementContext>(i);
}


size_t RuleLanguageParser::ProgramContext::getRuleIndex() const {
  return RuleLanguageParser::RuleProgram;
}

void RuleLanguageParser::ProgramContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterProgram(this);
}

void RuleLanguageParser::ProgramContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitProgram(this);
}


antlrcpp::Any RuleLanguageParser::ProgramContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitProgram(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::ProgramContext* RuleLanguageParser::program() {
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
    setState(23); 
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(22);
      statement();
      setState(25); 
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << RuleLanguageParser::T__0)
      | (1ULL << RuleLanguageParser::T__4)
      | (1ULL << RuleLanguageParser::T__5)
      | (1ULL << RuleLanguageParser::T__8))) != 0));
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- LiteralContext ------------------------------------------------------------------

RuleLanguageParser::LiteralContext::LiteralContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Traversal_opt() {
  return getToken(RuleLanguageParser::Traversal_opt, 0);
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Container_opt() {
  return getToken(RuleLanguageParser::Container_opt, 0);
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Load_estimator_opt() {
  return getToken(RuleLanguageParser::Load_estimator_opt, 0);
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Data_layout_opt() {
  return getToken(RuleLanguageParser::Data_layout_opt, 0);
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Newton3_opt() {
  return getToken(RuleLanguageParser::Newton3_opt, 0);
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Unsigned_val() {
  return getToken(RuleLanguageParser::Unsigned_val, 0);
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Bool_val() {
  return getToken(RuleLanguageParser::Bool_val, 0);
}


size_t RuleLanguageParser::LiteralContext::getRuleIndex() const {
  return RuleLanguageParser::RuleLiteral;
}

void RuleLanguageParser::LiteralContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterLiteral(this);
}

void RuleLanguageParser::LiteralContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitLiteral(this);
}


antlrcpp::Any RuleLanguageParser::LiteralContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitLiteral(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::LiteralContext* RuleLanguageParser::literal() {
  LiteralContext *_localctx = _tracker.createInstance<LiteralContext>(_ctx, getState());
  enterRule(_localctx, 2, RuleLanguageParser::RuleLiteral);
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
    setState(27);
    _la = _input->LA(1);
    if (!((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << RuleLanguageParser::Container_opt)
      | (1ULL << RuleLanguageParser::Traversal_opt)
      | (1ULL << RuleLanguageParser::Load_estimator_opt)
      | (1ULL << RuleLanguageParser::Data_layout_opt)
      | (1ULL << RuleLanguageParser::Newton3_opt)
      | (1ULL << RuleLanguageParser::Bool_val)
      | (1ULL << RuleLanguageParser::Unsigned_val))) != 0))) {
    _errHandler->recoverInline(this);
    }
    else {
      _errHandler->reportMatch(this);
      consume();
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Define_listContext ------------------------------------------------------------------

RuleLanguageParser::Define_listContext::Define_listContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* RuleLanguageParser::Define_listContext::Variable_name() {
  return getToken(RuleLanguageParser::Variable_name, 0);
}

std::vector<RuleLanguageParser::LiteralContext *> RuleLanguageParser::Define_listContext::literal() {
  return getRuleContexts<RuleLanguageParser::LiteralContext>();
}

RuleLanguageParser::LiteralContext* RuleLanguageParser::Define_listContext::literal(size_t i) {
  return getRuleContext<RuleLanguageParser::LiteralContext>(i);
}


size_t RuleLanguageParser::Define_listContext::getRuleIndex() const {
  return RuleLanguageParser::RuleDefine_list;
}

void RuleLanguageParser::Define_listContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterDefine_list(this);
}

void RuleLanguageParser::Define_listContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitDefine_list(this);
}


antlrcpp::Any RuleLanguageParser::Define_listContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitDefine_list(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Define_listContext* RuleLanguageParser::define_list() {
  Define_listContext *_localctx = _tracker.createInstance<Define_listContext>(_ctx, getState());
  enterRule(_localctx, 4, RuleLanguageParser::RuleDefine_list);
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
    match(RuleLanguageParser::T__0);
    setState(30);
    match(RuleLanguageParser::Variable_name);
    setState(31);
    match(RuleLanguageParser::T__1);
    setState(32);
    literal();
    setState(37);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == RuleLanguageParser::T__2) {
      setState(33);
      match(RuleLanguageParser::T__2);
      setState(34);
      literal();
      setState(39);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(40);
    match(RuleLanguageParser::T__3);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- DefineContext ------------------------------------------------------------------

RuleLanguageParser::DefineContext::DefineContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* RuleLanguageParser::DefineContext::Variable_name() {
  return getToken(RuleLanguageParser::Variable_name, 0);
}

RuleLanguageParser::LiteralContext* RuleLanguageParser::DefineContext::literal() {
  return getRuleContext<RuleLanguageParser::LiteralContext>(0);
}


size_t RuleLanguageParser::DefineContext::getRuleIndex() const {
  return RuleLanguageParser::RuleDefine;
}

void RuleLanguageParser::DefineContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterDefine(this);
}

void RuleLanguageParser::DefineContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitDefine(this);
}


antlrcpp::Any RuleLanguageParser::DefineContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitDefine(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::DefineContext* RuleLanguageParser::define() {
  DefineContext *_localctx = _tracker.createInstance<DefineContext>(_ctx, getState());
  enterRule(_localctx, 6, RuleLanguageParser::RuleDefine);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(42);
    match(RuleLanguageParser::T__4);
    setState(43);
    match(RuleLanguageParser::Variable_name);
    setState(44);
    match(RuleLanguageParser::T__1);
    setState(45);
    literal();
    setState(46);
    match(RuleLanguageParser::T__3);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- VariableContext ------------------------------------------------------------------

RuleLanguageParser::VariableContext::VariableContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* RuleLanguageParser::VariableContext::Variable_name() {
  return getToken(RuleLanguageParser::Variable_name, 0);
}


size_t RuleLanguageParser::VariableContext::getRuleIndex() const {
  return RuleLanguageParser::RuleVariable;
}

void RuleLanguageParser::VariableContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterVariable(this);
}

void RuleLanguageParser::VariableContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitVariable(this);
}


antlrcpp::Any RuleLanguageParser::VariableContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitVariable(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::VariableContext* RuleLanguageParser::variable() {
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
    setState(48);
    match(RuleLanguageParser::Variable_name);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ExpressionContext ------------------------------------------------------------------

RuleLanguageParser::ExpressionContext::ExpressionContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

RuleLanguageParser::VariableContext* RuleLanguageParser::ExpressionContext::variable() {
  return getRuleContext<RuleLanguageParser::VariableContext>(0);
}

RuleLanguageParser::LiteralContext* RuleLanguageParser::ExpressionContext::literal() {
  return getRuleContext<RuleLanguageParser::LiteralContext>(0);
}

std::vector<RuleLanguageParser::ExpressionContext *> RuleLanguageParser::ExpressionContext::expression() {
  return getRuleContexts<RuleLanguageParser::ExpressionContext>();
}

RuleLanguageParser::ExpressionContext* RuleLanguageParser::ExpressionContext::expression(size_t i) {
  return getRuleContext<RuleLanguageParser::ExpressionContext>(i);
}

tree::TerminalNode* RuleLanguageParser::ExpressionContext::Binary_op() {
  return getToken(RuleLanguageParser::Binary_op, 0);
}


size_t RuleLanguageParser::ExpressionContext::getRuleIndex() const {
  return RuleLanguageParser::RuleExpression;
}

void RuleLanguageParser::ExpressionContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterExpression(this);
}

void RuleLanguageParser::ExpressionContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitExpression(this);
}


antlrcpp::Any RuleLanguageParser::ExpressionContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitExpression(this);
  else
    return visitor->visitChildren(this);
}


RuleLanguageParser::ExpressionContext* RuleLanguageParser::expression() {
   return expression(0);
}

RuleLanguageParser::ExpressionContext* RuleLanguageParser::expression(int precedence) {
  ParserRuleContext *parentContext = _ctx;
  size_t parentState = getState();
  RuleLanguageParser::ExpressionContext *_localctx = _tracker.createInstance<ExpressionContext>(_ctx, parentState);
  RuleLanguageParser::ExpressionContext *previousContext = _localctx;
  (void)previousContext; // Silence compiler, in case the context is not used by generated code.
  size_t startState = 10;
  enterRecursionRule(_localctx, 10, RuleLanguageParser::RuleExpression, precedence);

    

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
    setState(53);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case RuleLanguageParser::Variable_name: {
        setState(51);
        variable();
        break;
      }

      case RuleLanguageParser::Container_opt:
      case RuleLanguageParser::Traversal_opt:
      case RuleLanguageParser::Load_estimator_opt:
      case RuleLanguageParser::Data_layout_opt:
      case RuleLanguageParser::Newton3_opt:
      case RuleLanguageParser::Bool_val:
      case RuleLanguageParser::Unsigned_val: {
        setState(52);
        literal();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
    _ctx->stop = _input->LT(-1);
    setState(60);
    _errHandler->sync(this);
    alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 3, _ctx);
    while (alt != 2 && alt != atn::ATN::INVALID_ALT_NUMBER) {
      if (alt == 1) {
        if (!_parseListeners.empty())
          triggerExitRuleEvent();
        previousContext = _localctx;
        _localctx = _tracker.createInstance<ExpressionContext>(parentContext, parentState);
        pushNewRecursionContext(_localctx, startState, RuleExpression);
        setState(55);

        if (!(precpred(_ctx, 3))) throw FailedPredicateException(this, "precpred(_ctx, 3)");
        setState(56);
        match(RuleLanguageParser::Binary_op);
        setState(57);
        expression(4); 
      }
      setState(62);
      _errHandler->sync(this);
      alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 3, _ctx);
    }
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }
  return _localctx;
}

//----------------- Property_valueContext ------------------------------------------------------------------

RuleLanguageParser::Property_valueContext::Property_valueContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* RuleLanguageParser::Property_valueContext::Variable_name() {
  return getToken(RuleLanguageParser::Variable_name, 0);
}

RuleLanguageParser::LiteralContext* RuleLanguageParser::Property_valueContext::literal() {
  return getRuleContext<RuleLanguageParser::LiteralContext>(0);
}


size_t RuleLanguageParser::Property_valueContext::getRuleIndex() const {
  return RuleLanguageParser::RuleProperty_value;
}

void RuleLanguageParser::Property_valueContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterProperty_value(this);
}

void RuleLanguageParser::Property_valueContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitProperty_value(this);
}


antlrcpp::Any RuleLanguageParser::Property_valueContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitProperty_value(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Property_valueContext* RuleLanguageParser::property_value() {
  Property_valueContext *_localctx = _tracker.createInstance<Property_valueContext>(_ctx, getState());
  enterRule(_localctx, 12, RuleLanguageParser::RuleProperty_value);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(65);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case RuleLanguageParser::Variable_name: {
        enterOuterAlt(_localctx, 1);
        setState(63);
        match(RuleLanguageParser::Variable_name);
        break;
      }

      case RuleLanguageParser::Container_opt:
      case RuleLanguageParser::Traversal_opt:
      case RuleLanguageParser::Load_estimator_opt:
      case RuleLanguageParser::Data_layout_opt:
      case RuleLanguageParser::Newton3_opt:
      case RuleLanguageParser::Bool_val:
      case RuleLanguageParser::Unsigned_val: {
        enterOuterAlt(_localctx, 2);
        setState(64);
        literal();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Configuration_patternContext ------------------------------------------------------------------

RuleLanguageParser::Configuration_patternContext::Configuration_patternContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<tree::TerminalNode *> RuleLanguageParser::Configuration_patternContext::Configuration_property() {
  return getTokens(RuleLanguageParser::Configuration_property);
}

tree::TerminalNode* RuleLanguageParser::Configuration_patternContext::Configuration_property(size_t i) {
  return getToken(RuleLanguageParser::Configuration_property, i);
}

std::vector<RuleLanguageParser::Property_valueContext *> RuleLanguageParser::Configuration_patternContext::property_value() {
  return getRuleContexts<RuleLanguageParser::Property_valueContext>();
}

RuleLanguageParser::Property_valueContext* RuleLanguageParser::Configuration_patternContext::property_value(size_t i) {
  return getRuleContext<RuleLanguageParser::Property_valueContext>(i);
}


size_t RuleLanguageParser::Configuration_patternContext::getRuleIndex() const {
  return RuleLanguageParser::RuleConfiguration_pattern;
}

void RuleLanguageParser::Configuration_patternContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterConfiguration_pattern(this);
}

void RuleLanguageParser::Configuration_patternContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitConfiguration_pattern(this);
}


antlrcpp::Any RuleLanguageParser::Configuration_patternContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitConfiguration_pattern(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Configuration_patternContext* RuleLanguageParser::configuration_pattern() {
  Configuration_patternContext *_localctx = _tracker.createInstance<Configuration_patternContext>(_ctx, getState());
  enterRule(_localctx, 14, RuleLanguageParser::RuleConfiguration_pattern);
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
    setState(67);
    match(RuleLanguageParser::T__5);

    setState(68);
    match(RuleLanguageParser::Configuration_property);
    setState(69);
    match(RuleLanguageParser::T__1);
    setState(70);
    property_value();
    setState(78);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == RuleLanguageParser::T__2) {
      setState(72);
      match(RuleLanguageParser::T__2);
      setState(73);
      match(RuleLanguageParser::Configuration_property);
      setState(74);
      match(RuleLanguageParser::T__1);
      setState(75);
      property_value();
      setState(80);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(81);
    match(RuleLanguageParser::T__6);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Configuration_orderContext ------------------------------------------------------------------

RuleLanguageParser::Configuration_orderContext::Configuration_orderContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<RuleLanguageParser::Configuration_patternContext *> RuleLanguageParser::Configuration_orderContext::configuration_pattern() {
  return getRuleContexts<RuleLanguageParser::Configuration_patternContext>();
}

RuleLanguageParser::Configuration_patternContext* RuleLanguageParser::Configuration_orderContext::configuration_pattern(size_t i) {
  return getRuleContext<RuleLanguageParser::Configuration_patternContext>(i);
}


size_t RuleLanguageParser::Configuration_orderContext::getRuleIndex() const {
  return RuleLanguageParser::RuleConfiguration_order;
}

void RuleLanguageParser::Configuration_orderContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterConfiguration_order(this);
}

void RuleLanguageParser::Configuration_orderContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitConfiguration_order(this);
}


antlrcpp::Any RuleLanguageParser::Configuration_orderContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitConfiguration_order(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Configuration_orderContext* RuleLanguageParser::configuration_order() {
  Configuration_orderContext *_localctx = _tracker.createInstance<Configuration_orderContext>(_ctx, getState());
  enterRule(_localctx, 16, RuleLanguageParser::RuleConfiguration_order);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(83);
    configuration_pattern();
    setState(84);
    match(RuleLanguageParser::T__7);
    setState(85);
    configuration_pattern();
    setState(86);
    match(RuleLanguageParser::T__3);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- StatementContext ------------------------------------------------------------------

RuleLanguageParser::StatementContext::StatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

RuleLanguageParser::Define_listContext* RuleLanguageParser::StatementContext::define_list() {
  return getRuleContext<RuleLanguageParser::Define_listContext>(0);
}

RuleLanguageParser::DefineContext* RuleLanguageParser::StatementContext::define() {
  return getRuleContext<RuleLanguageParser::DefineContext>(0);
}

RuleLanguageParser::If_statementContext* RuleLanguageParser::StatementContext::if_statement() {
  return getRuleContext<RuleLanguageParser::If_statementContext>(0);
}

RuleLanguageParser::Configuration_orderContext* RuleLanguageParser::StatementContext::configuration_order() {
  return getRuleContext<RuleLanguageParser::Configuration_orderContext>(0);
}


size_t RuleLanguageParser::StatementContext::getRuleIndex() const {
  return RuleLanguageParser::RuleStatement;
}

void RuleLanguageParser::StatementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterStatement(this);
}

void RuleLanguageParser::StatementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitStatement(this);
}


antlrcpp::Any RuleLanguageParser::StatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitStatement(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::StatementContext* RuleLanguageParser::statement() {
  StatementContext *_localctx = _tracker.createInstance<StatementContext>(_ctx, getState());
  enterRule(_localctx, 18, RuleLanguageParser::RuleStatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(92);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case RuleLanguageParser::T__0: {
        enterOuterAlt(_localctx, 1);
        setState(88);
        define_list();
        break;
      }

      case RuleLanguageParser::T__4: {
        enterOuterAlt(_localctx, 2);
        setState(89);
        define();
        break;
      }

      case RuleLanguageParser::T__8: {
        enterOuterAlt(_localctx, 3);
        setState(90);
        if_statement();
        break;
      }

      case RuleLanguageParser::T__5: {
        enterOuterAlt(_localctx, 4);
        setState(91);
        configuration_order();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- If_statementContext ------------------------------------------------------------------

RuleLanguageParser::If_statementContext::If_statementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

RuleLanguageParser::ExpressionContext* RuleLanguageParser::If_statementContext::expression() {
  return getRuleContext<RuleLanguageParser::ExpressionContext>(0);
}

std::vector<RuleLanguageParser::StatementContext *> RuleLanguageParser::If_statementContext::statement() {
  return getRuleContexts<RuleLanguageParser::StatementContext>();
}

RuleLanguageParser::StatementContext* RuleLanguageParser::If_statementContext::statement(size_t i) {
  return getRuleContext<RuleLanguageParser::StatementContext>(i);
}


size_t RuleLanguageParser::If_statementContext::getRuleIndex() const {
  return RuleLanguageParser::RuleIf_statement;
}

void RuleLanguageParser::If_statementContext::enterRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->enterIf_statement(this);
}

void RuleLanguageParser::If_statementContext::exitRule(tree::ParseTreeListener *listener) {
  auto parserListener = dynamic_cast<RuleLanguageListener *>(listener);
  if (parserListener != nullptr)
    parserListener->exitIf_statement(this);
}


antlrcpp::Any RuleLanguageParser::If_statementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitIf_statement(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::If_statementContext* RuleLanguageParser::if_statement() {
  If_statementContext *_localctx = _tracker.createInstance<If_statementContext>(_ctx, getState());
  enterRule(_localctx, 20, RuleLanguageParser::RuleIf_statement);
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
    setState(94);
    match(RuleLanguageParser::T__8);
    setState(95);
    expression(0);
    setState(96);
    match(RuleLanguageParser::T__9);
    setState(98); 
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(97);
      statement();
      setState(100); 
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << RuleLanguageParser::T__0)
      | (1ULL << RuleLanguageParser::T__4)
      | (1ULL << RuleLanguageParser::T__5)
      | (1ULL << RuleLanguageParser::T__8))) != 0));
    setState(102);
    match(RuleLanguageParser::T__10);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

bool RuleLanguageParser::sempred(RuleContext *context, size_t ruleIndex, size_t predicateIndex) {
  switch (ruleIndex) {
    case 5: return expressionSempred(dynamic_cast<ExpressionContext *>(context), predicateIndex);

  default:
    break;
  }
  return true;
}

bool RuleLanguageParser::expressionSempred(ExpressionContext *_localctx, size_t predicateIndex) {
  switch (predicateIndex) {
    case 0: return precpred(_ctx, 3);

  default:
    break;
  }
  return true;
}

// Static vars and initialization.
std::vector<dfa::DFA> RuleLanguageParser::_decisionToDFA;
atn::PredictionContextCache RuleLanguageParser::_sharedContextCache;

// We own the ATN which in turn owns the ATN states.
atn::ATN RuleLanguageParser::_atn;
std::vector<uint16_t> RuleLanguageParser::_serializedATN;

std::vector<std::string> RuleLanguageParser::_ruleNames = {
  "program", "literal", "define_list", "define", "variable", "expression", 
  "property_value", "configuration_pattern", "configuration_order", "statement", 
  "if_statement"
};

std::vector<std::string> RuleLanguageParser::_literalNames = {
  "", "'define_list'", "'='", "','", "';'", "'define'", "'['", "']'", "'>='", 
  "'if'", "':'", "'endif'", "", "", "'None'"
};

std::vector<std::string> RuleLanguageParser::_symbolicNames = {
  "", "", "", "", "", "", "", "", "", "", "", "", "Container_opt", "Traversal_opt", 
  "Load_estimator_opt", "Data_layout_opt", "Newton3_opt", "Bool_val", "Binary_op", 
  "Configuration_property", "DIGIT", "Unsigned_val", "Variable_name"
};

dfa::Vocabulary RuleLanguageParser::_vocabulary(_literalNames, _symbolicNames);

std::vector<std::string> RuleLanguageParser::_tokenNames;

RuleLanguageParser::Initializer::Initializer() {
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
    0x3, 0x608b, 0xa72a, 0x8133, 0xb9ed, 0x417c, 0x3be7, 0x7786, 0x5964, 
    0x3, 0x18, 0x6b, 0x4, 0x2, 0x9, 0x2, 0x4, 0x3, 0x9, 0x3, 0x4, 0x4, 0x9, 
    0x4, 0x4, 0x5, 0x9, 0x5, 0x4, 0x6, 0x9, 0x6, 0x4, 0x7, 0x9, 0x7, 0x4, 
    0x8, 0x9, 0x8, 0x4, 0x9, 0x9, 0x9, 0x4, 0xa, 0x9, 0xa, 0x4, 0xb, 0x9, 
    0xb, 0x4, 0xc, 0x9, 0xc, 0x3, 0x2, 0x6, 0x2, 0x1a, 0xa, 0x2, 0xd, 0x2, 
    0xe, 0x2, 0x1b, 0x3, 0x3, 0x3, 0x3, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 
    0x4, 0x3, 0x4, 0x3, 0x4, 0x7, 0x4, 0x26, 0xa, 0x4, 0xc, 0x4, 0xe, 0x4, 
    0x29, 0xb, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x5, 0x3, 0x5, 0x3, 0x5, 0x3, 
    0x5, 0x3, 0x5, 0x3, 0x5, 0x3, 0x6, 0x3, 0x6, 0x3, 0x7, 0x3, 0x7, 0x3, 
    0x7, 0x5, 0x7, 0x38, 0xa, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x7, 0x7, 
    0x3d, 0xa, 0x7, 0xc, 0x7, 0xe, 0x7, 0x40, 0xb, 0x7, 0x3, 0x8, 0x3, 0x8, 
    0x5, 0x8, 0x44, 0xa, 0x8, 0x3, 0x9, 0x3, 0x9, 0x3, 0x9, 0x3, 0x9, 0x3, 
    0x9, 0x3, 0x9, 0x3, 0x9, 0x3, 0x9, 0x3, 0x9, 0x7, 0x9, 0x4f, 0xa, 0x9, 
    0xc, 0x9, 0xe, 0x9, 0x52, 0xb, 0x9, 0x3, 0x9, 0x3, 0x9, 0x3, 0xa, 0x3, 
    0xa, 0x3, 0xa, 0x3, 0xa, 0x3, 0xa, 0x3, 0xb, 0x3, 0xb, 0x3, 0xb, 0x3, 
    0xb, 0x5, 0xb, 0x5f, 0xa, 0xb, 0x3, 0xc, 0x3, 0xc, 0x3, 0xc, 0x3, 0xc, 
    0x6, 0xc, 0x65, 0xa, 0xc, 0xd, 0xc, 0xe, 0xc, 0x66, 0x3, 0xc, 0x3, 0xc, 
    0x3, 0xc, 0x2, 0x3, 0xc, 0xd, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe, 0x10, 
    0x12, 0x14, 0x16, 0x2, 0x3, 0x4, 0x2, 0xe, 0x13, 0x17, 0x17, 0x2, 0x69, 
    0x2, 0x19, 0x3, 0x2, 0x2, 0x2, 0x4, 0x1d, 0x3, 0x2, 0x2, 0x2, 0x6, 0x1f, 
    0x3, 0x2, 0x2, 0x2, 0x8, 0x2c, 0x3, 0x2, 0x2, 0x2, 0xa, 0x32, 0x3, 0x2, 
    0x2, 0x2, 0xc, 0x37, 0x3, 0x2, 0x2, 0x2, 0xe, 0x43, 0x3, 0x2, 0x2, 0x2, 
    0x10, 0x45, 0x3, 0x2, 0x2, 0x2, 0x12, 0x55, 0x3, 0x2, 0x2, 0x2, 0x14, 
    0x5e, 0x3, 0x2, 0x2, 0x2, 0x16, 0x60, 0x3, 0x2, 0x2, 0x2, 0x18, 0x1a, 
    0x5, 0x14, 0xb, 0x2, 0x19, 0x18, 0x3, 0x2, 0x2, 0x2, 0x1a, 0x1b, 0x3, 
    0x2, 0x2, 0x2, 0x1b, 0x19, 0x3, 0x2, 0x2, 0x2, 0x1b, 0x1c, 0x3, 0x2, 
    0x2, 0x2, 0x1c, 0x3, 0x3, 0x2, 0x2, 0x2, 0x1d, 0x1e, 0x9, 0x2, 0x2, 
    0x2, 0x1e, 0x5, 0x3, 0x2, 0x2, 0x2, 0x1f, 0x20, 0x7, 0x3, 0x2, 0x2, 
    0x20, 0x21, 0x7, 0x18, 0x2, 0x2, 0x21, 0x22, 0x7, 0x4, 0x2, 0x2, 0x22, 
    0x27, 0x5, 0x4, 0x3, 0x2, 0x23, 0x24, 0x7, 0x5, 0x2, 0x2, 0x24, 0x26, 
    0x5, 0x4, 0x3, 0x2, 0x25, 0x23, 0x3, 0x2, 0x2, 0x2, 0x26, 0x29, 0x3, 
    0x2, 0x2, 0x2, 0x27, 0x25, 0x3, 0x2, 0x2, 0x2, 0x27, 0x28, 0x3, 0x2, 
    0x2, 0x2, 0x28, 0x2a, 0x3, 0x2, 0x2, 0x2, 0x29, 0x27, 0x3, 0x2, 0x2, 
    0x2, 0x2a, 0x2b, 0x7, 0x6, 0x2, 0x2, 0x2b, 0x7, 0x3, 0x2, 0x2, 0x2, 
    0x2c, 0x2d, 0x7, 0x7, 0x2, 0x2, 0x2d, 0x2e, 0x7, 0x18, 0x2, 0x2, 0x2e, 
    0x2f, 0x7, 0x4, 0x2, 0x2, 0x2f, 0x30, 0x5, 0x4, 0x3, 0x2, 0x30, 0x31, 
    0x7, 0x6, 0x2, 0x2, 0x31, 0x9, 0x3, 0x2, 0x2, 0x2, 0x32, 0x33, 0x7, 
    0x18, 0x2, 0x2, 0x33, 0xb, 0x3, 0x2, 0x2, 0x2, 0x34, 0x35, 0x8, 0x7, 
    0x1, 0x2, 0x35, 0x38, 0x5, 0xa, 0x6, 0x2, 0x36, 0x38, 0x5, 0x4, 0x3, 
    0x2, 0x37, 0x34, 0x3, 0x2, 0x2, 0x2, 0x37, 0x36, 0x3, 0x2, 0x2, 0x2, 
    0x38, 0x3e, 0x3, 0x2, 0x2, 0x2, 0x39, 0x3a, 0xc, 0x5, 0x2, 0x2, 0x3a, 
    0x3b, 0x7, 0x14, 0x2, 0x2, 0x3b, 0x3d, 0x5, 0xc, 0x7, 0x6, 0x3c, 0x39, 
    0x3, 0x2, 0x2, 0x2, 0x3d, 0x40, 0x3, 0x2, 0x2, 0x2, 0x3e, 0x3c, 0x3, 
    0x2, 0x2, 0x2, 0x3e, 0x3f, 0x3, 0x2, 0x2, 0x2, 0x3f, 0xd, 0x3, 0x2, 
    0x2, 0x2, 0x40, 0x3e, 0x3, 0x2, 0x2, 0x2, 0x41, 0x44, 0x7, 0x18, 0x2, 
    0x2, 0x42, 0x44, 0x5, 0x4, 0x3, 0x2, 0x43, 0x41, 0x3, 0x2, 0x2, 0x2, 
    0x43, 0x42, 0x3, 0x2, 0x2, 0x2, 0x44, 0xf, 0x3, 0x2, 0x2, 0x2, 0x45, 
    0x46, 0x7, 0x8, 0x2, 0x2, 0x46, 0x47, 0x7, 0x15, 0x2, 0x2, 0x47, 0x48, 
    0x7, 0x4, 0x2, 0x2, 0x48, 0x49, 0x5, 0xe, 0x8, 0x2, 0x49, 0x50, 0x3, 
    0x2, 0x2, 0x2, 0x4a, 0x4b, 0x7, 0x5, 0x2, 0x2, 0x4b, 0x4c, 0x7, 0x15, 
    0x2, 0x2, 0x4c, 0x4d, 0x7, 0x4, 0x2, 0x2, 0x4d, 0x4f, 0x5, 0xe, 0x8, 
    0x2, 0x4e, 0x4a, 0x3, 0x2, 0x2, 0x2, 0x4f, 0x52, 0x3, 0x2, 0x2, 0x2, 
    0x50, 0x4e, 0x3, 0x2, 0x2, 0x2, 0x50, 0x51, 0x3, 0x2, 0x2, 0x2, 0x51, 
    0x53, 0x3, 0x2, 0x2, 0x2, 0x52, 0x50, 0x3, 0x2, 0x2, 0x2, 0x53, 0x54, 
    0x7, 0x9, 0x2, 0x2, 0x54, 0x11, 0x3, 0x2, 0x2, 0x2, 0x55, 0x56, 0x5, 
    0x10, 0x9, 0x2, 0x56, 0x57, 0x7, 0xa, 0x2, 0x2, 0x57, 0x58, 0x5, 0x10, 
    0x9, 0x2, 0x58, 0x59, 0x7, 0x6, 0x2, 0x2, 0x59, 0x13, 0x3, 0x2, 0x2, 
    0x2, 0x5a, 0x5f, 0x5, 0x6, 0x4, 0x2, 0x5b, 0x5f, 0x5, 0x8, 0x5, 0x2, 
    0x5c, 0x5f, 0x5, 0x16, 0xc, 0x2, 0x5d, 0x5f, 0x5, 0x12, 0xa, 0x2, 0x5e, 
    0x5a, 0x3, 0x2, 0x2, 0x2, 0x5e, 0x5b, 0x3, 0x2, 0x2, 0x2, 0x5e, 0x5c, 
    0x3, 0x2, 0x2, 0x2, 0x5e, 0x5d, 0x3, 0x2, 0x2, 0x2, 0x5f, 0x15, 0x3, 
    0x2, 0x2, 0x2, 0x60, 0x61, 0x7, 0xb, 0x2, 0x2, 0x61, 0x62, 0x5, 0xc, 
    0x7, 0x2, 0x62, 0x64, 0x7, 0xc, 0x2, 0x2, 0x63, 0x65, 0x5, 0x14, 0xb, 
    0x2, 0x64, 0x63, 0x3, 0x2, 0x2, 0x2, 0x65, 0x66, 0x3, 0x2, 0x2, 0x2, 
    0x66, 0x64, 0x3, 0x2, 0x2, 0x2, 0x66, 0x67, 0x3, 0x2, 0x2, 0x2, 0x67, 
    0x68, 0x3, 0x2, 0x2, 0x2, 0x68, 0x69, 0x7, 0xd, 0x2, 0x2, 0x69, 0x17, 
    0x3, 0x2, 0x2, 0x2, 0xa, 0x1b, 0x27, 0x37, 0x3e, 0x43, 0x50, 0x5e, 0x66, 
  };

  atn::ATNDeserializer deserializer;
  _atn = deserializer.deserialize(_serializedATN);

  size_t count = _atn.getNumberOfDecisions();
  _decisionToDFA.reserve(count);
  for (size_t i = 0; i < count; i++) { 
    _decisionToDFA.emplace_back(_atn.getDecisionState(i), i);
  }
}

RuleLanguageParser::Initializer RuleLanguageParser::_init;

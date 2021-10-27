
// Generated from /home/tobias/AutoPas2/src/autopas/selectors/tuningStrategy/ruleBasedTuning/RuleLanguage.g4 by ANTLR 4.9.1


#include "RuleLanguageVisitor.h"

#include "RuleLanguageParser.h"


using namespace antlrcpp;
using namespace autopas_generated_rule_syntax;
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
    setState(25); 
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(24);
      statement();
      setState(27); 
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << RuleLanguageParser::T__1)
      | (1ULL << RuleLanguageParser::T__16)
      | (1ULL << RuleLanguageParser::T__17)
      | (1ULL << RuleLanguageParser::T__21))) != 0));
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- Unsigned_valContext ------------------------------------------------------------------

RuleLanguageParser::Unsigned_valContext::Unsigned_valContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* RuleLanguageParser::Unsigned_valContext::Unsigned_val() {
  return getToken(RuleLanguageParser::Unsigned_val, 0);
}

tree::TerminalNode* RuleLanguageParser::Unsigned_valContext::DIGIT() {
  return getToken(RuleLanguageParser::DIGIT, 0);
}


size_t RuleLanguageParser::Unsigned_valContext::getRuleIndex() const {
  return RuleLanguageParser::RuleUnsigned_val;
}


antlrcpp::Any RuleLanguageParser::Unsigned_valContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitUnsigned_val(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Unsigned_valContext* RuleLanguageParser::unsigned_val() {
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

RuleLanguageParser::Unsigned_valContext* RuleLanguageParser::LiteralContext::unsigned_val() {
  return getRuleContext<RuleLanguageParser::Unsigned_valContext>(0);
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Double_val() {
  return getToken(RuleLanguageParser::Double_val, 0);
}

tree::TerminalNode* RuleLanguageParser::LiteralContext::Bool_val() {
  return getToken(RuleLanguageParser::Bool_val, 0);
}


size_t RuleLanguageParser::LiteralContext::getRuleIndex() const {
  return RuleLanguageParser::RuleLiteral;
}


antlrcpp::Any RuleLanguageParser::LiteralContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitLiteral(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::LiteralContext* RuleLanguageParser::literal() {
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


antlrcpp::Any RuleLanguageParser::Define_listContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitDefine_list(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Define_listContext* RuleLanguageParser::define_list() {
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
    setState(66);
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

std::vector<RuleLanguageParser::ExpressionContext *> RuleLanguageParser::ExpressionContext::expression() {
  return getRuleContexts<RuleLanguageParser::ExpressionContext>();
}

RuleLanguageParser::ExpressionContext* RuleLanguageParser::ExpressionContext::expression(size_t i) {
  return getRuleContext<RuleLanguageParser::ExpressionContext>(i);
}

RuleLanguageParser::LiteralContext* RuleLanguageParser::ExpressionContext::literal() {
  return getRuleContext<RuleLanguageParser::LiteralContext>(0);
}

RuleLanguageParser::VariableContext* RuleLanguageParser::ExpressionContext::variable() {
  return getRuleContext<RuleLanguageParser::VariableContext>(0);
}


size_t RuleLanguageParser::ExpressionContext::getRuleIndex() const {
  return RuleLanguageParser::RuleExpression;
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
        dynamic_cast<ExpressionContext *>(_localctx)->op = match(RuleLanguageParser::T__11);
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
        if (!_parseListeners.empty())
          triggerExitRuleEvent();
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
          dynamic_cast<ExpressionContext *>(_localctx)->op = _input->LT(1);
          _la = _input->LA(1);
          if (!(_la == RuleLanguageParser::T__5

          || _la == RuleLanguageParser::T__6)) {
            dynamic_cast<ExpressionContext *>(_localctx)->op = _errHandler->recoverInline(this);
          }
          else {
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
          dynamic_cast<ExpressionContext *>(_localctx)->op = _input->LT(1);
          _la = _input->LA(1);
          if (!(_la == RuleLanguageParser::T__7

          || _la == RuleLanguageParser::T__8)) {
            dynamic_cast<ExpressionContext *>(_localctx)->op = _errHandler->recoverInline(this);
          }
          else {
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
          dynamic_cast<ExpressionContext *>(_localctx)->op = _input->LT(1);
          _la = _input->LA(1);
          if (!(_la == RuleLanguageParser::T__9

          || _la == RuleLanguageParser::T__10)) {
            dynamic_cast<ExpressionContext *>(_localctx)->op = _errHandler->recoverInline(this);
          }
          else {
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
          dynamic_cast<ExpressionContext *>(_localctx)->op = _input->LT(1);
          _la = _input->LA(1);
          if (!(_la == RuleLanguageParser::T__12

          || _la == RuleLanguageParser::T__13)) {
            dynamic_cast<ExpressionContext *>(_localctx)->op = _errHandler->recoverInline(this);
          }
          else {
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

RuleLanguageParser::ExpressionContext* RuleLanguageParser::DefineContext::expression() {
  return getRuleContext<RuleLanguageParser::ExpressionContext>(0);
}


size_t RuleLanguageParser::DefineContext::getRuleIndex() const {
  return RuleLanguageParser::RuleDefine;
}


antlrcpp::Any RuleLanguageParser::DefineContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitDefine(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::DefineContext* RuleLanguageParser::define() {
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


antlrcpp::Any RuleLanguageParser::Property_valueContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitProperty_value(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Property_valueContext* RuleLanguageParser::property_value() {
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


antlrcpp::Any RuleLanguageParser::Configuration_patternContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitConfiguration_pattern(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Configuration_patternContext* RuleLanguageParser::configuration_pattern() {
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

std::vector<tree::TerminalNode *> RuleLanguageParser::Configuration_orderContext::Configuration_property() {
  return getTokens(RuleLanguageParser::Configuration_property);
}

tree::TerminalNode* RuleLanguageParser::Configuration_orderContext::Configuration_property(size_t i) {
  return getToken(RuleLanguageParser::Configuration_property, i);
}


size_t RuleLanguageParser::Configuration_orderContext::getRuleIndex() const {
  return RuleLanguageParser::RuleConfiguration_order;
}


antlrcpp::Any RuleLanguageParser::Configuration_orderContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitConfiguration_order(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::Configuration_orderContext* RuleLanguageParser::configuration_order() {
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


antlrcpp::Any RuleLanguageParser::StatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitStatement(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::StatementContext* RuleLanguageParser::statement() {
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


antlrcpp::Any RuleLanguageParser::If_statementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<RuleLanguageVisitor*>(visitor))
    return parserVisitor->visitIf_statement(this);
  else
    return visitor->visitChildren(this);
}

RuleLanguageParser::If_statementContext* RuleLanguageParser::if_statement() {
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
    } while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << RuleLanguageParser::T__1)
      | (1ULL << RuleLanguageParser::T__16)
      | (1ULL << RuleLanguageParser::T__17)
      | (1ULL << RuleLanguageParser::T__21))) != 0));
    setState(152);
    match(RuleLanguageParser::T__23);
   
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
    case 0: return precpred(_ctx, 8);
    case 1: return precpred(_ctx, 7);
    case 2: return precpred(_ctx, 6);
    case 3: return precpred(_ctx, 4);

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
  "program", "unsigned_val", "literal", "define_list", "variable", "expression", 
  "define", "property_value", "configuration_pattern", "configuration_order", 
  "statement", "if_statement"
};

std::vector<std::string> RuleLanguageParser::_literalNames = {
  "", "'\"'", "'define_list'", "'='", "','", "';'", "'*'", "'/'", "'+'", 
  "'-'", "'>'", "'<'", "'not'", "'and'", "'or'", "'('", "')'", "'define'", 
  "'['", "']'", "'>='", "'with same'", "'if'", "':'", "'endif'", "", "", 
  "", "", "'None'"
};

std::vector<std::string> RuleLanguageParser::_symbolicNames = {
  "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
  "", "", "", "", "", "", "", "COMMENT", "WS", "Container_opt", "Traversal_opt", 
  "Load_estimator_opt", "Data_layout_opt", "Newton3_opt", "Bool_val", "Configuration_property", 
  "DIGIT", "Unsigned_val", "Double_val", "Variable_name"
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
    0x3, 0x27, 0x9d, 0x4, 0x2, 0x9, 0x2, 0x4, 0x3, 0x9, 0x3, 0x4, 0x4, 0x9, 
    0x4, 0x4, 0x5, 0x9, 0x5, 0x4, 0x6, 0x9, 0x6, 0x4, 0x7, 0x9, 0x7, 0x4, 
    0x8, 0x9, 0x8, 0x4, 0x9, 0x9, 0x9, 0x4, 0xa, 0x9, 0xa, 0x4, 0xb, 0x9, 
    0xb, 0x4, 0xc, 0x9, 0xc, 0x4, 0xd, 0x9, 0xd, 0x3, 0x2, 0x6, 0x2, 0x1c, 
    0xa, 0x2, 0xd, 0x2, 0xe, 0x2, 0x1d, 0x3, 0x3, 0x3, 0x3, 0x3, 0x4, 0x3, 
    0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 
    0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 
    0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x3, 0x4, 0x5, 0x4, 0x36, 0xa, 0x4, 
    0x3, 0x5, 0x3, 0x5, 0x3, 0x5, 0x3, 0x5, 0x3, 0x5, 0x3, 0x5, 0x7, 0x5, 
    0x3e, 0xa, 0x5, 0xc, 0x5, 0xe, 0x5, 0x41, 0xb, 0x5, 0x3, 0x5, 0x3, 0x5, 
    0x3, 0x6, 0x3, 0x6, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 
    0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x5, 0x7, 0x50, 0xa, 0x7, 0x3, 
    0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 
    0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x3, 0x7, 0x7, 0x7, 0x5e, 0xa, 0x7, 
    0xc, 0x7, 0xe, 0x7, 0x61, 0xb, 0x7, 0x3, 0x8, 0x3, 0x8, 0x3, 0x8, 0x3, 
    0x8, 0x3, 0x8, 0x3, 0x8, 0x3, 0x9, 0x3, 0x9, 0x5, 0x9, 0x6b, 0xa, 0x9, 
    0x3, 0xa, 0x3, 0xa, 0x3, 0xa, 0x3, 0xa, 0x3, 0xa, 0x3, 0xa, 0x3, 0xa, 
    0x3, 0xa, 0x3, 0xa, 0x7, 0xa, 0x76, 0xa, 0xa, 0xc, 0xa, 0xe, 0xa, 0x79, 
    0xb, 0xa, 0x3, 0xa, 0x3, 0xa, 0x3, 0xb, 0x3, 0xb, 0x3, 0xb, 0x3, 0xb, 
    0x3, 0xb, 0x3, 0xb, 0x3, 0xb, 0x7, 0xb, 0x84, 0xa, 0xb, 0xc, 0xb, 0xe, 
    0xb, 0x87, 0xb, 0xb, 0x5, 0xb, 0x89, 0xa, 0xb, 0x3, 0xb, 0x3, 0xb, 0x3, 
    0xc, 0x3, 0xc, 0x3, 0xc, 0x3, 0xc, 0x5, 0xc, 0x91, 0xa, 0xc, 0x3, 0xd, 
    0x3, 0xd, 0x3, 0xd, 0x3, 0xd, 0x6, 0xd, 0x97, 0xa, 0xd, 0xd, 0xd, 0xe, 
    0xd, 0x98, 0x3, 0xd, 0x3, 0xd, 0x3, 0xd, 0x2, 0x3, 0xc, 0xe, 0x2, 0x4, 
    0x6, 0x8, 0xa, 0xc, 0xe, 0x10, 0x12, 0x14, 0x16, 0x18, 0x2, 0x7, 0x3, 
    0x2, 0x24, 0x25, 0x3, 0x2, 0x8, 0x9, 0x3, 0x2, 0xa, 0xb, 0x3, 0x2, 0xc, 
    0xd, 0x3, 0x2, 0xf, 0x10, 0x2, 0xa8, 0x2, 0x1b, 0x3, 0x2, 0x2, 0x2, 
    0x4, 0x1f, 0x3, 0x2, 0x2, 0x2, 0x6, 0x35, 0x3, 0x2, 0x2, 0x2, 0x8, 0x37, 
    0x3, 0x2, 0x2, 0x2, 0xa, 0x44, 0x3, 0x2, 0x2, 0x2, 0xc, 0x4f, 0x3, 0x2, 
    0x2, 0x2, 0xe, 0x62, 0x3, 0x2, 0x2, 0x2, 0x10, 0x6a, 0x3, 0x2, 0x2, 
    0x2, 0x12, 0x6c, 0x3, 0x2, 0x2, 0x2, 0x14, 0x7c, 0x3, 0x2, 0x2, 0x2, 
    0x16, 0x90, 0x3, 0x2, 0x2, 0x2, 0x18, 0x92, 0x3, 0x2, 0x2, 0x2, 0x1a, 
    0x1c, 0x5, 0x16, 0xc, 0x2, 0x1b, 0x1a, 0x3, 0x2, 0x2, 0x2, 0x1c, 0x1d, 
    0x3, 0x2, 0x2, 0x2, 0x1d, 0x1b, 0x3, 0x2, 0x2, 0x2, 0x1d, 0x1e, 0x3, 
    0x2, 0x2, 0x2, 0x1e, 0x3, 0x3, 0x2, 0x2, 0x2, 0x1f, 0x20, 0x9, 0x2, 
    0x2, 0x2, 0x20, 0x5, 0x3, 0x2, 0x2, 0x2, 0x21, 0x22, 0x7, 0x3, 0x2, 
    0x2, 0x22, 0x23, 0x7, 0x1e, 0x2, 0x2, 0x23, 0x36, 0x7, 0x3, 0x2, 0x2, 
    0x24, 0x25, 0x7, 0x3, 0x2, 0x2, 0x25, 0x26, 0x7, 0x1d, 0x2, 0x2, 0x26, 
    0x36, 0x7, 0x3, 0x2, 0x2, 0x27, 0x28, 0x7, 0x3, 0x2, 0x2, 0x28, 0x29, 
    0x7, 0x1f, 0x2, 0x2, 0x29, 0x36, 0x7, 0x3, 0x2, 0x2, 0x2a, 0x2b, 0x7, 
    0x3, 0x2, 0x2, 0x2b, 0x2c, 0x7, 0x20, 0x2, 0x2, 0x2c, 0x36, 0x7, 0x3, 
    0x2, 0x2, 0x2d, 0x2e, 0x7, 0x3, 0x2, 0x2, 0x2e, 0x2f, 0x7, 0x21, 0x2, 
    0x2, 0x2f, 0x36, 0x7, 0x3, 0x2, 0x2, 0x30, 0x36, 0x5, 0x4, 0x3, 0x2, 
    0x31, 0x36, 0x7, 0x26, 0x2, 0x2, 0x32, 0x33, 0x7, 0x3, 0x2, 0x2, 0x33, 
    0x34, 0x7, 0x22, 0x2, 0x2, 0x34, 0x36, 0x7, 0x3, 0x2, 0x2, 0x35, 0x21, 
    0x3, 0x2, 0x2, 0x2, 0x35, 0x24, 0x3, 0x2, 0x2, 0x2, 0x35, 0x27, 0x3, 
    0x2, 0x2, 0x2, 0x35, 0x2a, 0x3, 0x2, 0x2, 0x2, 0x35, 0x2d, 0x3, 0x2, 
    0x2, 0x2, 0x35, 0x30, 0x3, 0x2, 0x2, 0x2, 0x35, 0x31, 0x3, 0x2, 0x2, 
    0x2, 0x35, 0x32, 0x3, 0x2, 0x2, 0x2, 0x36, 0x7, 0x3, 0x2, 0x2, 0x2, 
    0x37, 0x38, 0x7, 0x4, 0x2, 0x2, 0x38, 0x39, 0x7, 0x27, 0x2, 0x2, 0x39, 
    0x3a, 0x7, 0x5, 0x2, 0x2, 0x3a, 0x3f, 0x5, 0x6, 0x4, 0x2, 0x3b, 0x3c, 
    0x7, 0x6, 0x2, 0x2, 0x3c, 0x3e, 0x5, 0x6, 0x4, 0x2, 0x3d, 0x3b, 0x3, 
    0x2, 0x2, 0x2, 0x3e, 0x41, 0x3, 0x2, 0x2, 0x2, 0x3f, 0x3d, 0x3, 0x2, 
    0x2, 0x2, 0x3f, 0x40, 0x3, 0x2, 0x2, 0x2, 0x40, 0x42, 0x3, 0x2, 0x2, 
    0x2, 0x41, 0x3f, 0x3, 0x2, 0x2, 0x2, 0x42, 0x43, 0x7, 0x7, 0x2, 0x2, 
    0x43, 0x9, 0x3, 0x2, 0x2, 0x2, 0x44, 0x45, 0x7, 0x27, 0x2, 0x2, 0x45, 
    0xb, 0x3, 0x2, 0x2, 0x2, 0x46, 0x47, 0x8, 0x7, 0x1, 0x2, 0x47, 0x48, 
    0x7, 0xe, 0x2, 0x2, 0x48, 0x50, 0x5, 0xc, 0x7, 0x7, 0x49, 0x50, 0x5, 
    0x6, 0x4, 0x2, 0x4a, 0x50, 0x5, 0xa, 0x6, 0x2, 0x4b, 0x4c, 0x7, 0x11, 
    0x2, 0x2, 0x4c, 0x4d, 0x5, 0xc, 0x7, 0x2, 0x4d, 0x4e, 0x7, 0x12, 0x2, 
    0x2, 0x4e, 0x50, 0x3, 0x2, 0x2, 0x2, 0x4f, 0x46, 0x3, 0x2, 0x2, 0x2, 
    0x4f, 0x49, 0x3, 0x2, 0x2, 0x2, 0x4f, 0x4a, 0x3, 0x2, 0x2, 0x2, 0x4f, 
    0x4b, 0x3, 0x2, 0x2, 0x2, 0x50, 0x5f, 0x3, 0x2, 0x2, 0x2, 0x51, 0x52, 
    0xc, 0xa, 0x2, 0x2, 0x52, 0x53, 0x9, 0x3, 0x2, 0x2, 0x53, 0x5e, 0x5, 
    0xc, 0x7, 0xb, 0x54, 0x55, 0xc, 0x9, 0x2, 0x2, 0x55, 0x56, 0x9, 0x4, 
    0x2, 0x2, 0x56, 0x5e, 0x5, 0xc, 0x7, 0xa, 0x57, 0x58, 0xc, 0x8, 0x2, 
    0x2, 0x58, 0x59, 0x9, 0x5, 0x2, 0x2, 0x59, 0x5e, 0x5, 0xc, 0x7, 0x9, 
    0x5a, 0x5b, 0xc, 0x6, 0x2, 0x2, 0x5b, 0x5c, 0x9, 0x6, 0x2, 0x2, 0x5c, 
    0x5e, 0x5, 0xc, 0x7, 0x7, 0x5d, 0x51, 0x3, 0x2, 0x2, 0x2, 0x5d, 0x54, 
    0x3, 0x2, 0x2, 0x2, 0x5d, 0x57, 0x3, 0x2, 0x2, 0x2, 0x5d, 0x5a, 0x3, 
    0x2, 0x2, 0x2, 0x5e, 0x61, 0x3, 0x2, 0x2, 0x2, 0x5f, 0x5d, 0x3, 0x2, 
    0x2, 0x2, 0x5f, 0x60, 0x3, 0x2, 0x2, 0x2, 0x60, 0xd, 0x3, 0x2, 0x2, 
    0x2, 0x61, 0x5f, 0x3, 0x2, 0x2, 0x2, 0x62, 0x63, 0x7, 0x13, 0x2, 0x2, 
    0x63, 0x64, 0x7, 0x27, 0x2, 0x2, 0x64, 0x65, 0x7, 0x5, 0x2, 0x2, 0x65, 
    0x66, 0x5, 0xc, 0x7, 0x2, 0x66, 0x67, 0x7, 0x7, 0x2, 0x2, 0x67, 0xf, 
    0x3, 0x2, 0x2, 0x2, 0x68, 0x6b, 0x7, 0x27, 0x2, 0x2, 0x69, 0x6b, 0x5, 
    0x6, 0x4, 0x2, 0x6a, 0x68, 0x3, 0x2, 0x2, 0x2, 0x6a, 0x69, 0x3, 0x2, 
    0x2, 0x2, 0x6b, 0x11, 0x3, 0x2, 0x2, 0x2, 0x6c, 0x6d, 0x7, 0x14, 0x2, 
    0x2, 0x6d, 0x6e, 0x7, 0x23, 0x2, 0x2, 0x6e, 0x6f, 0x7, 0x5, 0x2, 0x2, 
    0x6f, 0x70, 0x5, 0x10, 0x9, 0x2, 0x70, 0x77, 0x3, 0x2, 0x2, 0x2, 0x71, 
    0x72, 0x7, 0x6, 0x2, 0x2, 0x72, 0x73, 0x7, 0x23, 0x2, 0x2, 0x73, 0x74, 
    0x7, 0x5, 0x2, 0x2, 0x74, 0x76, 0x5, 0x10, 0x9, 0x2, 0x75, 0x71, 0x3, 
    0x2, 0x2, 0x2, 0x76, 0x79, 0x3, 0x2, 0x2, 0x2, 0x77, 0x75, 0x3, 0x2, 
    0x2, 0x2, 0x77, 0x78, 0x3, 0x2, 0x2, 0x2, 0x78, 0x7a, 0x3, 0x2, 0x2, 
    0x2, 0x79, 0x77, 0x3, 0x2, 0x2, 0x2, 0x7a, 0x7b, 0x7, 0x15, 0x2, 0x2, 
    0x7b, 0x13, 0x3, 0x2, 0x2, 0x2, 0x7c, 0x7d, 0x5, 0x12, 0xa, 0x2, 0x7d, 
    0x7e, 0x7, 0x16, 0x2, 0x2, 0x7e, 0x88, 0x5, 0x12, 0xa, 0x2, 0x7f, 0x80, 
    0x7, 0x17, 0x2, 0x2, 0x80, 0x85, 0x7, 0x23, 0x2, 0x2, 0x81, 0x82, 0x7, 
    0x6, 0x2, 0x2, 0x82, 0x84, 0x7, 0x23, 0x2, 0x2, 0x83, 0x81, 0x3, 0x2, 
    0x2, 0x2, 0x84, 0x87, 0x3, 0x2, 0x2, 0x2, 0x85, 0x83, 0x3, 0x2, 0x2, 
    0x2, 0x85, 0x86, 0x3, 0x2, 0x2, 0x2, 0x86, 0x89, 0x3, 0x2, 0x2, 0x2, 
    0x87, 0x85, 0x3, 0x2, 0x2, 0x2, 0x88, 0x7f, 0x3, 0x2, 0x2, 0x2, 0x88, 
    0x89, 0x3, 0x2, 0x2, 0x2, 0x89, 0x8a, 0x3, 0x2, 0x2, 0x2, 0x8a, 0x8b, 
    0x7, 0x7, 0x2, 0x2, 0x8b, 0x15, 0x3, 0x2, 0x2, 0x2, 0x8c, 0x91, 0x5, 
    0x8, 0x5, 0x2, 0x8d, 0x91, 0x5, 0xe, 0x8, 0x2, 0x8e, 0x91, 0x5, 0x18, 
    0xd, 0x2, 0x8f, 0x91, 0x5, 0x14, 0xb, 0x2, 0x90, 0x8c, 0x3, 0x2, 0x2, 
    0x2, 0x90, 0x8d, 0x3, 0x2, 0x2, 0x2, 0x90, 0x8e, 0x3, 0x2, 0x2, 0x2, 
    0x90, 0x8f, 0x3, 0x2, 0x2, 0x2, 0x91, 0x17, 0x3, 0x2, 0x2, 0x2, 0x92, 
    0x93, 0x7, 0x18, 0x2, 0x2, 0x93, 0x94, 0x5, 0xc, 0x7, 0x2, 0x94, 0x96, 
    0x7, 0x19, 0x2, 0x2, 0x95, 0x97, 0x5, 0x16, 0xc, 0x2, 0x96, 0x95, 0x3, 
    0x2, 0x2, 0x2, 0x97, 0x98, 0x3, 0x2, 0x2, 0x2, 0x98, 0x96, 0x3, 0x2, 
    0x2, 0x2, 0x98, 0x99, 0x3, 0x2, 0x2, 0x2, 0x99, 0x9a, 0x3, 0x2, 0x2, 
    0x2, 0x9a, 0x9b, 0x7, 0x1a, 0x2, 0x2, 0x9b, 0x19, 0x3, 0x2, 0x2, 0x2, 
    0xe, 0x1d, 0x35, 0x3f, 0x4f, 0x5d, 0x5f, 0x6a, 0x77, 0x85, 0x88, 0x90, 
    0x98, 
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

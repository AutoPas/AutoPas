/**
 * @file TranslationVisitor.cpp
 * @author Manuel Lerchner
 * @date 09.05.24
 */

#include "TranslationVisitor.h"

#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzySetFactory.h"
#include "autopas/utils/ExceptionHandler.h"

antlrcpp::Any TranslationVisitor::visitRule_file(FuzzyLanguageParser::Rule_fileContext *context) {
  std::vector<std::shared_ptr<LinguisticVariable>> linguisticVariables;
  std::vector<FuzzyRule> fuzzyRules;

  for (auto *fuzzy_variableContext : context->linguistic_variable()) {
    std::shared_ptr<LinguisticVariable> lv = visit(fuzzy_variableContext).as<std::shared_ptr<LinguisticVariable>>();
    linguisticVariables.push_back(lv);

    _state._linguisticVariables[lv->getName()] = lv;
  }

  for (auto *fuzzy_ruleContext : context->fuzzy_rule()) {
    FuzzyRule rule = visit(fuzzy_ruleContext).as<FuzzyRule>();
    fuzzyRules.push_back(std::move(rule));
  }

  std::map<std::string, std::shared_ptr<FuzzyControlSystem>> fuzzyControlSystems;

  for (auto &rule : fuzzyRules) {
    auto dimensions = rule.getConsequent()->getCrispSet()->getDimensions();
    if (dimensions.size() != 1) {
      autopas::utils::ExceptionHandler::exception("Only rules with one dimensional output are supported! Rule: " +
                                                  std::string(rule));
    }
    std::string lvName = dimensions.begin()->first;

    if (fuzzyControlSystems.find(lvName) == fuzzyControlSystems.end()) {
      fuzzyControlSystems[lvName] = std::make_shared<FuzzyControlSystem>();
    }

    fuzzyControlSystems[lvName]->addRule(rule);
  }

  return std::make_pair(linguisticVariables, fuzzyControlSystems);
};

antlrcpp::Any TranslationVisitor::visitLinguistic_variable(FuzzyLanguageParser::Linguistic_variableContext *context) {
  // get the name and range of the linguistic variable
  std::string linguisticTerm = context->STRING()->getText();
  std::pair<double, double> range = {std::stod(context->NUMBER(0)->getText()),
                                     std::stod(context->NUMBER(1)->getText())};

  std::shared_ptr<LinguisticVariable> linguisticVariable = std::make_shared<LinguisticVariable>(linguisticTerm, range);

  // add all linguistic terms
  for (auto *fuzzy_termContext : context->fuzzy_term()) {
    std::shared_ptr<FuzzySet> term = visit(fuzzy_termContext).as<std::shared_ptr<FuzzySet>>();
    linguisticVariable->addLinguisticTerm(term);
  }

  return linguisticVariable;
};

antlrcpp::Any TranslationVisitor::visitFuzzy_term(FuzzyLanguageParser::Fuzzy_termContext *context) {
  std::string linguisticTerm = context->STRING()->getText();

  auto function = visit(context->function());
  auto [functionName, params] = function.as<std::pair<std::string, std::vector<double>>>();

  return FuzzySetFactory::makeFuzzySet(linguisticTerm, functionName, params);
};

antlrcpp::Any TranslationVisitor::visitFunction(FuzzyLanguageParser::FunctionContext *context) {
  std::string function = context->IDENTIFIER()->getText();
  std::vector<double> params;
  for (auto *number : context->NUMBER()) {
    params.push_back(std::stod(number->getText()));
  }

  return std::pair<std::string, std::vector<double>>(function, params);
};

antlrcpp::Any TranslationVisitor::visitFuzzy_rule(FuzzyLanguageParser::Fuzzy_ruleContext *context) {
  std::shared_ptr<FuzzySet> antecedent;
  std::shared_ptr<FuzzySet> consequent;

  antecedent = visit(context->fuzzy_set(0)).as<std::shared_ptr<FuzzySet>>();
  consequent = visit(context->fuzzy_set(1)).as<std::shared_ptr<FuzzySet>>();

  return FuzzyRule(antecedent, consequent);
};

antlrcpp::Any TranslationVisitor::visitOr(FuzzyLanguageParser::OrContext *context) {
  std::shared_ptr<FuzzySet> left = visit(context->fuzzy_set(0)).as<std::shared_ptr<FuzzySet>>();
  std::shared_ptr<FuzzySet> right = visit(context->fuzzy_set(1)).as<std::shared_ptr<FuzzySet>>();

  return left || right;
};

antlrcpp::Any TranslationVisitor::visitBrackets(FuzzyLanguageParser::BracketsContext *context) {
  return visit(context->fuzzy_set());
};

antlrcpp::Any TranslationVisitor::visitAnd(FuzzyLanguageParser::AndContext *context) {
  std::shared_ptr<FuzzySet> left = visit(context->fuzzy_set(0)).as<std::shared_ptr<FuzzySet>>();
  std::shared_ptr<FuzzySet> right = visit(context->fuzzy_set(1)).as<std::shared_ptr<FuzzySet>>();

  return left && right;
};

antlrcpp::Any TranslationVisitor::visitSelect(FuzzyLanguageParser::SelectContext *context) {
  std::string lvName = context->STRING(0)->getText();
  std::string termName = context->STRING(1)->getText();

  std::shared_ptr<LinguisticVariable> lv = _state._linguisticVariables[lvName];

  return lv->operator==(termName);
};

antlrcpp::Any TranslationVisitor::visitNegate(FuzzyLanguageParser::NegateContext *context) {
  std::shared_ptr<FuzzySet> fuzzySet = visit(context->fuzzy_set()).as<std::shared_ptr<FuzzySet>>();
  return !fuzzySet;
};

//
//
///**
// * Translates a program tree non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitProgram(RuleLanguageParser::ProgramContext *ctx) override {
//  std::vector<std::shared_ptr<Statement>> statements;
//  for (auto *statementContext : ctx->statement()) {
//    statements.push_back(getStatementType(visit(statementContext)));
//  }
//  return RuleBasedProgramTree{statements};
//}
//
///**
// * Translates a literal non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitLiteral(RuleLanguageParser::LiteralContext *ctx) override {
//  RuleVM::MemoryCell literal;
//  if (ctx->Bool_val()) {
//    literal = ctx->Bool_val()->getText() == "true";
//  } else if (ctx->Container_opt()) {
//    literal = ContainerOption::parseOptionExact(ctx->Container_opt()->getText());
//  } else if (ctx->Traversal_opt()) {
//    literal = TraversalOption::parseOptionExact(ctx->Traversal_opt()->getText());
//  } else if (ctx->Data_layout_opt()) {
//    literal = DataLayoutOption::parseOptionExact(ctx->Data_layout_opt()->getText());
//  } else if (ctx->Load_estimator_opt()) {
//    literal = LoadEstimatorOption::parseOptionExact(ctx->Load_estimator_opt()->getText());
//  } else if (ctx->Newton3_opt()) {
//    literal = Newton3Option::parseOptionExact(ctx->Newton3_opt()->getText());
//  } else if (ctx->unsigned_val()) {
//    literal = RuleVM::MemoryCell{static_cast<size_t>(std::stoull(ctx->unsigned_val()->getText()))};
//  } else if (ctx->Double_val()) {
//    literal = std::stod(ctx->Double_val()->getText());
//  } else {
//    throw std::runtime_error("literal '" + ctx->getText() + "' could not be parsed");
//  }
//
//  return std::make_shared<Literal>(literal);
//}
//
///**
// * Translates a define_list non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitDefine_list(RuleLanguageParser::Define_listContext *ctx) override {
//  std::vector<Literal> values;
//  for (auto *value : ctx->literal()) {
//    values.push_back(*(visit(value).as<std::shared_ptr<Literal>>().get()));
//  }
//  return std::make_shared<DefineList>(ctx->Variable_name()->getText(), values);
//}
//
///**
// * Translates a define non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitDefine(RuleLanguageParser::DefineContext *ctx) override {
//  return std::make_shared<Define>(ctx->Variable_name()->getText(), getExprType(visit(ctx->expression())));
//}
//
///**
// * Translates a variable non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitVariable(RuleLanguageParser::VariableContext *ctx) override {
//  auto varName = ctx->Variable_name()->getText();
//  auto it = parserContext.definitions.find(varName);
//  if (it != parserContext.definitions.end()) {
//    return std::make_shared<Variable>(it->second);
//  } else {
//    return std::make_shared<Variable>(context.definitionOf(varName));
//  }
//}
//
///**
// * Translates a expression non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitExpression(RuleLanguageParser::ExpressionContext *ctx) override {
//  if (ctx->expression().size() > 1) {
//    ;
//    static const std::map<std::string, BinaryOperator::Operator> opMap{
//        {"*", BinaryOperator::MUL},   {"/", BinaryOperator::DIV},     {"+", BinaryOperator::ADD},
//        {"-", BinaryOperator::SUB},   {">", BinaryOperator::GREATER}, {"<", BinaryOperator::LESS},
//        {"and", BinaryOperator::AND}, {"or", BinaryOperator::OR}};
//    return std::make_shared<BinaryOperator>(opMap.at(ctx->op->getText()), getExprType(visit(ctx->expression(0))),
//                                            getExprType(visit(ctx->expression(1))));
//  } else if (ctx->literal()) {
//    return visit(ctx->literal());
//  } else if (ctx->variable()) {
//    return visit(ctx->variable());
//  } else if (ctx->op->getText() == "not") {
//    return std::make_shared<UnaryOperator>(UnaryOperator::NOT, getExprType(visit(ctx->expression(0))));
//  } else {
//    return visit(ctx->expression(0));
//  }
//}
//
///**
// * Translates a property_value non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitProperty_value(RuleLanguageParser::Property_valueContext *ctx) override {
//  return visitChildren(ctx);
//}
//
///**
// * Returns the list with the given name in the already parsed program.
// * @param name The name of the list.
// * @return The list with the name if defined.
// */
//[[nodiscard]] const DefineList *resolveList(const std::string &name) const {
//  auto it = parserContext.lists.find(name);
//  if (it != parserContext.lists.end()) {
//    return it->second;
//  } else {
//    return context.getList(name);
//  }
//}
//
///**
// * Translates a configuration_pattern non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *ctx) override {
//  ConfigurationPattern pattern;
//  for (size_t i = 0; i < ctx->Configuration_property().size(); i++) {
//    auto property = ctx->Configuration_property(i)->getText();
//    auto *val = ctx->property_value(i);
//    std::vector<Literal> value;
//    if (val->literal()) {
//      value.push_back(*(visit(val->literal()).as<std::shared_ptr<Literal>>().get()));
//    } else {
//      value = resolveList(val->Variable_name()->getText())->values;
//    }
//
//    for (const auto &literal : value) {
//      pattern.add(literal.value);
//    }
//
//    // sanity check
//    const auto knownProperties = {
//        "container", "traversal", "dataLayout", "newton3", "loadEstimator", "cellSizeFactor",
//    };
//    if (std::find(knownProperties.begin(), knownProperties.end(), property) == knownProperties.end()) {
//      utils::ExceptionHandler::exception("RuleBasedProgramParser: Encountered unknown property! (" + property + ")");
//    }
//  }
//  return pattern;
//}
//
///**
// * Translates a configuration_order non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitConfiguration_order(RuleLanguageParser::Configuration_orderContext *ctx) override {
//  std::vector<ConfigurationOrder::SameProperty> sameProperties{};
//  for (size_t i = 0; i < ctx->Configuration_property().size(); i++) {
//    auto property = ctx->Configuration_property(i)->getText();
//    if (property == "container") {
//      sameProperties.push_back(ConfigurationOrder::SameProperty::container);
//    } else if (property == "traversal") {
//      sameProperties.push_back(ConfigurationOrder::SameProperty::traversal);
//    } else if (property == "dataLayout") {
//      sameProperties.push_back(ConfigurationOrder::SameProperty::dataLayout);
//    } else if (property == "newton3") {
//      sameProperties.push_back(ConfigurationOrder::SameProperty::newton3);
//    } else if (property == "loadEstimator") {
//      sameProperties.push_back(ConfigurationOrder::SameProperty::loadEstimator);
//    } else if (property == "cellSizeFactor") {
//      sameProperties.push_back(ConfigurationOrder::SameProperty::cellSizeFactor);
//    }
//  }
//  return std::make_shared<ConfigurationOrder>(visit(ctx->configuration_pattern(0)).as<ConfigurationPattern>(),
//                                              visit(ctx->configuration_pattern(1)).as<ConfigurationPattern>(),
//                                              sameProperties);
//}
//
///**
// * Translates a statement non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitStatement(RuleLanguageParser::StatementContext *ctx) override {
//  auto res = visitChildren(ctx);
//  if (res.is<std::shared_ptr<Define>>()) {
//    parserContext.definitions[res.as<std::shared_ptr<Define>>()->variable] = res.as<std::shared_ptr<Define>>().get();
//  } else if (res.is<std::shared_ptr<DefineList>>()) {
//    parserContext.lists[res.as<std::shared_ptr<DefineList>>()->listName] =
//    res.as<std::shared_ptr<DefineList>>().get();
//  }
//  return res;
//}
//
///**
// * Translates a if_statement non-literal.
// * @param ctx The parser context.
// * @return The antlrcpp::Any containing the parsed AST result.
// */
// antlrcpp::Any visitIf_statement(RuleLanguageParser::If_statementContext *ctx) override {
//  auto condition = getExprType(visit(ctx->expression()));
//  std::vector<std::shared_ptr<Statement>> statements;
//  for (auto *statementContext : ctx->statement()) {
//    statements.push_back(getStatementType(visit(statementContext)));
//  }
//  return std::make_shared<If>(condition, statements);
//}

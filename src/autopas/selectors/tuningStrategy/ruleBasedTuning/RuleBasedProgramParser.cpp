#include "RuleBasedProgramParser.h"

#include "antlr4-runtime.h"
#include "parser_generated/autopas_generated_rule_syntax/RuleLanguageBaseVisitor.h"
#include "parser_generated/autopas_generated_rule_syntax/RuleLanguageLexer.h"
#include "parser_generated/autopas_generated_rule_syntax/RuleLanguageParser.h"

namespace autopas::rule_syntax {
using namespace autopas_generated_rule_syntax;

/**
 * The antlr visitor that produces the ASt.
 */
class TranslationVisitor : public RuleLanguageBaseVisitor {
  struct ParserContext {
    std::map<std::string, const Define *> definitions;
    std::map<std::string, const DefineList *> lists;
  };

 public:
  explicit TranslationVisitor(CodeGenerationContext &context) : context(context) {}

  static std::shared_ptr<Expression> getExprType(const antlrcpp::Any &expr) {
    if (expr.is<std::shared_ptr<BinaryOperator>>()) {
      return expr.as<std::shared_ptr<BinaryOperator>>();
    }
    if (expr.is<std::shared_ptr<Variable>>()) {
      return expr.as<std::shared_ptr<Variable>>();
    }
    if (expr.is<std::shared_ptr<Literal>>()) {
      return expr.as<std::shared_ptr<Literal>>();
    }
    if (expr.is<std::shared_ptr<UnaryOperator>>()) {
      return expr.as<std::shared_ptr<UnaryOperator>>();
    }
    throw std::runtime_error("not an expression");
  }

  static std::shared_ptr<Statement> getStatementType(const antlrcpp::Any &statement) {
    if (statement.is<std::shared_ptr<DefineList>>()) {
      return statement.as<std::shared_ptr<DefineList>>();
    }
    if (statement.is<std::shared_ptr<Define>>()) {
      return statement.as<std::shared_ptr<Define>>();
    }
    if (statement.is<std::shared_ptr<If>>()) {
      return statement.as<std::shared_ptr<If>>();
    }
    if (statement.is<std::shared_ptr<ConfigurationOrder>>()) {
      return statement.as<std::shared_ptr<ConfigurationOrder>>();
    }
    throw std::runtime_error("not a statement");
  }

  antlrcpp::Any visitProgram(RuleLanguageParser::ProgramContext *ctx) override {
    std::vector<std::shared_ptr<Statement>> statements;
    for (auto *statementContext : ctx->statement()) {
      statements.push_back(getStatementType(visit(statementContext)));
    }
    return RuleBasedProgramTree{statements};
  }

  antlrcpp::Any visitLiteral(RuleLanguageParser::LiteralContext *ctx) override {
    RuleVM::MemoryCell literal;
    if (ctx->Bool_val()) {
      literal = ctx->Bool_val()->getText() == "true";
    } else if (ctx->Container_opt()) {
      literal = ContainerOption::parseOptionExact(ctx->Container_opt()->getText());
    } else if (ctx->Traversal_opt()) {
      literal = TraversalOption::parseOptionExact(ctx->Traversal_opt()->getText());
    } else if (ctx->Data_layout_opt()) {
      literal = DataLayoutOption::parseOptionExact(ctx->Data_layout_opt()->getText());
    } else if (ctx->Load_estimator_opt()) {
      literal = LoadEstimatorOption::parseOptionExact(ctx->Load_estimator_opt()->getText());
    } else if (ctx->Newton3_opt()) {
      literal = Newton3Option::parseOptionExact(ctx->Newton3_opt()->getText());
    } else if (ctx->unsigned_val()) {
      literal = RuleVM::MemoryCell{static_cast<size_t>(std::stoull(ctx->unsigned_val()->getText()))};
    } else if (ctx->Double_val()) {
      literal = std::stod(ctx->Double_val()->getText());
    } else {
      throw std::runtime_error("literal '" + ctx->getText() + "' could not be parsed");
    }

    return std::make_shared<Literal>(literal);
  }

  antlrcpp::Any visitDefine_list(RuleLanguageParser::Define_listContext *ctx) override {
    std::vector<Literal> values;
    for (auto *value : ctx->literal()) {
      values.push_back(*(visit(value).as<std::shared_ptr<Literal>>().get()));
    }
    return std::make_shared<DefineList>(ctx->Variable_name()->getText(), values);
  }

  antlrcpp::Any visitDefine(RuleLanguageParser::DefineContext *ctx) override {
    return std::make_shared<Define>(ctx->Variable_name()->getText(), getExprType(visit(ctx->expression())));
  }

  antlrcpp::Any visitVariable(RuleLanguageParser::VariableContext *ctx) override {
    auto varName = ctx->Variable_name()->getText();
    auto it = parserContext.definitions.find(varName);
    if (it != parserContext.definitions.end()) {
      return std::make_shared<Variable>(it->second);
    } else {
      return std::make_shared<Variable>(context.definitionOf(varName));
    }
  }

  antlrcpp::Any visitExpression(RuleLanguageParser::ExpressionContext *ctx) override {
    if (ctx->expression().size() > 1) {
      ;
      static const std::map<std::string, BinaryOperator::Operator> opMap{
          {"*", BinaryOperator::MUL},   {"/", BinaryOperator::DIV},     {"+", BinaryOperator::ADD},
          {"-", BinaryOperator::SUB},   {">", BinaryOperator::GREATER}, {"<", BinaryOperator::LESS},
          {"and", BinaryOperator::AND}, {"or", BinaryOperator::OR}};
      return std::make_shared<BinaryOperator>(opMap.at(ctx->op->getText()), getExprType(visit(ctx->expression(0))),
                                              getExprType(visit(ctx->expression(1))));
    } else if (ctx->literal()) {
      return visit(ctx->literal());
    } else if (ctx->variable()) {
      return visit(ctx->variable());
    } else if (ctx->op->getText() == "not") {
      return std::make_shared<UnaryOperator>(UnaryOperator::NOT, getExprType(visit(ctx->expression(0))));
    } else {
      return visit(ctx->expression(0));
    }
  }

  antlrcpp::Any visitProperty_value(RuleLanguageParser::Property_valueContext *ctx) override {
    return visitChildren(ctx);
  }

  [[nodiscard]] auto resolveList(const std::string &name) const {
    auto it = parserContext.lists.find(name);
    if (it != parserContext.lists.end()) {
      return it->second;
    } else {
      return context.getList(name);
    }
  }

  antlrcpp::Any visitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *ctx) override {
    ConfigurationPattern pattern;
    for (size_t i = 0; i < ctx->Configuration_property().size(); i++) {
      auto property = ctx->Configuration_property(i)->getText();
      auto val = ctx->property_value(i);
      std::vector<Literal> value;
      if (val->literal()) {
        value.push_back(*(visit(val->literal()).as<std::shared_ptr<Literal>>().get()));
      } else {
        value = resolveList(val->Variable_name()->getText())->values;
      }

      for (const auto &literal : value) {
        pattern.add(literal.value);
      }

      // TODO: Somehow check that type is correct. Maybe already in .g4 file during parsing
      if (property == "container") {
      } else if (property == "traversal") {
      } else if (property == "dataLayout") {
      } else if (property == "newton3") {
      } else if (property == "loadEstimator") {
      } else if (property == "cellSizeFactor") {
      }
    }
    return pattern;
  }

  antlrcpp::Any visitConfiguration_order(RuleLanguageParser::Configuration_orderContext *ctx) override {
    std::vector<ConfigurationOrder::SameProperty> sameProperties{};
    for (size_t i = 0; i < ctx->Configuration_property().size(); i++) {
      auto property = ctx->Configuration_property(i)->getText();
      if (property == "container") {
        sameProperties.push_back(ConfigurationOrder::SameProperty::container);
      } else if (property == "traversal") {
        sameProperties.push_back(ConfigurationOrder::SameProperty::traversal);
      } else if (property == "dataLayout") {
        sameProperties.push_back(ConfigurationOrder::SameProperty::dataLayout);
      } else if (property == "newton3") {
        sameProperties.push_back(ConfigurationOrder::SameProperty::newton3);
      } else if (property == "loadEstimator") {
        sameProperties.push_back(ConfigurationOrder::SameProperty::loadEstimator);
      } else if (property == "cellSizeFactor") {
        sameProperties.push_back(ConfigurationOrder::SameProperty::cellSizeFactor);
      }
    }
    return std::make_shared<ConfigurationOrder>(visit(ctx->configuration_pattern(0)).as<ConfigurationPattern>(),
                                                visit(ctx->configuration_pattern(1)).as<ConfigurationPattern>(),
                                                sameProperties);
  }

  antlrcpp::Any visitStatement(RuleLanguageParser::StatementContext *ctx) override {
    auto res = visitChildren(ctx);
    if (res.is<std::shared_ptr<Define>>()) {
      parserContext.definitions[res.as<std::shared_ptr<Define>>()->variable] = res.as<std::shared_ptr<Define>>().get();
    } else if (res.is<std::shared_ptr<DefineList>>()) {
      parserContext.lists[res.as<std::shared_ptr<DefineList>>()->listName] =
          res.as<std::shared_ptr<DefineList>>().get();
    }
    return res;
  }

  antlrcpp::Any visitIf_statement(RuleLanguageParser::If_statementContext *ctx) override {
    auto condition = getExprType(visit(ctx->expression()));
    std::vector<std::shared_ptr<Statement>> statements;
    for (auto *statementContext : ctx->statement()) {
      statements.push_back(getStatementType(visit(statementContext)));
    }
    return std::make_shared<If>(condition, statements);
  }

 private:
  CodeGenerationContext &context;
  ParserContext parserContext;
};

std::pair<RuleBasedProgramTree, CodeGenerationContext> RuleBasedProgramParser::parse(const std::string &programCode) {
  CodeGenerationContext context{{}};
  for (const auto &def : _initialDefinitions) {
    context.addGlobalVariable(def.second);
  }

  antlr4::ANTLRInputStream input{programCode};
  RuleLanguageLexer lexer(&input);
  antlr4::CommonTokenStream tokens(&lexer);

  tokens.fill();

  RuleLanguageParser parser(&tokens);
  antlr4::tree::ParseTree *tree = parser.program();

  TranslationVisitor visitor{context};
  auto ownTree = visitor.visit(tree).as<RuleBasedProgramTree>();

  return {ownTree, context};
}

}  // namespace autopas::rule_syntax

#include "RuleBasedProgramParser.h"

#include "antlr4-runtime.h"
#include "parser_generated/autopas_generated_rule_syntax/RuleLanguageBaseVisitor.h"
#include "parser_generated/autopas_generated_rule_syntax/RuleLanguageLexer.h"
#include "parser_generated/autopas_generated_rule_syntax/RuleLanguageParser.h"

namespace autopas::RuleSyntax {
using namespace AutopasGeneratedRuleSyntax;
using std::any_cast;

/**
 * The antlr visitor that produces the AST.
 */
class TranslationVisitor : public RuleLanguageBaseVisitor {
  struct ParserContext {
    /**
     * Constructor
     */
    ParserContext();

    /**
     * Destructor
     */
    ~ParserContext() noexcept;

    std::map<std::string, const Define *> definitions;
    std::map<std::string, const DefineList *> lists;
  };

 public:
  /**
   * Constructs the visitor with a given CodeGenerationContext.
   * @param context The CodeGenerationContext.
   */
  explicit TranslationVisitor(CodeGenerationContext &context) : context(context) {}

  /**
   * Helper function that returns a std::shared_ptr<Expression> from an std::any for multiple different expression
   * types.
   * @param expr The std::any.
   * @return The converted std::shared_ptr<Expression>.
   */
  static std::shared_ptr<Expression> getExprType(const std::any &expr) {
    if (expr.type() == typeid(std::shared_ptr<BinaryOperator>)) {
      return any_cast<std::shared_ptr<BinaryOperator>>(expr);
    }
    if (expr.type() == typeid(std::shared_ptr<Variable>)) {
      return any_cast<std::shared_ptr<Variable>>(expr);
    }
    if (expr.type() == typeid(std::shared_ptr<Literal>)) {
      return any_cast<std::shared_ptr<Literal>>(expr);
    }
    if (expr.type() == typeid(std::shared_ptr<UnaryOperator>)) {
      return any_cast<std::shared_ptr<UnaryOperator>>(expr);
    }
    throw std::runtime_error("not an expression");
  }

  /**
   * Helper function that returns a std::shared_ptr<Statement> from an std::any for multiple different statement
   * types.
   * @param statement The std::any.
   * @return The converted std::shared_ptr<Statement>.
   */
  static std::shared_ptr<Statement> getStatementType(const antlrcpp::Any &statement) {
    if (statement.type() == typeid(std::shared_ptr<DefineList>)) {
      return any_cast<std::shared_ptr<DefineList>>(statement);
    }
    if (statement.type() == typeid(std::shared_ptr<Define>)) {
      return any_cast<std::shared_ptr<Define>>(statement);
    }
    if (statement.type() == typeid(std::shared_ptr<If>)) {
      return any_cast<std::shared_ptr<If>>(statement);
    }
    if (statement.type() == typeid(std::shared_ptr<ConfigurationOrder>)) {
      return any_cast<std::shared_ptr<ConfigurationOrder>>(statement);
    }
    throw std::runtime_error("not a statement");
  }

  /**
   * Translates a program tree non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitProgram(RuleLanguageParser::ProgramContext *ctx) override {
    std::vector<std::shared_ptr<Statement>> statements;
    for (auto *statementContext : ctx->statement()) {
      statements.push_back(getStatementType(visit(statementContext)));
    }
    return RuleBasedProgramTree{statements};
  }

  /**
   * Translates a literal non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitLiteral(RuleLanguageParser::LiteralContext *ctx) override {
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

  /**
   * Translates a define_list non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitDefine_list(RuleLanguageParser::Define_listContext *ctx) override {
    std::vector<Literal> values;
    for (auto *value : ctx->literal()) {
      values.push_back(*any_cast<std::shared_ptr<Literal>>(visit(value)));
    }
    return std::make_shared<DefineList>(ctx->Variable_name()->getText(), values);
  }

  /**
   * Translates a define non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitDefine(RuleLanguageParser::DefineContext *ctx) override {
    return std::make_shared<Define>(ctx->Variable_name()->getText(), getExprType(visit(ctx->expression())));
  }

  /**
   * Translates a variable non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitVariable(RuleLanguageParser::VariableContext *ctx) override {
    const auto varName = ctx->Variable_name()->getText();
    const auto it = parserContext.definitions.find(varName);
    if (it != parserContext.definitions.end()) {
      return std::make_shared<Variable>(it->second);
    } else {
      return std::make_shared<Variable>(context.definitionOf(varName));
    }
  }

  /**
   * Translates a expression non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitExpression(RuleLanguageParser::ExpressionContext *ctx) override {
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

  /**
   * Translates a property_value non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitProperty_value(RuleLanguageParser::Property_valueContext *ctx) override { return visitChildren(ctx); }

  /**
   * Returns the list with the given name in the already parsed program.
   * @param name The name of the list.
   * @return The list with the name if defined.
   */
  [[nodiscard]] const DefineList *resolveList(const std::string &name) const {
    auto it = parserContext.lists.find(name);
    if (it != parserContext.lists.end()) {
      return it->second;
    } else {
      return context.getList(name);
    }
  }

  /**
   * Translates a configuration_pattern non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *ctx) override {
    ConfigurationPattern pattern;
    for (size_t i = 0; i < ctx->Configuration_property().size(); i++) {
      auto property = ctx->Configuration_property(i)->getText();
      auto *val = ctx->property_value(i);
      std::vector<Literal> value;
      if (val->literal()) {
        value.push_back(*any_cast<std::shared_ptr<Literal>>(visit(val->literal())));
      } else {
        value = resolveList(val->Variable_name()->getText())->values;
      }

      for (const auto &literal : value) {
        pattern.add(literal.value);
      }

      // sanity check
      const auto knownProperties = {
          "container", "traversal", "dataLayout", "newton3", "loadEstimator", "cellSizeFactor",
      };
      if (std::find(knownProperties.begin(), knownProperties.end(), property) == knownProperties.end()) {
        utils::ExceptionHandler::exception("RuleBasedProgramParser: Encountered unknown property! (" + property + ")");
      }
    }
    return pattern;
  }

  /**
   * Translates a configuration_order non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitConfiguration_order(RuleLanguageParser::Configuration_orderContext *ctx) override {
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
    return std::make_shared<ConfigurationOrder>(any_cast<ConfigurationPattern>(visit(ctx->configuration_pattern(0))),
                                                any_cast<ConfigurationPattern>(visit(ctx->configuration_pattern(1))),
                                                sameProperties);
  }

  /**
   * Translates a statement non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitStatement(RuleLanguageParser::StatementContext *ctx) override {
    auto res = visitChildren(ctx);
    if (res.type() == typeid(std::shared_ptr<Define>)) {
      parserContext.definitions[any_cast<std::shared_ptr<Define>>(res)->variable] =
          any_cast<std::shared_ptr<Define>>(res).get();
    } else if (res.type() == typeid(std::shared_ptr<DefineList>)) {
      parserContext.lists[any_cast<std::shared_ptr<DefineList>>(res)->listName] =
          any_cast<std::shared_ptr<DefineList>>(res).get();
    }
    return res;
  }

  /**
   * Translates a if_statement non-literal.
   * @param ctx The parser context.
   * @return The std::any containing the parsed AST result.
   */
  std::any visitIf_statement(RuleLanguageParser::If_statementContext *ctx) override {
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
TranslationVisitor::ParserContext::ParserContext() = default;

TranslationVisitor::ParserContext::~ParserContext() noexcept = default;

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
  auto ownTree = any_cast<RuleBasedProgramTree>(visitor.visit(tree));

  return {ownTree, context};
}

}  // namespace autopas::RuleSyntax

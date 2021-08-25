#include "RuleBasedProgramParser.h"

#include "antlr4-runtime.h"

#include "parser_generated/autopas_rule_syntax/RuleLanguageLexer.h"
#include "parser_generated/autopas_rule_syntax/RuleLanguageParser.h"
#include "parser_generated/autopas_rule_syntax/RuleLanguageBaseVisitor.h"

namespace autopas::rule_syntax {

std::pair<RuleBasedProgramTree, CodeGenerationContext> RuleBasedProgramParser::parse(const std::string &programCode) {
  CodeGenerationContext context{{}};
  for (const auto &def : _initialDefinitions) {
    context.addVariable(def.second);
  }

  RuleBasedProgramTree program;
  std::stringstream input{programCode};
  while (not input.eof()) {
    // program.statements.push_back(parseStatement(input));
  }

  return {program, context};
}

class TranslationVisitor : public RuleLanguageBaseVisitor {
  struct ParserContext {
    std::map<std::string, const Define*> definitions;
    std::map<std::string, const DefineList*> lists;
  };

 public:
  explicit TranslationVisitor(CodeGenerationContext& context) : context(context) {}

  antlrcpp::Any visitProgram(RuleLanguageParser::ProgramContext *ctx) override {
    std::vector<StatementVal> statements;
    for(auto* statementContext : ctx->statement()) {
      statements.push_back(StatementVal{visit(statementContext)});
    }
    return RuleBasedProgramTree{statements};
  }

  antlrcpp::Any visitLiteral(RuleLanguageParser::LiteralContext *ctx) override {
    RuleVM::MemoryCell literal;
    if(ctx->Bool_val()) {
      literal = ctx->Bool_val()->getText() == "true";
    }
    else if(ctx->Container_opt()) {
      literal = ContainerOption::parseOptionExact(ctx->Container_opt()->getText());
    }
    else if(ctx->Traversal_opt()) {
      literal = TraversalOption::parseOptionExact(ctx->Traversal_opt()->getText());
    }
    else if(ctx->Data_layout_opt()) {
      literal = DataLayoutOption::parseOptionExact(ctx->Data_layout_opt()->getText());
    }
    else if(ctx->Load_estimator_opt()) {
      literal = LoadEstimatorOption::parseOptionExact(ctx->Load_estimator_opt()->getText());
    }
    else if(ctx->Newton3_opt()) {
      literal = Newton3Option::parseOptionExact(ctx->Newton3_opt()->getText());
    }
    else if(ctx->Unsigned_val()) {
      literal = std::stoull(ctx->Unsigned_val()->getText());
    }
    else {
      throw std::runtime_error("literal could not be parsed");
    }

    return Literal{literal};
  }

  antlrcpp::Any visitDefine_list(RuleLanguageParser::Define_listContext *ctx) override {
    std::vector<Literal> values;
    for(auto* value : ctx->literal()) {
      values.push_back(visit(value));
    }
    return DefineList{ctx->Variable_name()->getText(), values};
  }

  antlrcpp::Any visitDefine(RuleLanguageParser::DefineContext *ctx) override {
    return Define{ctx->Variable_name()->getText(), visit(ctx->literal()).as<Literal>()};
  }

  antlrcpp::Any visitVariable(RuleLanguageParser::VariableContext *ctx) override {
    auto varName = ctx->Variable_name()->getText();
    auto it = parserContext.definitions.find(varName);
    if(it != parserContext.definitions.end()) {
      return Variable{it->second};
    } else {
      return Variable{context.definitionOf(varName)};
    }
  }

  antlrcpp::Any visitExpression(RuleLanguageParser::ExpressionContext *ctx) override {
    if(ctx->literal()) {
      return visit(ctx->literal());
    }
    else if(ctx->variable()) {
      return visit(ctx->variable());
    }
    else {
      auto opText = ctx->Binary_op()->getText();
      auto op = opText == "<" ? BinaryOperator::LESS :
                              opText == ">" ? BinaryOperator::GREATER :
                              opText == "and" ? BinaryOperator::AND : BinaryOperator::LESS;
      return BinaryOperator{op, visit(ctx->expression(0)), visit(ctx->expression(1))};
    }
  }

  antlrcpp::Any visitProperty_value(RuleLanguageParser::Property_valueContext *ctx) override {
    return visitChildren(ctx);
  }

  [[nodiscard]] auto resolveList(const std::string& name) const {
    auto it = parserContext.lists.find(name);
    if (it != parserContext.lists.end()) {
      return it->second;
    } else {
      return context.getList(name);
    }
  }

  antlrcpp::Any visitConfiguration_pattern(RuleLanguageParser::Configuration_patternContext *ctx) override {
    ConfigurationPattern pattern;
    for(size_t i = 0; i < ctx->Configuration_property().size(); i++) {
      auto property = ctx->Configuration_property(i)->getText();
      auto val = ctx->property_value(i);
      std::vector<Literal> value;
      if(val->literal()) {
        value.push_back(visit(val->literal()));
      }
      else {
        value = resolveList(val->Variable_name()->getText())->values;
      }

      for(const auto& literal : value) {
        pattern.add(literal.value);
      }

      // TODO: Somehow check that type is correct. Maybe already in .g4 file during parsing
      if(property == "container") {

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
    return ConfigurationOrder(visit(ctx->configuration_pattern(0)), visit(ctx->configuration_pattern(1)));
  }

  antlrcpp::Any visitStatement(RuleLanguageParser::StatementContext *ctx) override {
    auto res = visitChildren(ctx);
    if(res.is<Define>()) {
      parserContext.definitions[res.as<Define>().variable] = &(res.as<Define>());
    }
    else if(res.is<DefineList>()) {
      parserContext.lists[res.as<DefineList>().listName] =  &(res.as<DefineList>());
    }
    return res;
  }

  antlrcpp::Any visitIf_statement(RuleLanguageParser::If_statementContext *ctx) override {
    auto condition = visit(ctx->expression());
    std::vector<StatementVal> statements;
    for(auto* statementContext : ctx->statement()) {
      statements.push_back(visit(statementContext));
    }
    return If{condition, statements};
  }

 private:
  CodeGenerationContext& context;
  ParserContext parserContext;
};

void RuleBasedProgramParser::test() {
  std::ifstream ifs("/home/tobias/tmp/myfile.txt");
  std::string content( (std::istreambuf_iterator<char>(ifs) ),
                       (std::istreambuf_iterator<char>()) );

  antlrcpp::Any any;
  any = 10;
  int var;
  var = decltype(var){any};

  antlr4::ANTLRInputStream input{content};
  RuleLanguageLexer lexer(&input);
  antlr4::CommonTokenStream tokens(&lexer);

  tokens.fill();
  for (auto token : tokens.getTokens()) {
    std::cout << token->toString() << std::endl;
  }

  RuleLanguageParser parser(&tokens);
  antlr4::tree::ParseTree* tree = parser.program();

  CodeGenerationContext context{{}};

  TranslationVisitor visitor{context};
  auto ownTree = visitor.visit(tree).as<RuleBasedProgramTree>();

  std::cout << tree->toStringTree(&parser) << std::endl << std::endl;
}

} // namespace autopas_rule_syntax
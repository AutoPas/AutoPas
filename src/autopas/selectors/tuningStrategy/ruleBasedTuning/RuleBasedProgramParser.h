#pragma once

#include <any>
#include <boost/spirit/home/x3.hpp>

#include "RuleBasedProgramTree.h"

namespace autopas::rule_syntax {



namespace grammar {
struct VariableDefContext {
  std::map<std::string, Define*> definitions;
};

using namespace boost::spirit;

struct ContainerOption_ : x3::symbols<Literal> {
  ContainerOption_() {
    for (const auto &[option, name] : ContainerOption::getOptionNames()) {
      add(name, Literal{option});
    }
  }
};

struct TraversalOption_ : x3::symbols<Literal> {
  TraversalOption_() {
    for (const auto& [option, name] : TraversalOption::getOptionNames()) {
      add(name, Literal{option});
    }
  }
};

struct LoadEstimatorOption_ : x3::symbols<Literal> {
  LoadEstimatorOption_() {
    for (const auto& [option, name] : LoadEstimatorOption::getOptionNames()) {
      add(name, Literal{option});
    }
  }
};

struct DataLayoutOption_ : x3::symbols<Literal> {
  DataLayoutOption_() {
    for (const auto& [option, name] : DataLayoutOption::getOptionNames()) {
      add(name, Literal{option});
    }
  }
};

struct Newton3Option_ : x3::symbols<Literal> {
  Newton3Option_() {
    for(const auto& [option, name] : Newton3Option::getOptionNames()) {
      add(name, Literal{option});
    }
  }
};

struct Bool_ : x3::symbols<Literal> {
  Bool_() {
    add("false", Literal{false}) ("true", Literal{true});
  }
};

struct BinaryOperator_ : x3::symbols<BinaryOperator::Operator> {
  BinaryOperator_() {
    add ("<", BinaryOperator::LESS) (">", BinaryOperator::GREATER) ("and", BinaryOperator::AND);
  }
};

struct ConfigurationProperty_ : x3::symbols<Type> {
  ConfigurationProperty_() {
    add ("container", Type::CONTAINER)
        ("traversal", Type::TRAVERSAL)
        ("dataLayout", Type::DATA_LAYOUT)
        ("newton3", Type::NEWTON3)
        ("loadEstimator", Type::LOAD_ESTIMATOR)
        ("cellSizeFactor", Type::CELL_SIZE_FACTOR)
      ;
  }
};

inline auto toLiteral = [](const auto& ctx) {
  x3::_val(ctx) = Literal{x3::_attr(ctx)};
};

/*inline auto toSharedPtr = [](const auto& ctx) {
  using value_t = std::remove_reference_t<decltype(x3::_attr(ctx))>;
  x3::_val(ctx) = std::make_shared<value_t>(x3::_attr(ctx));
};*/

inline auto variableDefineForName = [](const auto& ctx) {
  x3::_val(ctx) = x3::get<VariableDefContext>(ctx).get().definitions.at(x3::_attr(ctx));
};

inline auto addToVariableDefines = [](const auto& ctx) {
  x3::get<VariableDefContext>(ctx).get().definitions[x3::_attr(ctx)->variable] = &x3::_attr(ctx).get();
};

inline x3::rule <class unsigned_val, Literal> unsigned_val = "unsigned val";
auto const unsigned_val_def = x3::ulong_[toLiteral];

inline x3::rule<class literal_val, Literal> literal_val = "literal_val";
auto const literal_val_def = (TraversalOption_{} | ContainerOption_{} | LoadEstimatorOption_{} | DataLayoutOption_{}
                          | Newton3Option_{} | unsigned_val | Bool_{});

inline x3::rule<class literal, Literal> literal = "literal";
auto const literal_def = literal_val;

inline x3::rule<class variable_name, std::string> variable_name = "variable name";
auto const variable_name_def = x3::ascii::alpha >> *(x3::ascii::alnum);

inline x3::rule<class define_list, DefineList> define_list = "define_list";
auto const define_list_def = "define_list" >> variable_name >> '=' >> (literal % ',') >> ';';

inline x3::rule<class define, Define> define = "define";
auto const define_def = "define" >> variable_name >> '=' >> literal >> ';';

inline x3::rule<class variable, Variable> variable = "variable";
auto const variable_def = variable_name[variableDefineForName];

inline x3::rule<class expression, ExpressionVal> expression = "expression";
inline x3::rule<class binary_operator, BinaryOperator> binary_operator = "binary operator";

inline x3::rule<class property_value, std::vector<Literal>> property_value = "property value";
auto const property_value_def = variable_name | literal;

auto const expression_def = binary_operator | variable | literal;
auto const binary_operator_def = expression >> BinaryOperator_{} >> expression;

inline x3::rule<class configuration_pattern, ConfigurationPattern> configuration_pattern = "configuration pattern";
auto const configuration_pattern_def = '[' >> (ConfigurationProperty_{} >> '=' >> property_value) % ',' >> ']' ;

inline x3::rule<class configuration_order, ConfigurationOrder> configuration_order = "configuration order";
auto const configuration_order_def = configuration_pattern >> ">=" >> configuration_pattern >> ';';

inline x3::rule<class statement, StatementVal> statement = "statement";
inline x3::rule<class if_statement, If> if_statement = "if";

auto const if_statement_def = "if" >> expression >> ":" >> statement % x3::space >> "endif";
auto const statement_def = define_list | define[addToVariableDefines] | if_statement | configuration_order;

inline x3::rule<class program, RuleBasedProgramTree> program = "program";
auto const program_def = statement % x3::space;

BOOST_SPIRIT_DEFINE(unsigned_val, literal_val, literal, variable_name, define, define_list, variable, expression,
                    binary_operator, property_value, configuration_pattern, configuration_order,
                    statement, if_statement, program);
}

class RuleBasedProgramParser {
 public:
  explicit RuleBasedProgramParser(std::map<std::string, Define> initialDefinitions)
      : _initialDefinitions(std::move(initialDefinitions)) {}

  std::shared_ptr<Statement> parseStatement(std::stringstream &input) {
    switch (input.peek()) { }
    return {};
  }

  static void test() {
    using namespace boost::spirit;


    RuleBasedProgramTree program;
    std::any testRes;
    Define defi{"", Literal()};

    std::string test{"lc_c08"};
    grammar::VariableDefContext varContext;
    auto first = test.begin();
    bool matches = x3::phrase_parse(first, test.end(),
                                    x3::with<grammar::VariableDefContext>(std::ref(varContext))[grammar::expression],
                                    x3::space, testRes);
    std::cout << "Matches: " << matches << std::endl;
    std::cout << "Full Match: " << (first == test.end()) << std::endl;
    //for(const auto& val : values) {
      //std::cout << ContainerOption::getOptionNames()[val] << std::endl;
    //}
  }

  std::pair<RuleBasedProgramTree, CodeGenerationContext> parse(const std::string &programCode) {
    CodeGenerationContext context{{}};
    for (const auto &def : _initialDefinitions) {
      context.addVariable(def.second.variable);
    }

    RuleBasedProgramTree program;
    std::stringstream input{programCode};
    while (not input.eof()) {
      //program.statements.push_back(parseStatement(input));
    }

    return {{}, context};
  }

 private:
  std::map<std::string, Define> _initialDefinitions;
  std::map<std::string, DefineList> _lists;
};

} // namespace autopas
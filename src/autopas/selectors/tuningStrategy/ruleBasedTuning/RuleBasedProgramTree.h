#pragma once

#include "RuleVM.h"

namespace autopas {

/**
 * Enum describing all types that are allowed in a rule program.
 */
enum class Type { BOOL, DOUBLE, SIZE_T, CONTAINER, TRAVERSAL, LOAD_ESTIMATOR, DATA_LAYOUT, NEWTON3, CELL_SIZE_FACTOR };

/**
 * A configuration pattern that matches zero or more actual configurations, used in the rule language. A configuration
 * matches a pattern if all its components are allowed in the pattern. (If a component is empty in the pattern, it means
 * that all possible values are allowed).
 */
struct ConfigurationPattern {
  /**
   * The allowed containers.
   */
  std::set<ContainerOption> _containers;
  /**
   * The allowed traversals.
   */
  std::set<TraversalOption> _traversals;
  /**
   * The allowed load estimators.
   */
  std::set<LoadEstimatorOption> _loadEstimators;
  /**
   * The allowed data layouts.
   */
  std::set<DataLayoutOption> _dataLayouts;
  /**
   * The allowed newton3 options.
   */
  std::set<Newton3Option> _newton3Options;
  /**
   * The allowed cell size factors.
   */
  std::set<double> _cellSizeFactors;

  /**
   * Constructs a pattern that matches all configurations.
   */
  ConfigurationPattern() = default;

  /**
   * Adds a value to the allowed values of this configuration pattern. The component that is modified depends on the
   * contained value of the passed variant. (E.g. if a Newton3Option is passed in the variant, _newton3Options is
   * modified.
   * @param value The value to add to the allowed values of the configuration pattern (a variant).
   */
  void add(const RuleVM::MemoryCell &value) {
    std::visit(
        [&](const auto &val) {
          addHelper(val, _containers, _traversals, _loadEstimators, _dataLayouts, _newton3Options, _cellSizeFactors);
        },
        value);
  }

  /**
   * @param configuration the configuration to check
   * @return If the configuration matches this pattern.
   */
  [[nodiscard]] bool matches(const Configuration &configuration) const {
    return (_containers.empty() or contains(_containers, configuration.container)) and
           (_traversals.empty() or contains(_traversals, configuration.traversal)) and
           (_loadEstimators.empty() or contains(_loadEstimators, configuration.loadEstimator)) and
           (_dataLayouts.empty() or contains(_dataLayouts, configuration.dataLayout)) and
           (_newton3Options.empty() or contains(_newton3Options, configuration.newton3)) and
           (_cellSizeFactors.empty() or contains(_cellSizeFactors, configuration.cellSizeFactor));
  }

  /**
   * @returns A string that represents this configuration pattern.
   */
  [[nodiscard]] std::string toString() const {
    std::string res = optionSetToString(_containers) + "; " + optionSetToString(_traversals) + "; " +
                      optionSetToString(_loadEstimators) + "; " + optionSetToString(_loadEstimators) + "; " +
                      optionSetToString(_dataLayouts) + "; " + optionSetToString(_newton3Options) + "; ";
    if (not _cellSizeFactors.empty()) {
      res += std::accumulate(std::next(_cellSizeFactors.begin()), _cellSizeFactors.end(),
                             std::to_string(*_cellSizeFactors.begin()),
                             [](std::string s, double d) { return std::move(s) + ',', std::to_string(d); });
    }
    return res;
  }

 private:
  template <class T>
  [[nodiscard]] static std::string optionSetToString(const std::set<T> &set) {
    if (not set.empty()) {
      auto comma_fold = [](std::string a, T b) { return std::move(a) + ',' + b.to_string(); };

      return std::accumulate(std::next(set.begin()), set.end(), set.begin()->to_string(), comma_fold);
    }
    return {};
  }

  template <typename T>
  [[nodiscard]] static bool contains(const std::set<T> &set, const T &option) {
    return set.find(option) != set.end();
  }

  template <class T, class SetT>
  static void addHelper2(T value, SetT &set) {
    if constexpr (std::is_same_v<T, typename SetT::value_type>) {
      set.insert(value);
    }
  }

  template <class T, class... SetT>
  static void addHelper(T value, SetT &... sets) {
    (addHelper2(value, sets), ...);
  }
};

/**
 * Contains the AST of a rule program.
 */
namespace rule_syntax {

class CodeGenerationContext;

struct Define;

/**
 * A statement in the rule language.
 */
struct Statement {
  /**
   * Virtual default constructor.
   */
  virtual ~Statement() = default;
  /**
   * Generates code for this statement.
   * @param context The context to generate code in.
   * @param program The program to append the code to.
   */
  virtual void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const = 0;
};

/**
 * Converts the type of a RuleVM memory cell to the type in the AST. Depends on MemoryCell and Type defining the types
 * in the same order.
 * @param memoryCell The memory cell to get the type of.
 * @return The type of the passed MemoryCell in the AST.
 */
inline Type typeOf(const RuleVM::MemoryCell &memoryCell) { return static_cast<Type>(memoryCell.index()); }

/**
 * An expression in the rule language.
 */
struct Expression {
  /**
   * Virtual default constructor.
   */
  virtual ~Expression() = default;
  /**
   * @returns the type of the expression.
   */
  [[nodiscard]] virtual Type getType() const = 0;
  /**
   * Generates code for evaluating this expression at runtime. Value is put on top of the stack. May consume part of the
   * stack.
   * @param context The context to generate the code in.
   * @param program The program to append the code to.
   */
  virtual void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const = 0;
};

struct DefineList;

/**
 * A ConfigurationOrder statement in the rule language e.g.
 *     '[container=LinkedCells, newton3=enabled] >= [container=DirectSum] with same cellSizeFactor'
 * This order says that all configurations that match the left pattern (configuration has LinkedCells as container and
 * newton3 enabled) are better that all configurations that match the right pattern (configuration has container
 * DirectSum), provided they have the the same cellSizeFactor.
 */
struct ConfigurationOrder : public Statement {
  /**
   * The greater/better configuration pattern.
   */
  ConfigurationPattern greater;
  /**
   * The smaller/worse configuration pattern.
   */
  ConfigurationPattern smaller;

  /**
   * Enum containing all possible values of 'with same' options.
   */
  enum class SameProperty { container, traversal, dataLayout, newton3, loadEstimator, cellSizeFactor };

  /**
   * Vector containing all set same properties.
   */
  std::vector<SameProperty> sameProperties;

  /**
   * Constructs a configuration order that says all configurations are better than all configurations. Will never be
   * true.
   */
  ConfigurationOrder() = default;

  /**
   * Construct a configuration order.
   * @param greater The pattern of better configuration orders.
   * @param smaller The pattern of worse configuration orders.
   * @param same The same properties of the order.
   */
  ConfigurationOrder(ConfigurationPattern greater, ConfigurationPattern smaller, std::vector<SameProperty> same = {})
      : greater(std::move(greater)), smaller(std::move(smaller)), sameProperties(std::move(same)) {}

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override;

  /**
   * @returns A string describing this order.
   */
  [[nodiscard]] std::string toString() const {
    std::string sameString = sameProperties.empty() ? "" : "with same ";
    for (auto same : sameProperties) {
      sameString += samePropertyToString(same) + ", ";
    }
    return greater.toString() + " >= " + smaller.toString() + sameString;
  }

  /**
   * Checks if two configurations have the equal same properties as required by this configuration order.
   * @param conf1 One configuration.
   * @param conf2 Another configuration.
   * @return True, if both configurations have equal required same properties.
   */
  [[nodiscard]] bool haveEqualSameProperties(const Configuration &conf1, const Configuration &conf2) const {
    for (auto same : sameProperties) {
      switch (same) {
        case SameProperty::container:
          if (conf1.container != conf2.container) return false;
          break;
        case SameProperty::traversal:
          if (conf1.traversal != conf2.traversal) return false;
          break;
        case SameProperty::dataLayout:
          if (conf1.dataLayout != conf2.dataLayout) return false;
          break;
        case SameProperty::newton3:
          if (conf1.newton3 != conf2.newton3) return false;
          break;
        case SameProperty::loadEstimator:
          if (conf1.loadEstimator != conf2.loadEstimator) return false;
          break;
        case SameProperty::cellSizeFactor:
          if (conf1.cellSizeFactor != conf2.cellSizeFactor) return false;
          break;
      }
    }
    return true;
  }

 private:
  static std::string samePropertyToString(SameProperty same) {
    switch (same) {
      case SameProperty::container:
        return "container";
      case SameProperty::traversal:
        return "traversal";
      case SameProperty::dataLayout:
        return "dataLayout";
      case SameProperty::newton3:
        return "newton3";
      case SameProperty::loadEstimator:
        return "loadEstimator";
      case SameProperty::cellSizeFactor:
        return "cellSizeFactor";
    }
    return {};
  }
};

/**
 * A CodeGenerationContext keeps track of the needed stack size of the program, the mapping from variables to stack
 * ddresses, the currently known lists, and the configuration patterns with they indices in the program encountered so
 * far.
 */
class CodeGenerationContext {
 public:
  /**
   * Constructs a CodeGenerationContext.
   * @param initialAddressEnvironment The predefined variables with their definition and stack address.
   */
  explicit CodeGenerationContext(std::map<std::string, std::pair<const Define *, size_t>> initialAddressEnvironment)
      : currentNeededStack(0),
        maxNeededStack(0),
        addressEnvironment(std::move(initialAddressEnvironment)),
        initialNumVariables(addressEnvironment.size()) {}

  /**
   * Allocates the given number of MemoryCells on the stack.
   * @param num The number of cells to allocate.
   */
  void allocateStack(size_t num) {
    currentNeededStack += num;
    maxNeededStack = std::max(maxNeededStack, currentNeededStack);
  }

  /**
   * Frees the given number of MemoryCells on the stack.
   * @param num The number of cells to free.
   */
  void freeStack(size_t num) { currentNeededStack -= num; }

  /**
   * Adds a local variable definition that can be used in the program later.
   * @param definition The definition to add.
   */
  void addLocalVariable(const Define &definition);

  /**
   * Adds a global variable definition that can be used in the program later.
   * @param definition The definition to add.
   */
  void addGlobalVariable(const Define &definition);

  /**
   * @param name The name of a variable.
   * @returns the address of the variable on the stack.
   */
  [[nodiscard]] size_t addressOf(const std::string &name) const;

  /**
   * @param name The name of a variable.
   * @return The definition of the variable.
   */
  [[nodiscard]] const Define *definitionOf(const std::string &name) const;

  /**
   * @returns the maximum needed stack size of the program.
   */
  [[nodiscard]] auto getMaxStackSize() const { return maxNeededStack; }

  /**
   * Adds a list definition that can be used later in the program.
   * @param name The name of the list.
   * @param list The list definition.
   */
  void addList(const std::string &name, const DefineList *list) { lists.insert({name, list}); }

  /**
   * @param name The name of a defined list.
   * @returns the definition of the list with the given name.
   */
  auto getList(const std::string &name) { return lists.at(name); }

  /**
   * Adds a configuration order and assigns it an index that can be used to refer to it later.
   * @param configurationOrder The new configuration order.
   * @return The index of the new configuration order.
   */
  auto addConfigurationOrder(const ConfigurationOrder &configurationOrder) {
    configurationOrders.emplace_back(configurationOrder);
    return configurationOrders.size() - 1;
  }

  /**
   * @returns All configuration orders.
   */
  [[nodiscard]] const auto &getConfigurationOrders() const { return configurationOrders; }

  /**
   * @param idx An index of a configuration order.
   * @returns The smaller pattern of the configuration order with the given index.
   */
  [[nodiscard]] auto smallerConfigurationPatternByIndex(size_t idx) const {
    return configurationOrders.at(idx).smaller;
  }

  /**
   * @returns The number of local variables.
   */
  [[nodiscard]] auto getNumLocalVariables() const { return addressEnvironment.size() - initialNumVariables; }

 private:
  size_t currentNeededStack;
  size_t maxNeededStack;

  std::map<std::string, std::pair<const Define *, size_t>> addressEnvironment;
  size_t initialNumVariables;
  std::map<std::string, const DefineList *> lists;
  std::vector<ConfigurationOrder> configurationOrders;
};

/**
 * A literal in the rule language.
 */
struct Literal : public Expression {
  /**
   * The value of the literal.
   */
  RuleVM::MemoryCell value;
  /**
   * The type of the literal in AST.
   */
  Type type;

  /**
   * Constructs a literal from a MemoryCell.
   * @param value The memory cell to construct from.
   */
  explicit Literal(RuleVM::MemoryCell value) : value(value), type(typeOf(value)) {}

  /**
   * @returns the type of the literal.
   */
  [[nodiscard]] Type getType() const override { return type; }

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override {
    program.instructions.emplace_back(RuleVM::LOADC, value);
    context.allocateStack(1);
  }
};

/**
 * A list definition in the rule language.
 */
struct DefineList : public Statement {
  /**
   * The name of the list.
   */
  std::string listName;
  /**
   * The values in the list.
   */
  std::vector<Literal> values;

  /**
   * Constructs a list definition.
   * @param listName The list name.
   * @param values The values in the list.
   */
  DefineList(std::string listName, std::vector<Literal> values)
      : listName(std::move(listName)), values(std::move(values)) {}

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override {
    context.addList(listName, this);
  }
};

struct Variable : public Expression {
  const Define *definition;

  explicit Variable(const Define *definition) : definition(definition) {}

  [[nodiscard]] Type getType() const override;

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override;
};

struct UnaryOperator : public Expression {
  enum Operator { NOT };

  std::shared_ptr<Expression> child;
  Operator op;

  UnaryOperator(Operator op, std::shared_ptr<Expression> child)
    : child(std::move(child)), op(op) {}

  [[nodiscard]] Type getType() const override {
    return Type::BOOL;
  }

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override;
};

struct BinaryOperator : public Expression {
  enum Operator { LESS, GREATER, AND, OR, ADD, SUB, MUL, DIV };

  std::shared_ptr<Expression> left;
  Operator op;
  std::shared_ptr<Expression> right;

  BinaryOperator(Operator op, std::shared_ptr<Expression> left, std::shared_ptr<Expression> right)
      : left(std::move(left)), op(op), right(std::move(right)) {}

  [[nodiscard]] Type getType() const override {
    return op == LESS or op == GREATER or op == AND or op == OR
               ? Type::BOOL
               : (left->getType() == Type::SIZE_T and right->getType() == Type::SIZE_T ? Type::SIZE_T : Type::DOUBLE);
  }

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override;
};

struct Define : public Statement {
  std::string variable;
  std::shared_ptr<Expression> value;

  Define(std::string variable, std::shared_ptr<Expression> value)
      : variable(std::move(variable)), value(std::move(value)) {}

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override {
    value->generateCode(context, program);

    context.addLocalVariable(*this);
    program.instructions.emplace_back(RuleVM::STOREA, context.addressOf(variable));
    context.freeStack(1);
  }
};

struct If : public Statement {
  std::shared_ptr<Expression> condition;
  std::vector<std::shared_ptr<Statement>> consequences;

  If(std::shared_ptr<Expression> condition, std::vector<std::shared_ptr<Statement>> consequences)
      : condition(std::move(condition)), consequences(std::move(consequences)) {}

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override {
    condition->generateCode(context, program);
    program.instructions.emplace_back(RuleVM::JUMPZERO);
    context.freeStack(1);
    auto jumpInstrIdx = program.instructions.size() - 1;

    for (const auto &statement : consequences) {
      statement->generateCode(context, program);
    }

    program.instructions[jumpInstrIdx].payload = program.instructions.size();
  }
};

/**
 * The AST of a rule program.
 */
struct RuleBasedProgramTree {
  /**
   * All statements in the program.
   */
  std::vector<std::shared_ptr<Statement>> statements;

  /**
   * Generates code for the whole program described by this AST.
   * @param context The context to generate in the code in.
   * @returns The program for this AST.
   */
  [[nodiscard]] RuleVM::Program generateCode(CodeGenerationContext &context) const {
    RuleVM::Program program;
    program.instructions.emplace_back(RuleVM::RESERVE, 0ul);
    auto reserveInstructionIdx = program.instructions.size() - 1;
    for (const auto &statement : statements) {
      statement->generateCode(context, program);
    }
    program.instructions[reserveInstructionIdx].payload = context.getNumLocalVariables();
    program.neededStackSize = context.getNumLocalVariables() + context.getMaxStackSize();
    program.instructions.emplace_back(RuleVM::HALT);
    return program;
  }
};

}  // namespace rule_syntax
}  // namespace autopas

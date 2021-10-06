#pragma once

#include "RuleVM.h"

namespace autopas {

enum class Type { BOOL, DOUBLE, SIZE_T, CONTAINER, TRAVERSAL, LOAD_ESTIMATOR, DATA_LAYOUT, NEWTON3, CELL_SIZE_FACTOR };

struct ConfigurationPattern {
  std::set<ContainerOption> _containers;
  std::set<TraversalOption> _traversals;
  std::set<LoadEstimatorOption> _loadEstimators;
  std::set<DataLayoutOption> _dataLayouts;
  std::set<Newton3Option> _newton3Options;
  std::set<double> _cellSizeFactors;

  ConfigurationPattern() = default;

  void add(const RuleVM::MemoryCell &value) {
    std::visit(
        [&](const auto &val) {
          addHelper(val, _containers, _traversals, _loadEstimators, _dataLayouts, _newton3Options, _cellSizeFactors);
        },
        value);
  }

  [[nodiscard]] bool matches(const Configuration &configuration) const {
    return (_containers.empty() or contains(_containers, configuration.container)) and
           (_traversals.empty() or contains(_traversals, configuration.traversal)) and
           (_loadEstimators.empty() or contains(_loadEstimators, configuration.loadEstimator)) and
           (_dataLayouts.empty() or contains(_dataLayouts, configuration.dataLayout)) and
           (_newton3Options.empty() or contains(_newton3Options, configuration.newton3)) and
           (_cellSizeFactors.empty() or contains(_cellSizeFactors, configuration.cellSizeFactor));
  }

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

namespace rule_syntax {

class CodeGenerationContext;

struct Define;

struct Statement {
  virtual void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const = 0;
};

inline Type typeOf(const RuleVM::MemoryCell &memoryCell) { return static_cast<Type>(memoryCell.index()); }

struct Expression {
  [[nodiscard]] virtual Type getType() const = 0;
  virtual void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const = 0;
};

struct DefineList;

struct ConfigurationOrder : public Statement {
  ConfigurationPattern greater;
  ConfigurationPattern smaller;

  enum class SameProperty { container, traversal, dataLayout, newton3, loadEstimator, cellSizeFactor };

  std::vector<SameProperty> sameProperties;

  ConfigurationOrder() = default;

  ConfigurationOrder(ConfigurationPattern greater, ConfigurationPattern smaller, std::vector<SameProperty> same = {})
      : greater(std::move(greater)), smaller(std::move(smaller)), sameProperties(std::move(same)) {}

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override;

  [[nodiscard]] std::string toString() const {
    std::string sameString = sameProperties.empty() ? "" : "with same ";
    for (auto same : sameProperties) {
      sameString += samePropertyToString(same) + ", ";
    }
    return greater.toString() + " >= " + smaller.toString() + sameString;
  }

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

class CodeGenerationContext {
 public:
  explicit CodeGenerationContext(std::map<std::string, std::pair<const Define *, size_t>> initialAddressEnvironment)
      : currentNeededStack(0),
        maxNeededStack(0),
        addressEnvironment(std::move(initialAddressEnvironment)),
        initialNumVariables(addressEnvironment.size()) {}

  void allocateStack(size_t num) {
    currentNeededStack += num;
    maxNeededStack = std::max(maxNeededStack, currentNeededStack);
  }

  void freeStack(size_t num) { currentNeededStack -= num; }

  void addLocalVariable(const Define &definition);

  void addGlobalVariable(const Define &definition);

  [[nodiscard]] size_t addressOf(const std::string &name) const;

  [[nodiscard]] const Define *definitionOf(const std::string &name) const;

  [[nodiscard]] auto getMaxStackSize() const { return maxNeededStack; }

  void addList(const std::string &name, const DefineList *list) { lists.insert({name, list}); }

  auto getList(const std::string &name) { return lists.at(name); }

  auto addConfigurationOrder(const ConfigurationOrder &configurationOrder) {
    configurationOrders.emplace_back(configurationOrder);
    return configurationOrders.size() - 1;
  }

  [[nodiscard]] const auto &getConfigurationOrders() const { return configurationOrders; }

  [[nodiscard]] auto smallerConfigurationPatternByIndex(size_t idx) const {
    return configurationOrders.at(idx).smaller;
  }

  [[nodiscard]] auto getNumLocalVariables() const { return addressEnvironment.size() - initialNumVariables; }

 private:
  size_t currentNeededStack;
  size_t maxNeededStack;

  std::map<std::string, std::pair<const Define *, size_t>> addressEnvironment;
  size_t initialNumVariables;
  std::map<std::string, const DefineList *> lists;
  std::vector<ConfigurationOrder> configurationOrders;
};

struct Literal : public Expression {
  RuleVM::MemoryCell value;
  Type type;

  explicit Literal(RuleVM::MemoryCell value) : value(value), type(typeOf(value)) {}

  [[nodiscard]] Type getType() const override { return type; }

  void generateCode(CodeGenerationContext &context, RuleVM::Program &program) const override {
    program.instructions.push_back({RuleVM::LOADC, value});
    context.allocateStack(1);
  }
};

struct DefineList : public Statement {
  std::string listName;
  std::vector<Literal> values;

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

struct BinaryOperator : public Expression {
  enum Operator { LESS, GREATER, AND };

  std::shared_ptr<Expression> left;
  Operator op;
  std::shared_ptr<Expression> right;

  BinaryOperator(Operator op, std::shared_ptr<Expression> left, std::shared_ptr<Expression> right)
      : left(std::move(left)), op(op), right(std::move(right)) {}

  [[nodiscard]] Type getType() const override { return Type::BOOL; }

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
    program.instructions.push_back({RuleVM::STOREA, context.addressOf(variable)});
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
    program.instructions.push_back({RuleVM::JUMPZERO});
    context.freeStack(1);
    auto jumpInstrIdx = program.instructions.size() - 1;

    for (const auto &statement : consequences) {
      statement->generateCode(context, program);
    }

    program.instructions[jumpInstrIdx].payload = program.instructions.size();
  }
};

struct RuleBasedProgramTree {
  std::vector<std::shared_ptr<Statement>> statements;

  [[nodiscard]] RuleVM::Program generateCode(CodeGenerationContext &context) const {
    RuleVM::Program program;
    program.instructions.push_back({RuleVM::RESERVE, 0ul});
    auto reserveInstructionIdx = program.instructions.size() - 1;
    for (const auto &statement : statements) {
      statement->generateCode(context, program);
    }
    program.instructions[reserveInstructionIdx].payload = context.getNumLocalVariables();
    program.neededStackSize = context.getNumLocalVariables() + context.getMaxStackSize();
    program.instructions.push_back({RuleVM::HALT});
    return program;
  }
};

}  // namespace rule_syntax
}  // namespace autopas

#include "RuleBasedProgramTree.h"

namespace autopas {
ConfigurationPattern::ConfigurationPattern() = default;

ConfigurationPattern::ConfigurationPattern(const ConfigurationPattern &configurationPattern) = default;

ConfigurationPattern::~ConfigurationPattern() noexcept = default;

namespace RuleSyntax {
void Variable::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
  program.instructions.emplace_back(RuleVM::LOADA, context.addressOf(definition->variable));
  context.allocateStack(1);
}

Type Variable::getType() const { return definition->value->getType(); }

UnaryOperator::UnaryOperator(UnaryOperator::Operator op, std::shared_ptr<Expression> child)
    : child(std::move(child)), op(op) {}

void UnaryOperator::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
  child->generateCode(context, program);

  program.instructions.emplace_back(RuleVM::NOT);
}

Type UnaryOperator::getType() const { return Type::BOOL; }

void BinaryOperator::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
  left->generateCode(context, program);
  right->generateCode(context, program);

  static const std::map<Operator, RuleVM::CMD> opMap{
      {LESS, RuleVM::LESS}, {GREATER, RuleVM::GREATER}, {AND, RuleVM::AND}, {OR, RuleVM::OR},
      {ADD, RuleVM::ADD},   {SUB, RuleVM::SUB},         {MUL, RuleVM::MUL}, {DIV, RuleVM::DIV}};
  auto opCode = opMap.at(op);
  program.instructions.emplace_back(opCode);

  context.freeStack(1);
}

BinaryOperator::BinaryOperator(BinaryOperator::Operator op, std::shared_ptr<Expression> left,
                               std::shared_ptr<Expression> right)
    : left(std::move(left)), op(op), right(std::move(right)) {}

Type BinaryOperator::getType() const {
  return op == LESS or op == GREATER or op == AND or op == OR
             ? Type::BOOL
             : (left->getType() == Type::SIZE_T and right->getType() == Type::SIZE_T ? Type::SIZE_T : Type::DOUBLE);
}

ConfigurationOrder::ConfigurationOrder(ConfigurationPattern greater, ConfigurationPattern smaller,
                                       std::vector<SameProperty> same)
    : greater(std::move(greater)), smaller(std::move(smaller)), sameProperties(std::move(same)) {}

void ConfigurationOrder::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
  auto idx = context.addConfigurationOrder(*this);
  program.instructions.emplace_back(RuleVM::OUTPUTC, idx);
}

std::string ConfigurationOrder::toString() const {
  std::string sameString = sameProperties.empty() ? "" : "with same ";
  for (auto same : sameProperties) {
    sameString += samePropertyToString(same) + ", ";
  }
  return greater.toString() + " >= " + smaller.toString() + sameString;
}

bool ConfigurationOrder::haveEqualSameProperties(const Configuration &conf1, const Configuration &conf2) const {
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

std::string ConfigurationOrder::samePropertyToString(ConfigurationOrder::SameProperty same) {
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

CodeGenerationContext::CodeGenerationContext(
    std::map<std::string, std::pair<const Define *, size_t>> initialAddressEnvironment)
    : currentNeededStack(0),
      maxNeededStack(0),
      addressEnvironment(std::move(initialAddressEnvironment)),
      initialNumVariables(addressEnvironment.size()) {}

CodeGenerationContext::~CodeGenerationContext() = default;

void CodeGenerationContext::addLocalVariable(const Define &definition) {
  addressEnvironment.insert({definition.variable, {&definition, addressEnvironment.size()}});
}

void CodeGenerationContext::addGlobalVariable(const Define &definition) {
  initialNumVariables++;
  addLocalVariable(definition);
}

[[nodiscard]] size_t CodeGenerationContext::addressOf(const std::string &name) const {
  return addressEnvironment.at(name).second;
}

const Define *CodeGenerationContext::definitionOf(const std::string &name) const {
  return addressEnvironment.at(name).first;
}

void CodeGenerationContext::allocateStack(size_t num) {
  currentNeededStack += num;
  maxNeededStack = std::max(maxNeededStack, currentNeededStack);
}

void CodeGenerationContext::freeStack(size_t num) { currentNeededStack -= num; }

void CodeGenerationContext::addList(const std::string &name, const DefineList *list) { lists.insert({name, list}); }

const DefineList *CodeGenerationContext::getList(const std::string &name) { return lists.at(name); }

size_t CodeGenerationContext::addConfigurationOrder(const ConfigurationOrder &configurationOrder) {
  configurationOrders.emplace_back(configurationOrder);
  return configurationOrders.size() - 1;
}

const std::vector<ConfigurationOrder> &CodeGenerationContext::getConfigurationOrders() const {
  return configurationOrders;
}

ConfigurationPattern CodeGenerationContext::smallerConfigurationPatternByIndex(size_t idx) const {
  return configurationOrders.at(idx).smaller;
}

size_t CodeGenerationContext::getNumLocalVariables() const { return addressEnvironment.size() - initialNumVariables; }

Literal::Literal(RuleVM::MemoryCell value) : value(value), type(typeOf(value)) {}

Type Literal::getType() const { return type; }

void Literal::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
  program.instructions.emplace_back(RuleVM::LOADC, value);
  context.allocateStack(1);
}

DefineList::DefineList(std::string listName, std::vector<Literal> values)
    : listName(std::move(listName)), values(std::move(values)) {}

void DefineList::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
  context.addList(listName, this);
}

Define::Define(std::string variable, std::shared_ptr<Expression> value)
    : variable(std::move(variable)), value(std::move(value)) {}

Define::~Define() {}

void Define::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
  value->generateCode(context, program);

  context.addLocalVariable(*this);
  program.instructions.emplace_back(RuleVM::STOREA, context.addressOf(variable));
  context.freeStack(1);
}

If::If(std::shared_ptr<Expression> condition, std::vector<std::shared_ptr<Statement>> consequences)
    : condition(std::move(condition)), consequences(std::move(consequences)) {}

void If::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
  condition->generateCode(context, program);
  program.instructions.emplace_back(RuleVM::JUMPZERO);
  context.freeStack(1);
  auto jumpInstrIdx = program.instructions.size() - 1;

  for (const auto &statement : consequences) {
    statement->generateCode(context, program);
  }

  program.instructions[jumpInstrIdx].payload = program.instructions.size();
}

RuleVM::Program RuleBasedProgramTree::generateCode(CodeGenerationContext &context) const {
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
}  // namespace RuleSyntax
}  // namespace autopas
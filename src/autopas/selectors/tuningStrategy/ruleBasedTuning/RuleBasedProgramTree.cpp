#include "RuleBasedProgramTree.h"

namespace autopas::rule_syntax {
  void Variable::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
    program.instructions.push_back({RuleVM::LOADA, context.addressOf(definition->variable)});
    context.allocateStack(1);
  }

  Type Variable::getType() const {
    return std::visit([this](const auto& expr) {return expr.getType(); }, definition->value);
  }

  void BinaryOperator::generateCode(CodeGenerationContext &context, RuleVM::Program &program) const {
    std::visit([&](const auto& expr) {expr.generateCode(context, program); }, *left);
    std::visit([&](const auto& expr) {expr.generateCode(context, program); }, *right);

    auto opCode = op == LESS ? RuleVM::LESS : (op == GREATER ? RuleVM::GREATER : RuleVM::AND);
    program.instructions.push_back({opCode});

    context.freeStack(1);
  }
} // namespace autopas::rule_syntax
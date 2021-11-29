#pragma once

#include <cstddef>
#include <variant>

#include "autopas/selectors/Configuration.h"

namespace autopas {
class RuleVM {
 public:
  using MemoryCell = std::variant<bool, double, size_t, ContainerOption, TraversalOption, LoadEstimatorOption,
                                  DataLayoutOption, Newton3Option>;

  enum CMD {
    LOADC,
    LOADA,
    STOREA,
    RESERVE,
    LESS,
    GREATER,
    EQUAL,
    JUMPZERO,
    OUTPUTC,
    CONDOUTPUTC,
    HALT,
    AND,
    OR,
    POP,
    MUL,
    DIV,
    ADD,
    SUB,
    NOT
  };

  struct Instruction {
    CMD cmd;
    MemoryCell payload;

    explicit Instruction(CMD cmd, MemoryCell payload = MemoryCell{0ul}) : cmd(cmd), payload(payload) {}
  };

  struct Program {
    std::vector<Instruction> instructions;
    size_t neededStackSize;
  };

  std::vector<size_t> execute(const Program &program, const std::vector<MemoryCell> &initialStack) {
    _programCounter = 0;
    _removedPatterns.clear();
    _stack = initialStack;
    _stack.resize(program.neededStackSize + initialStack.size());
    _stackPointer = initialStack.size() - 1;
    _halt = false;

    while (not _halt) {
      executeInstruction(program.instructions.at(_programCounter++));
    }

    return _removedPatterns;
  }

 private:
  void executeInstruction(Instruction instruction) {
    switch (instruction.cmd) {
      case LOADC:
        _stack.at(++_stackPointer) = instruction.payload;
        break;
      case LOADA:
        _stack.at(++_stackPointer) = _stack.at(std::get<size_t>(instruction.payload));
        break;
      case STOREA:
        _stack.at(std::get<size_t>(instruction.payload)) = _stack.at(_stackPointer--);
        break;
      case RESERVE:
        _stackPointer += std::get<size_t>(instruction.payload);
        break;
      case LESS: {
        bool res = compare<std::less>();
        _stackPointer--;
        _stack.at(_stackPointer) = res;
        break;
      }
      case GREATER: {
        bool res = compare<std::greater>();
        _stackPointer--;
        _stack.at(_stackPointer) = res;
        break;
      }
      case EQUAL: {
        bool res = compare<std::equal_to>();
        _stackPointer--;
        _stack.at(_stackPointer) = res;
        break;
      }
      case JUMPZERO: {
        bool shouldJump = not std::get<bool>(_stack.at(_stackPointer--));
        if (shouldJump) {
          _programCounter = std::get<size_t>(instruction.payload);
        }
        break;
      }
      case OUTPUTC:
        _removedPatterns.push_back(std::get<size_t>(instruction.payload));
        break;
      case CONDOUTPUTC:
        if (std::get<bool>(_stack.at(_stackPointer))) {
          _removedPatterns.push_back(std::get<size_t>(instruction.payload));
        }
        break;
      case HALT:
        _halt = true;
        break;
      case AND: {
        bool res = std::get<bool>(_stack.at(_stackPointer)) and std::get<bool>(_stack.at(_stackPointer - 1));
        _stack.at(--_stackPointer) = res;
        break;
      } break;
      case OR: {
        bool res = std::get<bool>(_stack.at(_stackPointer)) or std::get<bool>(_stack.at(_stackPointer - 1));
        _stack.at(--_stackPointer) = res;
        break;
      }
      case POP:
        _stackPointer--;
        break;
      case MUL: {
        auto res =
            computeBinary(_stack.at(_stackPointer - 1), _stack.at(_stackPointer), [](auto l, auto r) { return l * r; });
        _stack.at(--_stackPointer) = res;
        break;
      }
      case DIV: {
        auto res =
            computeBinary(_stack.at(_stackPointer - 1), _stack.at(_stackPointer), [](auto l, auto r) { return l / r; });
        _stack.at(--_stackPointer) = res;
        break;
      }
      case ADD: {
        auto res =
            computeBinary(_stack.at(_stackPointer - 1), _stack.at(_stackPointer), [](auto l, auto r) { return l + r; });
        _stack.at(--_stackPointer) = res;
        break;
      }
      case SUB: {
        auto res =
            computeBinary(_stack.at(_stackPointer - 1), _stack.at(_stackPointer), [](auto l, auto r) { return l - r; });
        _stack.at(--_stackPointer) = res;
        break;
      }
      case NOT:
        _stack.at(_stackPointer) = not std::get<bool>(_stack.at(_stackPointer));
        break;
    }
  }

  template <class T>
  static constexpr auto isNumericVal() {
    using type = std::remove_cv_t<std::remove_reference_t<T>>;
    return std::is_same_v<type, double> or std::is_same_v<type, size_t>;
  }

  template <template <class> typename Compare>
  [[nodiscard]] bool compare() {
    bool res = std::visit(
        [](auto &&left, auto &&right) {
          if constexpr (std::is_same_v<decltype(left), decltype(right)>) {
            return Compare{}(left, right);
          } else if constexpr (isNumericVal<decltype(left)>() and isNumericVal<decltype(right)>()) {
            return Compare{}(left, right);
          } else {
            throw std::runtime_error("RuleVM: cannot compare");
            return false;
          }
        },
        _stack.at(_stackPointer - 1), _stack.at(_stackPointer));
    return res;
  }

  template <class Functor>
  [[nodiscard]] MemoryCell computeBinary(const MemoryCell &left, const MemoryCell &right, Functor op) {
    return std::visit(
        [&op](auto &&l, auto &&r) {
          if constexpr (isNumericVal<decltype(l)>() and isNumericVal<decltype(r)>()) {
            return MemoryCell{op(l, r)};
          } else {
            throw std::runtime_error("RuleVM: cannot compute binary operator with these operators");
            return MemoryCell{};
          }
        },
        left, right);
  }

 private:
  size_t _programCounter;
  size_t _stackPointer;
  bool _halt;

  std::vector<MemoryCell> _stack;
  std::vector<size_t> _removedPatterns;
};
}  // namespace autopas

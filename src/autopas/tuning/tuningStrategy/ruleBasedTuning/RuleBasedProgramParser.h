#pragma once

#include <any>

#include "RuleBasedProgramTree.h"

namespace autopas::rule_syntax {

/**
 * Parses a rule program and produces an AST (RuleBasedProgramTree) and a corresponding CodeGenerationContext.
 */
class RuleBasedProgramParser {
 public:
  /**
   * Creates a RuleBasedProgramParser with some predefined variables that can be used in the program.
   * @param initialDefinitions The predefined variables that can be used in the program. Typically, this is live info.
   */
  explicit RuleBasedProgramParser(std::vector<std::pair<std::string, Define>> &initialDefinitions)
      : _initialDefinitions(initialDefinitions) {}

  /**
   * Parses a rule based program given as a string.
   * @param programCode The rule based program to parse.
   * @return The AST and a corresponding CodeGenerationContext.
   */
  std::pair<RuleBasedProgramTree, CodeGenerationContext> parse(const std::string &programCode);

 private:
  /**
   * The predefined variables that can be used in a parsed rule program.
   */
  std::vector<std::pair<std::string, Define>> &_initialDefinitions;
  /**
   * The defined lists with their names in the parsed program.
   */
  std::map<std::string, DefineList> _lists;
};

}  // namespace autopas::rule_syntax
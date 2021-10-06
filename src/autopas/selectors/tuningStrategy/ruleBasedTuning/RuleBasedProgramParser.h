#pragma once

#include <any>

#include "RuleBasedProgramTree.h"

namespace autopas::rule_syntax {

class RuleBasedProgramParser {
 public:
  explicit RuleBasedProgramParser(std::vector<std::pair<std::string, Define>> &initialDefinitions)
      : _initialDefinitions(initialDefinitions) {}

  std::pair<RuleBasedProgramTree, CodeGenerationContext> parse(const std::string &programCode);

 private:
  std::vector<std::pair<std::string, Define>> &_initialDefinitions;
  std::map<std::string, DefineList> _lists;
};

}  // namespace autopas::rule_syntax
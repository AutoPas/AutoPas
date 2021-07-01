/**
 * @file RuleBasedTuningTest.cpp
 * @author humig
 * @date 30.06.2021
 */

#include "RuleBasedTuningTest.h"

TEST_F(RuleBasedTuningTest, testParser) {
  autopas::rule_syntax::RuleBasedProgramParser::test();
}

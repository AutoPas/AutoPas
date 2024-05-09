/**
 * @file FuzzyRuleErrorListener.h
 * @author Manuel Lerchner
 * @date 09.05.24
 */

#pragma once

#include "antlr4-runtime.h"

namespace autopas_generated_fuzzy_rule_syntax {

using namespace antlr4;

/**
 * This class implements an error listener for the fuzzy rule parser.
 * It catches all errors and throws a runtime error.
 */
class FuzzyRuleErrorListener : public antlr4::BaseErrorListener {
 public:
  void syntaxError(antlr4::Recognizer *recognizer, antlr4::Token *offendingSymbol, size_t line,
                   size_t charPositionInLine, const std::string &msg, std::exception_ptr e) override;

  void reportAmbiguity(antlr4::Parser *recognizer, const antlr4::dfa::DFA &dfa, size_t startIndex, size_t stopIndex,
                       bool exact, const antlrcpp::BitSet &ambigAlts, antlr4::atn::ATNConfigSet *configs) override;

  void reportAttemptingFullContext(antlr4::Parser *recognizer, const antlr4::dfa::DFA &dfa, size_t startIndex,
                                   size_t stopIndex, const antlrcpp::BitSet &conflictingAlts,
                                   antlr4::atn::ATNConfigSet *configs) override;

  void reportContextSensitivity(antlr4::Parser *recognizer, const antlr4::dfa::DFA &dfa, size_t startIndex,
                                size_t stopIndex, size_t prediction, antlr4::atn::ATNConfigSet *configs) override;
};

}  // namespace autopas_generated_fuzzy_rule_syntax
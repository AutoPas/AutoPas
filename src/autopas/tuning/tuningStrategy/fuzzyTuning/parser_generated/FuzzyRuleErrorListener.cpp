/**
 * @file FuzzyRuleErrorListener.cpp
 * @author Manuel Lerchner
 * @date 09.05.24
 */

#include "FuzzyRuleErrorListener.h"

namespace AutopasGeneratedFuzzyRuleSyntax {

void FuzzyRuleErrorListener::syntaxError(antlr4::Recognizer *recognizer, antlr4::Token *offendingSymbol, size_t line,
                                         size_t charPositionInLine, const std::string &msg, std::exception_ptr e) {
  std::string error = "line " + std::to_string(line) + ":" + std::to_string(charPositionInLine) + " " + msg;

  throw std::runtime_error(error);
}

void FuzzyRuleErrorListener::reportAmbiguity(antlr4::Parser *recognizer, const antlr4::dfa::DFA &dfa, size_t startIndex,
                                             size_t stopIndex, bool exact, const antlrcpp::BitSet &ambigAlts,
                                             antlr4::atn::ATNConfigSet *configs) {
  throw std::runtime_error("Ambiguity in fuzzy rule");
}

void FuzzyRuleErrorListener::reportAttemptingFullContext(antlr4::Parser *recognizer, const antlr4::dfa::DFA &dfa,
                                                         size_t startIndex, size_t stopIndex,
                                                         const antlrcpp::BitSet &conflictingAlts,
                                                         antlr4::atn::ATNConfigSet *configs) {
  throw std::runtime_error("Attempting full context in fuzzy rule");
}

void FuzzyRuleErrorListener::reportContextSensitivity(antlr4::Parser *recognizer, const antlr4::dfa::DFA &dfa,
                                                      size_t startIndex, size_t stopIndex, size_t prediction,
                                                      antlr4::atn::ATNConfigSet *configs) {
  throw std::runtime_error("Context sensitivity in fuzzy rule");
}

}  // namespace AutopasGeneratedFuzzyRuleSyntax

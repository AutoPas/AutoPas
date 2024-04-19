/**
 * @file LinguisticVariable.h
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#pragma once

#include <memory>
#include <vector>

#include "CrispSet.h"
#include "FuzzySet.h"

namespace autopas::fuzzy_logic {

class LinguisticVariable {
 public:
  /**
   * Constructs a LinguisticVariable with the given name.
   * @param name
   */
  explicit LinguisticVariable(const std::string &name, const std::pair<double, double> &range);

  /**
   * Adds a FuzzySet to the LinguisticVariable.
   * @param fuzzySet
   */
  void addLinguisticTerm(const std::shared_ptr<FuzzySet> &fuzzySet);

  /**
   * Overload of the operator== to get a FuzzySet by its name. This allows a very concise syntax to create fuzzy rules.
   * @param linguisticTerm
   * @return The FuzzySet with the given name.
   */
  std::shared_ptr<FuzzySet> operator==(const std::string &linguisticTerm) const;

 private:
  /**
   * The CrispSet on which the LinguisticVariable is defined.
   */
  std::shared_ptr<CrispSet> _crispSet;

  /**
   * All linguistic terms of the LinguisticVariable organized by their name.
   */
  std::map<std::string, std::shared_ptr<FuzzySet>> _linguisticTerms;
};
}  // namespace autopas::fuzzy_logic
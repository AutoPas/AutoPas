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

namespace autopas::FuzzyLogic {

/**
 * A class representing a LinguisticVariable. A LinguisticVariable is defined on a CrispSet and consists of several
 * FuzzySets, which are the linguistic terms of the LinguisticVariable.
 */
class LinguisticVariable {
 public:
  /**
   * Constructs a LinguisticVariable with the given name and range.
   * @param name The name of the LinguisticVariable.
   * @param range The range of the LinguisticVariable in the form [min, max].
   */
  explicit LinguisticVariable(const std::string &name, std::pair<double, double> range);

  /**
   * Adds a new linguistic term to the LinguisticVariable.
   * Additionally updates the CrispSet of the linguistic term to the current CrispSet.
   * @param linguisticTerm The linguistic term to add.
   */
  void addLinguisticTerm(const std::shared_ptr<FuzzySet> &linguisticTerm);

  /**
   * Overload of the operator== where the left-hand side is a linguistic variable, and the right-hand side is a
   * linguistic term. Returns the fuzzy set corresponding to a linguistic term. This allows a very concise syntax to
   * create fuzzy rules.
   * @param linguisticTerm The name of the linguistic term to extract.
   * @return The FuzzySet with the given name.
   */
  std::shared_ptr<FuzzySet> operator==(const std::string &linguisticTerm) const;

  /**
   * Getter for the name of the LinguisticVariable.
   * @return The name of the LinguisticVariable.
   */
  [[nodiscard]] const std::string &getName() const;

  /**
   *  Returns a string representation of the LinguisticVariable.
   */
  explicit operator std::string() const;

 private:
  /**
   * The name of the LinguisticVariable.
   */
  std::string _name;

  /**
   * The CrispSet on which the LinguisticVariable is defined.
   */
  std::shared_ptr<CrispSet> _crispSet;

  /**
   * All linguistic terms of the LinguisticVariable organized by their name.
   */
  std::map<std::string, std::shared_ptr<FuzzySet>> _linguisticTerms;
};
}  // namespace autopas::FuzzyLogic
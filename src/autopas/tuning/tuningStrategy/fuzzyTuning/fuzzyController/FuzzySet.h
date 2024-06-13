/**
 * @file FuzzySet.h
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <variant>

#include "CrispSet.h"

namespace autopas::fuzzy_logic {

/**
 * Used to represent a mathematical Fuzzy-Set.
 */
class FuzzySet {
 public:
  /**
   * A map of the form {dimension_name: value}. Used to represent higher dimensional function inputs.
   */
  using Data = std::map<std::string, double>;

  /**
   * A composed membership function is a membership function that relies on different FuzzySets to calculate the
   * membership value. This is important as it allows for straightforward implementation of fuzzy logic operators.
   * Composed membership represent higher dimensional membership functions (m=f(x1,x2,...,xn)). The necessary data is
   * passed as a map of the form {dimension_name: value}.
   */
  using ComposedMembershipFunction = std::function<double(const Data &)>;

  /**
   * A base membership function calculates the membership value based on a direct function evaluation.
   * Base membership function take a single value and returns the corresponding membership value. (m=f(x))
   * Additionally, the base membership function stores information about the membership function and its parameters.
   */
  using BaseMembershipFunction = std::tuple<std::string, std::vector<double>, std::function<double(double)>>;

  /**
   * Constructs a FuzzySet with a given linguistic term and base membership function.
   * @param linguisticTerm
   * @param baseMembershipFunction
   */
  FuzzySet(std::string linguisticTerm, BaseMembershipFunction &&baseMembershipFunction);

  /**
   * Constructs a FuzzySet with a given linguistic term, composed membership function and crisp set.
   * @param linguisticTerm
   * @param membershipFunction
   * @param crispSet
   */
  FuzzySet(std::string linguisticTerm, ComposedMembershipFunction &&membershipFunction,
           const std::shared_ptr<CrispSet> &crispSet);

  /**
   * Evaluates the membership function of this FuzzySet at the given value.
   * @param data A map of the form {dimension_name: value}.
   * @return The membership value of the given value in this FuzzySet.
   */
  [[nodiscard]] double evaluate_membership(const Data &data) const;

  /**
   * Calculates the x-coordinate of the center of gravity of this FuzzySet.
   * @param numSamples The number of samples to use for the numerical CoG calculation.
   * @return The x-coordinate of the CoG of this FuzzySet.
   */
  [[nodiscard]] double CoG(size_t numSamples) const;

  /**
   * Calculates the mean of the meanOfMaximum of this FuzzySet.
   * @param numSamples The number of samples to use for the numerical mean of meanOfMaximum calculation.
   * @return The mean of the meanOfMaximum of this FuzzySet.
   */
  [[nodiscard]] double meanOfMaximum(size_t numSamples) const;

  /**
   * Returns a string representation of the BaseMembershipFunction.
   * @return A string representation of the BaseMembershipFunction.
   */
  [[nodiscard]] std::string printBaseMembershipFunction() const;

  /**
   * Returns a string representation of the FuzzySet.
   */
  explicit operator std::string() const;

  /**
   * Returns the linguistic term of the FuzzySet.
   * @return The linguistic term of the FuzzySet.
   */
  [[nodiscard]] const std::string &getLinguisticTerm() const;

  /**
   * Returns the crisp set of the FuzzySet.
   * @return The crisp set of the FuzzySet.
   */
  [[nodiscard]] const std::shared_ptr<CrispSet> &getCrispSet() const;

  /**
   * Sets the crisp set of the FuzzySet.
   * @param crispSet
   */
  void setCrispSet(const std::shared_ptr<CrispSet> &crispSet);

  /**
   * Make the Union Operator a friend of the class.
   */
  friend std::shared_ptr<FuzzySet> operator||(const std::shared_ptr<FuzzySet> &lhs,
                                              const std::shared_ptr<FuzzySet> &rhs);

  /**
   * Make the Intersection Operator a friend of the class.
   */
  friend std::shared_ptr<FuzzySet> operator&&(const std::shared_ptr<FuzzySet> &lhs,
                                              const std::shared_ptr<FuzzySet> &rhs);

  /**
   * Make the Complement Operator a friend of the class.
   */
  friend std::shared_ptr<FuzzySet> operator!(const std::shared_ptr<FuzzySet> &fuzzySet);

 private:
  /**
   * The linguistic term of the FuzzySet.
   */
  std::string _linguisticTerm;

  /**
   * The membership function of the FuzzySet. Can be a base or composed membership function.
   */
  std::variant<BaseMembershipFunction, ComposedMembershipFunction> _membershipFunction;

  /**
   * The crisp set on which the FuzzySet is defined.
   */
  std::shared_ptr<CrispSet> _crispSet;
};

/**
 * Returns the union of two FuzzySets.
 * @param lhs
 * @param rhs
 * @return The union of lhs and rhs.
 */
std::shared_ptr<FuzzySet> operator||(const std::shared_ptr<FuzzySet> &lhs, const std::shared_ptr<FuzzySet> &rhs);

/**
 * Returns the intersection of two FuzzySets.
 * @param lhs
 * @param rhs
 * @return The intersection of lhs and rhs.
 */
std::shared_ptr<FuzzySet> operator&&(const std::shared_ptr<FuzzySet> &lhs, const std::shared_ptr<FuzzySet> &rhs);

/**
 * Returns the complement of a FuzzySet.
 * @param fuzzySet
 * @return The complement of fuzzySet.
 */
std::shared_ptr<FuzzySet> operator!(const std::shared_ptr<FuzzySet> &fuzzySet);

}  // namespace autopas::fuzzy_logic

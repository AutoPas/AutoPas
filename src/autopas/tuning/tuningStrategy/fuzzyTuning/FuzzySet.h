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
  using ComposedMembershipFunction = std::function<double(const std::map<std::string, double> &)>;
  using BaseMembershipFunction = std::function<double(double)>;

  /**
   * Constructs a FuzzySet with the given linguistic term and membership function.
   * @param linguisticTerm
   * @param membershipFunction
   */
  FuzzySet(std::string linguisticTerm, const std::shared_ptr<ComposedMembershipFunction> &membershipFunction);

  /**
   * Constructs a FuzzySet with the given linguistic term and membership function.
   * @param linguisticTerm
   * @param membershipFunction
   */
  FuzzySet(std::string linguisticTerm, const std::shared_ptr<BaseMembershipFunction> &baseMembershipFunction);

  /**
   * Constructs a FuzzySet with the given linguistic term and crisp set.
   * @param linguisticTerm
   * @param membershipFunction
   * @param crispSet
   */
  FuzzySet(std::string linguisticTerm, const std::shared_ptr<ComposedMembershipFunction> &membershipFunction,
           const std::shared_ptr<CrispSet> &crispSet);

  /**
   * Evaluates the membership function of this FuzzySet at the given value.
   * @param data A map of the form {dimension_name: value}.
   * @return The membership value of the given value in this FuzzySet.
   */
  [[nodiscard]] double evaluate_membership(const std::map<std::string, double> &data) const;

  /**
   * Calculates the x-coordinate of the centroid of this FuzzySet.
   * @return The x-coordinate of the centroid of this FuzzySet.
   */
  [[nodiscard]] double centroid(size_t numSamples = 100) const;

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
  friend std::shared_ptr<FuzzySet> operator!(const std::shared_ptr<FuzzySet> &lhs);

 private:
  std::string _linguisticTerm;
  std::shared_ptr<ComposedMembershipFunction> _membershipFunction;
  std::optional<std::shared_ptr<BaseMembershipFunction>> _baseMembershipFunction;
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
 * @param set
 * @return The complement of lhs.
 */
std::shared_ptr<FuzzySet> operator!(const std::shared_ptr<FuzzySet> &set);

}  // namespace autopas::fuzzy_logic

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
  using MembershipFunction = std::function<double(std::map<std::string, double>)>;
  using BaseMembershipFunction = std::function<double(double)>;

  /**
   * Constructs a FuzzySet with the given linguistic term and membership function.
   * @param linguisticTerm
   * @param membershipFunction
   */
  FuzzySet(const std::string &linguisticTerm, const MembershipFunction &membershipFunction);

  /**
   * Constructs a FuzzySet with the given linguistic term and membership function.
   * @param linguisticTerm
   * @param membershipFunction
   */
  FuzzySet(const std::string &linguisticTerm, const BaseMembershipFunction &membershipFunction);

  /**
   * Constructs a FuzzySet with the given linguistic term and crisp set.
   * @param linguisticTerm
   * @param membershipFunction
   * @param crispSet
   */
  FuzzySet(const std::string &linguisticTerm, const MembershipFunction &membershipFunction,
           const std::shared_ptr<CrispSet> &crispSet);

  /**
   * Evaluates the membership function of this FuzzySet at the given value.
   * @param data A map of the form {dimension_name: value}.
   * @return The membership value of the given value in this FuzzySet.
   */
  double evaluate_membership(const std::map<std::string, double> &data) const;

  /**
   * Calculates the x-coordinate of the centroid of this FuzzySet.
   * @return The x-coordinate of the centroid of this FuzzySet.
   */
  double centroid(size_t numSamples = 100) const;

  /**
   * Calculates the intersection of two FuzzySets.
   * @param rhs
   * @return A new FuzzySet, which is the intersection of this and rhs.
   */
  FuzzySet operator&&(const FuzzySet &rhs) const;

  /**
   * Calculates the union of two FuzzySets.
   * @param rhs
   * @return A new FuzzySet, which is the union of this and rhs.
   */
  FuzzySet operator||(const FuzzySet &rhs) const;

  /**
   * Calculates the complement of this FuzzySet.
   * @return A new FuzzySet, which is the complement of this.
   */
  FuzzySet operator!() const;

  /**
   * Returns the linguistic term of the FuzzySet.
   * @return The linguistic term of the FuzzySet.
   */
  [[nodiscard]] const std::string &getLinguisticTerm() const;

  /**
   * Returns the membership function of the FuzzySet.
   * @return The membership function of the FuzzySet.
   */
  [[nodiscard]] const MembershipFunction &getMembershipFunction() const;

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
   * Sets the base membership function of the FuzzySet.
   * @param baseMembershipFunction
   */
  void setBaseMembershipFunction(const BaseMembershipFunction &baseMembershipFunction);

 private:
  std::string _linguisticTerm;
  MembershipFunction _membershipFunction;
  std::optional<BaseMembershipFunction> _baseMembershipFunction;
  std::shared_ptr<CrispSet> _crispSet;
};

}  // namespace autopas::fuzzy_logic

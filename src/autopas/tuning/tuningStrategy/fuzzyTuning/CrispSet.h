/**
 * @file CrispSet.h
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#pragma once

#include <map>
#include <memory>
#include <string>

namespace autopas::fuzzy_logic {

/**
 * Used to represent arbitrary Crisp-Sets, on which Fuzzy-Sets can be defined.
 * In this implementation, a Crisp-Set is represented by a map of dimensions, where each dimension states a
 * interval of possible values.
 */
class CrispSet {
 public:
  /**
   * Constructs a CrispSet with empty dimensions.
   */
  explicit CrispSet() = default;

  /**
   * Constructs a CrispSet with the given dimensions.
   * @param dimensions
   */
  explicit CrispSet(const std::string &name, const std::pair<double, double> &range);

  /**
   * Calculates the cartesian product of two CrispSets. The resulting CrispSet will contain dimensions from
   * both input CrispSets.
   * @param rhs
   * @return A new CrispSet, which is the cartesian product of this and rhs.
   */
  std::shared_ptr<CrispSet> operator*(CrispSet &rhs) const;

  /**
   * Returns the dimensions of the CrispSet.
   * @return The dimensions of the CrispSet.
   */
  [[nodiscard]] std::map<std::string, std::pair<double, double>> &getDimensions();

 private:
  /**
   * All dimensions of the CrispSet. Each dimension is represented by a tuple of the form of [min, max], which
   * specifies the interval of possible values for this dimension.
   */
  std::map<std::string, std::pair<double, double>> _dimensions;
};

}  // namespace autopas::fuzzy_logic

/**
 * @file VectorizationPatternOption.h
 * @author L. Gall
 * @date 23.7.24
 */

#pragma once

#include <set>
#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {

/**
 * Class representing the choices of possible vectorization patterns for the Pairwise Functors
 */
class VectorizationPatternOption : public Option<VectorizationPatternOption> {
 public:
  /**
   * Possible choices for the vector patterns
   */
  enum Value {
    /**
     * Interact one particle from the first list with the full vector length from the second list
     */
    p1xVec,
    /**
     * Interact two particles from the first list with half the vector length from the second list
     */
    p2xVecDiv2,
    /**
     * Interact half the vector length from the first list with two particles from the second list
     */
    pVecDiv2x2,
    /**
     * Interact the full vector length from the first list with one particle from the second list
     */
    pVecx1
  };

  /**
   * Constructor
   */
  VectorizationPatternOption() = default;

  /**
   * Constructor from value
   * @param option
   */
  constexpr VectorizationPatternOption(Value option) : _value(option) {}

  /**
   * Cast to value
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting
   * @return
   */
  static std::set<VectorizationPatternOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of Vectorization Patterns
   * @return map option -> string representation
   */
  static std::map<VectorizationPatternOption, std::string> getOptionNames() {
    return {
        {VectorizationPatternOption::p1xVec, "1xVectorLength"},
        {VectorizationPatternOption::p2xVecDiv2, "2xVectorLengthDiv2"},
        {VectorizationPatternOption::pVecDiv2x2, "VectorLengthDiv2x2"},
        {VectorizationPatternOption::pVecx1, "VectorLengthx1"},
    };
  }

 private:
  Value _value{Value(-1)};
};

}  // namespace options
}  // namespace autopas
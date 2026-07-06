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
 *
 * Vectorization patterns refer to the specific strategies used to load and
 * arrange particle data into SIMD registers. Different patterns can lead to different
 * SIMD efficiency due to variation in memory access behavior and data locality.
 *
 * For details and benchmarking results, see:
 * https://doi.org/10.48550/arXiv.2512.03565
 */
class VectorizationPatternOption : public Option<VectorizationPatternOption> {
 public:
  /**
   * Possible choices for the vector patterns
   */
  enum Value {
    /**
     * Interact one particle from the first list with the full vector length from the second list
     *
     * ---------------------------------
     * | i | i | i | i | i | i | i | i |
     * ---------------------------------
     * | j |j+1|j+2|j+3|j+4|j+5|j+6|j+7|
     *  --------------------------------
     */
    p1xVec,
    /**
     * Interact two particles from the first list with half the vector length from the second list
     *
     * ---------------------------------
     * | i | i | i | i |i+1|i+1|i+1|i+1|
     * ---------------------------------
     * | j |j+1|j+2|j+3| j |j+1|j+2|j+3|
     *  --------------------------------
     */
    p2xVecDiv2,
    /**
     * Interact half the vector length from the first list with two particles from the second list
     *
     * ---------------------------------
     * | i |i+1|i+2|i+3| i |i+1|i+2|i+3|
     * ---------------------------------
     * | j | j | j | j |j+1|j+1|j+1|j+1|
     *  --------------------------------
     */
    pVecDiv2x2,
    /**
     * Interact the full vector length from the first list with one particle from the second list
     *
     * ---------------------------------
     * | i |i+1|i+2|i+3|i+4|i+5|i+6|i+7|
     * ---------------------------------
     * | j | j | j | j | j | j | j | j |
     *  --------------------------------
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
        {p1xVec, "1xVectorLength"},
        {p2xVecDiv2, "2xVectorLengthDiv2"},
        {pVecDiv2x2, "VectorLengthDiv2x2"},
        {pVecx1, "VectorLengthx1"},
    };
  }

 private:
  Value _value{Value(-1)};
};

}  // namespace options
}  // namespace autopas
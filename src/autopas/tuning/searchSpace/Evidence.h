/**
 * @file Evidence.h
 * @author F. Gratl
 * @date 23.06.23
 */

#pragma once

#include <cstddef>

namespace autopas {
/**
 * Helper class that associates a measurement with the iteration when it was taken.
 */
class Evidence {
 public:
  /**
   * Iteration in which the measurement was taken.
   */
  std::size_t iteration;
  /**
   * Tuning phase in which the measurement was taken.
   */
  std::size_t tuningPhase;
  /**
   * Value of the measurement (time, energy, ...).
   */
  long value;

  /**
   * Equality operator.
   * @param rhs
   * @return
   */
  bool operator==(const Evidence &rhs) const;
  /**
   * Inequality operator.
   * @param rhs
   * @return
   */
  bool operator!=(const Evidence &rhs) const;
};
}  // namespace autopas
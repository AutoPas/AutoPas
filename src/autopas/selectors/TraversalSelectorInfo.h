/**
 * @file TraversalSelectorInfo.h
 * @author humig
 * @date 28.05.19
 */

#pragma once
#include <array>

namespace autopas {
/**
 * Info for traversals of a specific container.
 * @tparam ParticleCell
 */
class TraversalSelectorInfo {
 public:
  /**
   * Dummy constructor such that this class can be used in maps
   */
  TraversalSelectorInfo() : dims({0, 0, 0}), interactionLength(1.0), cellLength({0.0, 0.0, 0.0}) {}

  /**
   * Constructor of the TraversalSelector class.
   * @param dims Array with the dimension lengths of the domain
   * @param interactionLength Interaction length (cutoff radius + skin)
   * @param cellLength cell length.
   */
  explicit TraversalSelectorInfo(const std::array<unsigned long, 3> &dims, const double interactionLength = 1.0,
                                 const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
      : dims(dims), interactionLength(interactionLength), cellLength(cellLength) {}

  /**
   * Array with the dimension lengths of the domain
   * indicating whether or not the optimalTraversalOption is already initialized
   */
  const std::array<unsigned long, 3> dims;

  /**
   * Interaction length
   */
  const double interactionLength;

  /**
   * cell length
   */
  const std::array<double, 3> cellLength;
};
}  // namespace autopas

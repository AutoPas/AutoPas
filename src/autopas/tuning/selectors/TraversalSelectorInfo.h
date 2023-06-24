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
 */
class TraversalSelectorInfo {
 public:
  /**
   * Dummy constructor such that this class can be used in maps
   */
  TraversalSelectorInfo()
      : cellsPerDim({0, 0, 0}), interactionLength(1.0), cellLength({0.0, 0.0, 0.0}), clusterSize{0} {}

  /**
   * Constructor of the TraversalSelector class.
   * @param cellsPerDim Array with the dimension lengths of the domain in number of cells.
   * @param interactionLength Interaction length (cutoff radius + skin)
   * @param cellLength cell length.
   * @param clusterSize The size of a cluster (set this to 0 if not applicable).
   */
  explicit TraversalSelectorInfo(const std::array<unsigned long, 3> &cellsPerDim, const double interactionLength,
                                 const std::array<double, 3> &cellLength, const unsigned int clusterSize)
      : cellsPerDim(cellsPerDim),
        interactionLength(interactionLength),
        cellLength(cellLength),
        clusterSize{clusterSize} {}

  /**
   * Number of cells in the domain per dimension.
   */
  const std::array<unsigned long, 3> cellsPerDim;

  /**
   * Interaction length
   */
  const double interactionLength;

  /**
   * cell length
   */
  const std::array<double, 3> cellLength;

  /**
   * This specifies the size of verlet clusters.
   */
  const unsigned int clusterSize;
};
}  // namespace autopas

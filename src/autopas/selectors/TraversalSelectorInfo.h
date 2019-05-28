/**
 * @file TraversalSelectorInfo.h
 * @author humig
 * @date 28.05.19
 */

#pragma once

namespace autopas {

template <class ParticleCell>
class TraversalSelector;

/**
 * Info for traversals of a specific container.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class TraversalSelectorInfo {
  friend TraversalSelector<ParticleCell>;

 public:
  /**
   * Dummy constructor such that this class can be used in maps
   */
  TraversalSelectorInfo() : _dims({0, 0, 0}), _cutoff(1.0), _cellLength({0.0, 0.0, 0.0}) {}
  /**
   * Constructor of the TraversalSelector class.
   * @param dims Array with the dimension lengths of the domain
   * @param cutoff Cutoff radius
   * @param cellLength cell length.
   */
  explicit TraversalSelectorInfo(const std::array<unsigned long, 3> &dims, const double cutoff = 1.0,
                                 const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
      : _dims(dims), _cutoff(cutoff), _cellLength(cellLength) {}

 private:
  /**
   * indicating whether or not the optimalTraversalOption is already initialized
   */
  const std::array<unsigned long, 3> _dims;

  const double _cutoff;

  const std::array<double, 3> _cellLength;
};
}  // namespace autopas
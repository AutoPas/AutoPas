/**
 * @file TraversalSelectorInfoAdaptive.h
 * @author C.Menges
 * @date 5.07.2019
 */

#pragma once

#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/adaptiveLinkedCells/Octree.h"
#include "autopas/particles/Particle.h"
#include "autopas/selectors/TraversalSelectorInfo.h"

#include <array>

namespace autopas {
/**
 * Info for traversals of the adaptive LC container.
 * @tparam ParticleCell
 */
class TraversalSelectorInfoAdaptive : public TraversalSelectorInfo {
 public:
  /**
   * Constructor of the TraversalSelectorInfoAdaptive class.
   * @param dims Array with the dimension lengths of the domain.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   * @param octree Octree used to organize adaptive cells.
   */
  explicit TraversalSelectorInfoAdaptive(const std::array<unsigned long, 3> &dims, const double cutoff,
                                         const std::array<double, 3> &cellLength/*,
                                         Octree<Particle, ParticleCell<Particle>> &octree*/)
      : TraversalSelectorInfo(dims, cutoff, cellLength)/*, octree(octree) */{}

  /**
   * Octree used to organize adaptive cells.
   */
  // Octree<Particle, ParticleCell<Particle>> &octree;
};
}  // namespace autopas

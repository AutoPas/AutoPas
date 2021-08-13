/**
 * @file LeavingParticleCollector.h
 *
 * @author seckler
 * @date 13.08.2021
 */

#pragma once

#include "autopas/options/IteratorBehavior.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/inBox.h"

namespace autopas {
/**
 * Namespace to collect leaving particles from a container.
 */
namespace LeavingParticleCollector {
/**
 * Collects leaving particles and marks halo particles as dummy.
 * @note This function does not move or actually delete any particles!
 * @return Returns a vector of leaving particles.
 * @todo: openmp parallelization, general optimizations (do not iterate over all particles, but only over boundary,
 * eventually merge both loops?
 */
template <class ContainerType>
std::vector<typename ContainerType::ParticleType> collectParticlesAndMarkNonOwnedAsDummy(ContainerType &container) {
  // First mark halo particles as dummy!
  for (auto iter = container.begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
    iter->setOwnershipState(OwnershipState::dummy);
  }

  std::vector<typename ContainerType::ParticleType> leavingParticles;
  // Collect leaving particles
  for (auto iter = container.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    if (not utils::inBox(iter->getR(), container.getBoxMin(), container.getBoxMax())) {
      leavingParticles.push_back(*iter);
      iter->setOwnershipState(OwnershipState::dummy);
    }
  }
  return leavingParticles;
}
}  // namespace LeavingParticleCollector
}  // namespace autopas
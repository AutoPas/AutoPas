/**
 * @file LeavingParticleCollector.h
 *
 * @author seckler
 * @date 13.08.2021
 */

#pragma once

#include "autopas/options/IteratorBehavior.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/inBox.h"

/**
 * Namespace to collect leaving particles from a container.
 */
namespace autopas::LeavingParticleCollector {
/**
 * Collects leaving particles and marks halo particles as dummy.
 * @note This function does not move or actually delete any particles!
 * @tparam ContainerType The type of the container.
 * @param container The container from which particles should be collected.
 * @return Returns a vector of leaving particles.
 */
template <class ContainerType>
std::vector<typename ContainerType::ParticleType> collectParticlesAndMarkNonOwnedAsDummy(ContainerType &container) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto upperBoundForMisplacement = container.getCutoff();
  const auto containerDimensions = container.getBoxMax() - container.getBoxMin();
  const auto lowerHaloCorner = container.getBoxMin() - upperBoundForMisplacement;
  const auto upperHaloCorner = container.getBoxMax() + upperBoundForMisplacement;
  // We need to extend these boundaries into the container because we need to find particles that are stored in the
  // container but have since moved out. Extending by skin is not enough since we are interested in particles that have
  // travelled beyond that and skin can be 0.
  const auto lowerInnerCorner = container.getBoxMin() + upperBoundForMisplacement;
  const auto upperInnerCorner = container.getBoxMax() - upperBoundForMisplacement;

  // volumes describing the halo regions on the six faces of the cuboid shaped container.
  const std::array<std::tuple<decltype(lowerHaloCorner), decltype(upperHaloCorner)>, 6> haloVolumes{{
      // left + lower edge
      {{lowerHaloCorner[0], lowerInnerCorner[1], lowerHaloCorner[2]},
       {lowerInnerCorner[0], upperInnerCorner[1], upperInnerCorner[2]}},
      // right + upper edge
      {{upperInnerCorner[0], lowerInnerCorner[1], lowerInnerCorner[2]},
       {upperHaloCorner[0], upperInnerCorner[1], upperHaloCorner[2]}},
      // top + left edge
      {{lowerHaloCorner[0], lowerInnerCorner[1], upperInnerCorner[2]},
       {upperInnerCorner[0], upperInnerCorner[1], upperHaloCorner[2]}},
      // bottom + right edge
      {{lowerInnerCorner[0], lowerInnerCorner[1], lowerHaloCorner[2]},
       {upperHaloCorner[0], upperInnerCorner[1], lowerInnerCorner[2]}},
      // front with all edges and corners
      {lowerHaloCorner, {upperHaloCorner[0], lowerInnerCorner[1], upperHaloCorner[2]}},
      // back with all edges and corners
      {{lowerHaloCorner[0], upperInnerCorner[1], lowerHaloCorner[2]}, upperHaloCorner},
  }};

  std::vector<typename ContainerType::ParticleType> leavingParticles{};
#ifdef AUTOPAS_OPENMP
#pragma omp declare reduction(vecMergeParticle : std::vector<typename ContainerType::ParticleType> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  // this has to be a `parallel` and not `parallel for` because we are actually interested in parallelizing the inner
  // loops without having a barrier between outer loop iterations.
  // The way this works is that all threads iterate the outer loop but due to the way the iterators are parallelized,
  // each thread only picks its iterations of the inner loops, hence nothing is done multiple times.
#pragma omp parallel reduction(vecMergeParticle : leavingParticles)
#endif
  for (const auto &[regionStart, regionEnd] : haloVolumes) {
    for (auto iter = container.getRegionIterator(regionStart, regionEnd, autopas::IteratorBehavior::ownedOrHalo);
         iter.isValid(); ++iter) {
      // Mark halo particles for deletion
      if (iter->isHalo()) {
        iter->setOwnershipState(OwnershipState::dummy);
        // Collect leaving particles and mark them for deletion
      } else if (not utils::inBox(iter->getR(), container.getBoxMin(), container.getBoxMax())) {
        leavingParticles.push_back(*iter);
        iter->setOwnershipState(OwnershipState::dummy);
      }
    }
  }
  return leavingParticles;
}
}  // namespace autopas::LeavingParticleCollector

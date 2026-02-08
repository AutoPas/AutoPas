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
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

/**
 * Namespace to collect leaving particles from a container.
 */
namespace autopas::LeavingParticleCollector {
/**
 * For a given container that supports `getCutoff()`, `getBoxMin()` and `getBoxMax()`,
 * calculates a set of 6 non-overlapping boxes that fully cover the halo.
 * These boxes extend into the container to catch halo particles which have moved inside the container.
 *
 * @note this function is public for testing purposes.
 *
 * @tparam ContainerType
 * @param container
 * @return Vector of halo boxes. vector<tuple<lower_corner, upper_corner>>
 */
template <class ContainerType>
std::array<std::tuple<const std::array<double, 3>, const std::array<double, 3>>, 6> calculateHaloVolumes(
    const ContainerType &container) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto upperBoundForMisplacement = container.getVerletSkin() / 2.0;
  const auto interactionLength = container.getInteractionLength();
  const auto &boxMin = container.getBoxMin();
  const auto &boxMax = container.getBoxMax();
  // Halo extends one interaction length away from the edge of the container. On top of that particles can drift.
  const auto lowerHaloCorner = boxMin - interactionLength - upperBoundForMisplacement;
  const auto upperHaloCorner = boxMax + interactionLength + upperBoundForMisplacement;
  // We need to extend these boundaries into the container because we need to include cells that still hold particles,
  // which now have left the cells region towards the outside of the container.
  // travelled beyond that and skin can be 0.
  const auto lowerInnerCorner = boxMin + upperBoundForMisplacement;
  const auto upperInnerCorner = boxMax - upperBoundForMisplacement;

  // volumes describing the halo regions on the six faces of the cuboid shaped container.
  // The total volume can be thought of as:
  // 6 faces
  // 8 corners
  // 12 edges
  // These 26 pieces are combined to six volumes, the minimal number to cover it with region iterators.
  const std::array<std::tuple<decltype(lowerHaloCorner), decltype(upperHaloCorner)>, 6> haloVolumes{{
      // left face + edge shared by left and bottom faces
      {{lowerHaloCorner[0], lowerInnerCorner[1], lowerHaloCorner[2]},
       {lowerInnerCorner[0], upperInnerCorner[1], upperInnerCorner[2]}},
      // right face + edge shared by right and upper faces
      {{upperInnerCorner[0], lowerInnerCorner[1], lowerInnerCorner[2]},
       {upperHaloCorner[0], upperInnerCorner[1], upperHaloCorner[2]}},
      // top face + edge shared by left and top faces
      {{lowerHaloCorner[0], lowerInnerCorner[1], upperInnerCorner[2]},
       {upperInnerCorner[0], upperInnerCorner[1], upperHaloCorner[2]}},
      // bottom face + edge shared by right and bottom faces
      {{lowerInnerCorner[0], lowerInnerCorner[1], lowerHaloCorner[2]},
       {upperHaloCorner[0], upperInnerCorner[1], lowerInnerCorner[2]}},
      // front face with the all adjacent four edges and four corners
      {lowerHaloCorner, {upperHaloCorner[0], lowerInnerCorner[1], upperHaloCorner[2]}},
      // back face with the all adjacent four edges and four corners
      {{lowerHaloCorner[0], upperInnerCorner[1], lowerHaloCorner[2]}, upperHaloCorner},
  }};

  return haloVolumes;
}

/**
 * Collects leaving particles and marks halo particles as dummy.
 * @note This function does not move or actually delete any particles!
 * @tparam ContainerType The type of the container.
 * @param container The container from which particles should be collected.
 * @return Returns a vector of leaving particles.
 */
template <class ContainerType>
std::vector<typename ContainerType::ParticleType> collectParticlesAndMarkNonOwnedAsDummy(ContainerType &container) {
  const auto haloVolumes = calculateHaloVolumes(container);

  std::vector<typename ContainerType::ParticleType> leavingParticles{};
  // Random number which is better than 0
  leavingParticles.reserve(10);

  // aliases for less getter accesses
  const auto &boxMin = container.getBoxMin();
  const auto &boxMax = container.getBoxMax();
  // custom openmp reduction to concatenate all local vectors to one at the end of a parallel region
  AUTOPAS_OPENMP(declare reduction(vecMergeParticle :                                                 \
                                   std::vector<typename ContainerType::ParticleType> :                \
                                       omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())))
  // this has to be a `parallel` and not `parallel for` because we are actually interested in parallelizing the inner
  // loops without having a barrier between outer loop iterations.
  // The way this works is that all threads iterate the outer loop but due to the way the iterators are parallelized,
  // each thread only picks its iterations of the inner loops, hence nothing is done multiple times.
  AUTOPAS_OPENMP(parallel reduction(vecMergeParticle : leavingParticles) num_threads(autopas_get_preferred_num_threads()))
  for (const auto &[regionStart, regionEnd] : haloVolumes) {
    for (auto iter = container.getRegionIterator(regionStart, regionEnd, autopas::IteratorBehavior::ownedOrHalo);
         iter.isValid(); ++iter) {
      // Mark halo particles for deletion
      if (iter->isHalo()) {
        iter->setOwnershipState(OwnershipState::dummy);
        // Collect leaving particles and mark them for deletion
      } else if (not utils::inBox(iter->getR(), boxMin, boxMax)) {
        leavingParticles.push_back(*iter);
        iter->setOwnershipState(OwnershipState::dummy);
      }
    }
  }
  return leavingParticles;
}
}  // namespace autopas::LeavingParticleCollector

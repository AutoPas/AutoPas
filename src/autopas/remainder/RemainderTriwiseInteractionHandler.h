/**
 * @file RemainderTriwiseInteractionHandler.h
 * @author muehlhaeusser
 * @date 12.02.1997
 */

#pragma once

#include <cmath>
#include <memory>
#include <mutex>
#include <numeric>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/cells/SortedCellView.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

/**
 * Handles triwise interactions involving particle buffers (particles not yet inserted into the main container).
 * This includes Buffer <-> Container and Buffer <-> Buffer interactions.
 *
 * @tparam Particle_T
 */
template <typename Particle_T>
class RemainderTriwiseInteractionHandler {
 public:
  /**
   * Constructor for RemainderTriwiseInteractionHandler
   * @param spatialLocks  Domain locks passed by LogicHandler
   */
  explicit RemainderTriwiseInteractionHandler(
      std::vector<std::vector<std::vector<std::unique_ptr<std::mutex>>>> &spatialLocks)
      : _spatialLocks(spatialLocks) {}

  /**
   * Performs the interactions ParticleContainer::computeInteractions() did not cover.
   *
   * These interactions are:
   *  - particleBuffer    <-> container
   *  - haloParticleBuffer -> container
   *  - particleBuffer    <-> particleBuffer
   *  - haloParticleBuffer -> particleBuffer
   *
   * @note Buffers need to have at least one (potentially itself empty) cell. They must not be an empty vector of cells.
   *
   * @tparam newton3
   * @tparam ContainerType Type of the particle container.
   * @tparam TriwiseFunctor Functor type to use for the interactions.
   * @param f
   * @param container Reference to the container. Preferably pass the container with the actual static type.
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class ContainerType, class TriwiseFunctor>
  void computeRemainderInteractions(TriwiseFunctor *f, ContainerType &container,
                                    std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                    std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers) {
    // Vector to collect pointers to all buffer particles
    std::vector<Particle_T *> bufferParticles;
    const auto numOwnedBufferParticles = collectBufferParticles(bufferParticles, particleBuffers, haloParticleBuffers);

    // The following part performs the main remainder traversal.

    // only activate time measurements if it will actually be logged
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    autopas::utils::Timer timerBufferBufferBuffer;
    autopas::utils::Timer timerBufferBufferContainer;
    autopas::utils::Timer timerBufferContainerContainer;
    timerBufferBufferBuffer.start();
#endif

    // Step 1: Triwise interactions of all particles in the buffers (owned and halo)
    remainderHelper3bBufferBufferBufferAoS(bufferParticles, numOwnedBufferParticles, f);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerBufferBufferBuffer.stop();
    timerBufferBufferContainer.start();
#endif

    // Step 2: Triwise interactions of 2 buffer particles with 1 container particle
    remainderHelper3bBufferBufferContainerAoS(bufferParticles, numOwnedBufferParticles, container, f);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerBufferBufferContainer.stop();
    timerBufferContainerContainer.start();
#endif

    // Step 3: Triwise interactions of 1 buffer particle and 2 container particles
    remainderHelper3bBufferContainerContainerAoS<newton3>(bufferParticles, container, f);
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerBufferContainerContainer.stop();
#endif

    AutoPasLog(TRACE, "Timer Buffer <-> Buffer <-> Buffer       : {}", timerBufferBufferBuffer.getTotalTime());
    AutoPasLog(TRACE, "Timer Buffer <-> Buffer <-> Container    : {}", timerBufferBufferContainer.getTotalTime());
    AutoPasLog(TRACE, "Timer Buffer <-> Container <-> Container : {}", timerBufferContainerContainer.getTotalTime());
  }

 private:
  /**
   * Reference to spatial locks to lock regional access to the container.
   */
  std::vector<std::vector<std::vector<std::unique_ptr<std::mutex>>>> &_spatialLocks;

  /**
   * Helper Method for computeRemainderInteractions. This method collects pointers to all owned halo particle buffers.
   *
   * @param bufferParticles Reference to a vector of particle pointers which will be filled with pointers to all buffer
   * particles (owned and halo).
   * @param particleBuffers Vector of particle buffers.
   * @param haloParticleBuffers Vector of halo particle buffers.
   * @return Number of owned buffer particles. The bufferParticles vector is two-way split, with the first half being
   * owned buffer particles and the second half being halo buffer particles.
   */
  size_t collectBufferParticles(std::vector<Particle_T *> &bufferParticles,
                                std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers) {
    // Reserve the needed amount
    auto cellToVec = [](auto &cell) -> std::vector<Particle_T> & { return cell._particles; };

    const size_t numOwnedBufferParticles =
        std::transform_reduce(particleBuffers.begin(), particleBuffers.end(), 0, std::plus<>(),
                              [&](auto &vec) { return cellToVec(vec).size(); });

    const size_t numHaloBufferParticles =
        std::transform_reduce(haloParticleBuffers.begin(), haloParticleBuffers.end(), 0, std::plus<>(),
                              [&](auto &vec) { return cellToVec(vec).size(); });

    bufferParticles.reserve(numOwnedBufferParticles + numHaloBufferParticles);

    // Collect owned buffer particles
    for (auto &buffer : particleBuffers) {
      for (auto &p : buffer._particles) {
        bufferParticles.push_back(&p);
      }
    }
    // Collect halo buffer particles
    for (auto &buffer : haloParticleBuffers) {
      for (auto &p : buffer._particles) {
        bufferParticles.push_back(&p);
      }
    }
    return numOwnedBufferParticles;
  }

  /**
   * Helper Method for computeRemainderInteractions. This method calculates all interactions between all triplets
   * within the buffer particles.
   *
   * @tparam TriwiseFunctor Functor type to use for the interactions.
   * @param bufferParticles Vector of Particle pointers containing pointers to all buffer particles (owned and halo).
   * @param numOwnedBufferParticles
   * @param f
   */
  template <class TriwiseFunctor>
  void remainderHelper3bBufferBufferBufferAoS(const std::vector<Particle_T *> &bufferParticles,
                                              const size_t numOwnedBufferParticles, TriwiseFunctor *f) {
    AUTOPAS_OPENMP(parallel for)
    for (auto i = 0; i < numOwnedBufferParticles; ++i) {
      Particle_T &p1 = *bufferParticles[i];

      for (auto j = 0; j < bufferParticles.size(); ++j) {
        if (i == j) continue;
        Particle_T &p2 = *bufferParticles[j];

        for (auto k = j + 1; k < bufferParticles.size(); ++k) {
          if (k == i) continue;
          Particle_T &p3 = *bufferParticles[k];

          f->AoSFunctor(p1, p2, p3, false);
        }
      }
    }
  }

  /**
   * Helper Method for computeRemainderInteractions. This method calculates all interactions between two buffer
   * particles and one particle from the container.
   *
   * @tparam ContainerType Type of the particle container.
   * @tparam TriwiseFunctor Functor type to use for the interactions.
   * @param bufferParticles Vector of Particle pointers containing pointers to all buffer particles (owned and halo).
   * @param numOwnedBufferParticles Number of owned particles. (First half of the bufferParticles).
   * @param container
   * @param f
   */
  template <class ContainerType, class TriwiseFunctor>
  void remainderHelper3bBufferBufferContainerAoS(const std::vector<Particle_T *> &bufferParticles,
                                                 size_t numOwnedBufferParticles, ContainerType &container,
                                                 TriwiseFunctor *f) {
    using utils::ArrayUtils::static_cast_copy_array;
    using namespace autopas::utils::ArrayMath::literals;

    const auto haloBoxMin = container.getBoxMin() - container.getInteractionLength();
    const auto interactionLengthInv = 1. / container.getInteractionLength();
    const double cutoff = container.getCutoff();

    AUTOPAS_OPENMP(parallel for)
    for (auto i = 0; i < bufferParticles.size(); ++i) {
      Particle_T &p1 = *bufferParticles[i];
      const auto min = p1.getR() - cutoff;
      const auto max = p1.getR() + cutoff;

      for (auto j = 0; j < bufferParticles.size(); ++j) {
        if (j == i) continue;
        Particle_T &p2 = *bufferParticles[j];
        container.forEachInRegion(
            [&](auto &p3) {
              const auto lockCoords = static_cast_copy_array<size_t>((p3.getR() - haloBoxMin) * interactionLengthInv);
              if (i < numOwnedBufferParticles) f->AoSFunctor(p1, p2, p3, false);
              if (!p3.isHalo() && i < j) {
                const std::lock_guard<std::mutex> lock(*_spatialLocks[lockCoords[0]][lockCoords[1]][lockCoords[2]]);
                f->AoSFunctor(p3, p1, p2, false);
              }
            },
            min, max, IteratorBehavior::ownedOrHalo);
      }
    }
  }

  /**
   * Helper Method for computeRemainderInteractions. This method calculates all interactions between one buffer
   * particle and two particles from the container.
   *
   * @tparam newton3
   * @tparam ContainerType Type of the particle container.
   * @tparam TriwiseFunctor Functor type to use for the interactions.
   * @param bufferParticles Vector of Particle pointers containing pointers to all buffer particles (owned and halo).
   * @param numOwnedBufferParticles Number of owned particles. (First half of the bufferParticles).
   * @param container
   * @param f
   */

  template <bool newton3, class ContainerType, class TriwiseFunctor>
  void remainderHelper3bBufferContainerContainerAoS(const std::vector<Particle_T *> &bufferParticles,
                                                    ContainerType &container, TriwiseFunctor *f) {
    using namespace autopas::utils::ArrayMath::literals;
    if (bufferParticles.empty()) {
      return;
    }
    const size_t numBufferParticles = bufferParticles.size();
    const double cutoff = container.getCutoff();
    const double cutoffSquared = cutoff * cutoff;

    // Step 1: Create a sorted cell view on all buffer particles
    const auto domainDiagonal = container.getBoxMax() - container.getBoxMin();
    ReferenceParticleCell<Particle_T> cellView;
    cellView.reserve(numBufferParticles);
    for (auto pRef : bufferParticles) {
      cellView.addParticleReference(pRef);
    }
    SortedCellView<ReferenceParticleCell<Particle_T>> sortedCellView(cellView, domainDiagonal);

    // Step 2: Find all valid particle pairs from the container for each of the buffer particles. This can be done in
    // parallel since all particles are accessed in read-only manner. The data structure is a vector, which holds a list
    // of all valid container particle pairs for each buffer particle.
    std::vector<std::vector<std::pair<Particle_T *, Particle_T *>>> bufferParticleNeighborLists(numBufferParticles);

    // Estimate the number of pairs per buffer particle to reserve memory
    const double searchVolume = 8. * cutoff * cutoff * cutoff;
    const double domainVolume = utils::ArrayMath::prod(container.getBoxMax() - container.getBoxMin());
    const double numParticlesPerSearchVolume =
        container.getNumberOfParticles(IteratorBehavior::owned) * (searchVolume / domainVolume);
    const auto estimatedPairsPerParticle =
        static_cast<size_t>(numParticlesPerSearchVolume * numParticlesPerSearchVolume / 2);

  AUTOPAS_OPENMP(parallel for schedule(static))
  for (size_t i = 0; i < numBufferParticles; i++) {
    bufferParticleNeighborLists[i].reserve(estimatedPairsPerParticle);
    auto &[pProjection, p1Ptr] = sortedCellView._particles[i];
    const auto pos1 = p1Ptr->getR();
    const auto boxMin = pos1 - cutoff;
    const auto boxMax = pos1 + cutoff;

    auto p2Iter = container.getRegionIterator(
        boxMin, boxMax, IteratorBehavior::ownedOrHalo | IteratorBehavior::forceSequential, std::nullopt);
    for (; p2Iter.isValid(); ++p2Iter) {
      const auto pos2 = p2Iter->getR();
      const auto p1p2Dist = pos2 - pos1;
      if (utils::ArrayMath::dot(p1p2Dist, p1p2Dist) > cutoffSquared) {
        continue;
      }
      auto p3Iter = p2Iter;
      ++p3Iter;
      for (; p3Iter.isValid(); ++p3Iter) {
        const auto pos3 = p3Iter->getR();
        const auto p2p3Dist = pos3 - pos2;
        if (utils::ArrayMath::dot(p2p3Dist, p2p3Dist) > cutoffSquared) {
          continue;
        }
        const auto p1p3Dist = pos3 - pos1;
        if (utils::ArrayMath::dot(p1p3Dist, p1p3Dist) > cutoffSquared) {
          continue;
        }
        // All distances are below the cutoff, so append to the neighbor list.
        bufferParticleNeighborLists[i].emplace_back(&(*p2Iter), &(*p3Iter));
      }
    }
  }

  // Step 3: Bin buffer particles into slices along the domain diagonal of size cutoff
  // Every "bin" contains pointers to the respective buffer particles alongside their index in the sortedCellView to
  // simplify mapping back later.
  std::vector<std::vector<std::pair<size_t, Particle_T *>>> binnedBufferParticles;
  // const auto estimatedSizePerBin = numBufferParticles / (utils::ArrayMath::L2Norm(domainDiagonal) / cutoff);
  binnedBufferParticles.emplace_back();
  double binStart = sortedCellView._particles[0].first;  // Projection of the first particle
  for (size_t i = 0; i < numBufferParticles; ++i) {
    auto &[pProjection, pPtr] = sortedCellView._particles[i];

    // The slice thickness is at least the cutoff. The next slice starts with the position of the next particle.
    if (pProjection - binStart > cutoff) {
      binStart = pProjection;
      binnedBufferParticles.emplace_back();
    }
    binnedBufferParticles.back().emplace_back(i, pPtr);
  }

  // Step 4: Iterate over the binned buffer particles in a 3-colored parallel manner to avoid race conditions.
  for (auto color = 0; color < 3; color++) {
    AUTOPAS_OPENMP(parallel for schedule(dynamic))
    for (size_t i = 0 + color; i < binnedBufferParticles.size(); i += 3) {
      for (auto &[p1Index, p1Ptr] : binnedBufferParticles[i]) {
        for (auto &[p2Ptr, p3Ptr] : bufferParticleNeighborLists[p1Index]) {
          if constexpr (newton3) {
            f->AoSFunctor(*p1Ptr, *p2Ptr, *p3Ptr, true);
          } else {
            if (p1Ptr->isOwned()) {
              f->AoSFunctor(*p1Ptr, *p2Ptr, *p3Ptr, false);
            }
            if (p2Ptr->isOwned()) {
              f->AoSFunctor(*p2Ptr, *p1Ptr, *p3Ptr, false);
            }
            if (p3Ptr->isOwned()) {
              f->AoSFunctor(*p3Ptr, *p1Ptr, *p2Ptr, false);
            }
          }
        }
      }
    }
  }
  }
};
}  // namespace autopas
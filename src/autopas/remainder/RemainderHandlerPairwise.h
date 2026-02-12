/**
 * @file RemainderHandlerPairwise.h
 * @author muehlhaeusser
 * @date 12.02.1997
 */

#pragma once

#include <array>
#include <cmath>
#include <memory>
#include <mutex>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/checkFunctorType.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

/**
 * Handles pairwise interactions involving particle buffers (particles not yet inserted into the main container).
 * This includes Buffer <-> Container and Buffer <-> Buffer interactions.
 *
 * @tparam Particle_T
 */
template <typename Particle_T>
class RemainderHandlerPairwise {
 public:
  explicit RemainderHandlerPairwise(std::vector<std::vector<std::vector<std::unique_ptr<std::mutex>>>> &spatialLocks)
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
   * @note Buffers need to have at least one (empty) cell. They must not be empty.
   *
   * @tparam newton3
   * @tparam ContainerType Type of the particle container.
   * @tparam PairwiseFunctor
   * @param f
   * @param container Reference the container. Preferably pass the container with the actual static type.
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   * @param useSoA Use SoAFunctor calls instead of AoSFunctor calls wherever it makes sense.
   */
  template <bool newton3, class ContainerType, class PairwiseFunctor>
  void computeRemainderInteractions(PairwiseFunctor *f, ContainerType &container,
                                    std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                    std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers, bool useSoA) {
    // Balance buffers. This makes processing them with static scheduling quite efficient.
    // Also, if particles were not inserted in parallel, this enables us to process them in parallel now.
    // Cost is at max O(2N) worst O(N) per buffer collection and negligible compared to interacting them.
    auto cellToVec = [](auto &cell) -> std::vector<Particle_T> & { return cell._particles; };
    utils::ArrayUtils::balanceVectors(particleBuffers, cellToVec);
    utils::ArrayUtils::balanceVectors(haloParticleBuffers, cellToVec);

    // The following part performs the main remainder traversal. The actual calculation is done in 4 steps carried out
    // in three helper functions.

    // only activate time measurements if it will actually be logged
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    autopas::utils::Timer timerBufferContainer, timerPBufferPBuffer, timerPBufferHBuffer, timerBufferSoAConversion;
    timerBufferContainer.start();
#endif
    // steps 1 & 2.
    // particleBuffer with all particles close in container
    // and haloParticleBuffer with owned, close particles in container.
    // This is always AoS-based because the container particles are found with region iterators,
    // which don't have an SoA interface.
    remainderHelperBufferContainerAoS<newton3>(f, container, particleBuffers, haloParticleBuffers);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerBufferContainer.stop();
    timerBufferSoAConversion.start();
#endif

    if (useSoA) {
      // All (halo-)buffer interactions shall happen vectorized, hence, load all buffer data into SoAs
      for (auto &buffer : particleBuffers) {
        f->SoALoader(buffer, buffer._particleSoABuffer, 0, false);
      }
      for (auto &buffer : haloParticleBuffers) {
        f->SoALoader(buffer, buffer._particleSoABuffer, 0, false);
      }
    }
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerBufferSoAConversion.stop();
    timerPBufferPBuffer.start();
#endif

    // step 3. particleBuffer with itself and all other buffers
    remainderHelperBufferBuffer<newton3>(f, particleBuffers, useSoA);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerPBufferPBuffer.stop();
    timerPBufferHBuffer.start();
#endif

    // step 4. particleBuffer with haloParticleBuffer
    remainderHelperBufferHaloBuffer(f, particleBuffers, haloParticleBuffers, useSoA);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerPBufferHBuffer.stop();
    timerBufferSoAConversion.start();
#endif

    // unpack particle SoAs. Halo data is not interesting
    if (useSoA) {
      for (auto &buffer : particleBuffers) f->SoAExtractor(buffer, buffer._particleSoABuffer, 0);
    }

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerBufferSoAConversion.stop();
#endif

    AutoPasLog(TRACE, "Timer Buffers <-> Container  (1+2): {}", timerBufferContainer.getTotalTime());
    AutoPasLog(TRACE, "Timer PBuffers<-> PBuffer    (  3): {}", timerPBufferPBuffer.getTotalTime());
    AutoPasLog(TRACE, "Timer PBuffers<-> HBuffer    (  4): {}", timerPBufferHBuffer.getTotalTime());
    AutoPasLog(TRACE, "Timer Load and extract SoA buffers: {}", timerBufferSoAConversion.getTotalTime());

    // Note: haloParticleBuffer with itself is NOT needed, as interactions between halo particles are unneeded!
  }

 private:
  /**
   * Reference to spatial locks to lock regional access to the container.
   */
  std::vector<std::vector<std::vector<std::unique_ptr<std::mutex>>>> &_spatialLocks;

  /**
   * Helper Method for computeRemainderInteractions.
   * This method calculates all interactions between buffers and containers.
   *
   * @tparam newton3
   * @tparam ContainerType Type of the particle container.
   * @tparam PairwiseFunctor
   * @param f
   * @param container Reference the container. Preferably pass the container with the actual static type.
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class ContainerType, class PairwiseFunctor>
  void remainderHelperBufferContainerAoS(PairwiseFunctor *f, ContainerType &container,
                                         std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                         std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers) {
    using utils::ArrayUtils::static_cast_copy_array;
    using namespace autopas::utils::ArrayMath::literals;

    // Bunch of shorthands
    const auto cutoff = container.getCutoff();
    const auto interactionLength = container.getInteractionLength();
    const auto haloBoxMin = container.getBoxMin() - interactionLength;
    const auto totalBoxLengthInv = 1. / (container.getBoxMax() + interactionLength - haloBoxMin);
    const std::array<size_t, 3> spacialLocksPerDim{_spatialLocks.size(), _spatialLocks[0].size(),
                                                   _spatialLocks[0][0].size()};

    // Helper function to obtain the lock responsible for a given position.
    // Implemented as lambda because it can reuse a lot of information that is known in this context.
    const auto getSpacialLock = [&](const std::array<double, 3> &pos) -> std::mutex & {
      const auto posDistFromLowerCorner = pos - haloBoxMin;
      const auto relativePos = posDistFromLowerCorner * totalBoxLengthInv;
      // Lock coordinates are the position scaled by the number of locks
      const auto lockCoords =
          static_cast_copy_array<size_t>(static_cast_copy_array<double>(spacialLocksPerDim) * relativePos);
      return *_spatialLocks[lockCoords[0]][lockCoords[1]][lockCoords[2]];
    };

    // one halo and particle buffer pair per thread
    AUTOPAS_OPENMP(parallel for schedule(static, 1) default(shared))
    for (int bufferId = 0; bufferId < particleBuffers.size(); ++bufferId) {
      auto &particleBuffer = particleBuffers[bufferId];
      auto &haloParticleBuffer = haloParticleBuffers[bufferId];

      // 1. particleBuffer with all close particles in container
      for (auto &&p1 : particleBuffer) {
        const auto min = p1.getR() - cutoff;
        const auto max = p1.getR() + cutoff;
        container.forEachInRegion(
            [&](auto &p2) {
              if constexpr (newton3) {
                const std::lock_guard<std::mutex> lock(getSpacialLock(p2.getR()));
                f->AoSFunctor(p1, p2, true);
              } else {
                f->AoSFunctor(p1, p2, false);
                // no need to calculate force enacted on a halo
                if (not p2.isHalo()) {
                  const std::lock_guard<std::mutex> lock(getSpacialLock(p2.getR()));
                  f->AoSFunctor(p2, p1, false);
                }
              }
            },
            min, max, IteratorBehavior::ownedOrHalo);
      }

      // 2. haloParticleBuffer with owned, close particles in container
      for (auto &&p1halo : haloParticleBuffer) {
        const auto min = p1halo.getR() - cutoff;
        const auto max = p1halo.getR() + cutoff;
        container.forEachInRegion(
            [&](auto &p2) {
              // No need to apply anything to p1halo
              //   -> AoSFunctor(p1, p2, false) not needed as it neither adds force nor Upot (potential energy)
              //   -> newton3 argument needed for correct globals
              const std::lock_guard<std::mutex> lock(getSpacialLock(p2.getR()));
              f->AoSFunctor(p2, p1halo, newton3);
            },
            min, max, IteratorBehavior::owned);
      }
    }
  }

  /**
   * Helper Method for computeRemainderInteractions.
   * This method calculates all interactions between buffers and buffers
   * @tparam newton3
   * @tparam PairwiseFunctor
   * @param f
   * @param particleBuffers Vector of particle buffers. These particles' force vectors will be updated.
   * @param useSoA Use SoA based interactions instead of AoS.
   */
  template <bool newton3, class PairwiseFunctor>
  void remainderHelperBufferBuffer(PairwiseFunctor *f, std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                   bool useSoA) {
    if (useSoA)
      remainderHelperBufferBufferSoA<newton3>(f, particleBuffers);
    else
      remainderHelperBufferBufferAoS<newton3>(f, particleBuffers);
  }

  /**
   * AoS implementation for remainderHelperBufferBuffer().
   * @tparam newton3
   * @tparam PairwiseFunctor
   * @param f
   * @param particleBuffers
   */
  template <bool newton3, class PairwiseFunctor>
  void remainderHelperBufferBufferAoS(PairwiseFunctor *f, std::vector<FullParticleCell<Particle_T>> &particleBuffers) {
    // For all interactions between different buffers we turn newton3 always off,
    // which ensures that only one thread at a time is writing to a buffer.
    // This saves expensive locks.

    // We can not use collapse here without locks, otherwise races would occur.
    AUTOPAS_OPENMP(parallel for)
    for (size_t bufferIdxI = 0; bufferIdxI < particleBuffers.size(); ++bufferIdxI) {
      for (size_t bufferIdxJOffset = 0; bufferIdxJOffset < particleBuffers.size(); ++bufferIdxJOffset) {
        // Let each bufferI use a different starting point for bufferJ to minimize false sharing
        const auto bufferIdxJ = (bufferIdxI + bufferIdxJOffset) % particleBuffers.size();

        // interact the two buffers
        if (bufferIdxI == bufferIdxJ) {
          // CASE Same buffer
          // Only use Newton3 if it is allowed, and we are working on only one buffer. This avoids data races.
          const bool useNewton3 = newton3;
          auto &bufferRef = particleBuffers[bufferIdxI];
          const auto bufferSize = bufferRef.size();
          for (auto i = 0; i < bufferSize; ++i) {
            auto &p1 = bufferRef[i];
            // If Newton3 is disabled run over the whole buffer, otherwise only what is ahead
            for (auto j = useNewton3 ? i + 1 : 0; j < bufferSize; ++j) {
              if (i == j) {
                continue;
              }
              auto &p2 = bufferRef[j];
              f->AoSFunctor(p1, p2, useNewton3);
            }
          }
        } else {
          // CASE: Two buffers
          for (auto &p1 : particleBuffers[bufferIdxI]) {
            for (auto &p2 : particleBuffers[bufferIdxJ]) {
              f->AoSFunctor(p1, p2, false);
            }
          }
        }
      }
    }
  }

  /**
   * SoA implementation for remainderHelperBufferBuffer().
   * @tparam newton3
   * @tparam PairwiseFunctor
   * @param f
   * @param particleBuffers
   */
  template <bool newton3, class PairwiseFunctor>
  void remainderHelperBufferBufferSoA(PairwiseFunctor *f, std::vector<FullParticleCell<Particle_T>> &particleBuffers) {
    // we can not use collapse here without locks, otherwise races would occur.
    AUTOPAS_OPENMP(parallel for)
    for (size_t i = 0; i < particleBuffers.size(); ++i) {
      for (size_t jj = 0; jj < particleBuffers.size(); ++jj) {
        auto *particleBufferSoAA = &particleBuffers[i]._particleSoABuffer;
        // instead of starting every (parallel) iteration i at j == 0 offset them by i to minimize waiting times at
        // locks
        const auto j = (i + jj) % particleBuffers.size();
        if (i == j) {
          // For buffer interactions where bufferA == bufferB we can always enable newton3 if it is allowed.
          f->SoAFunctorSingle(*particleBufferSoAA, newton3);
        } else {
          // For all interactions between different buffers we turn newton3 always off,
          // which ensures that only one thread at a time is writing to a buffer. This saves expensive locks.
          auto *particleBufferSoAB = &particleBuffers[j]._particleSoABuffer;
          f->SoAFunctorPair(*particleBufferSoAA, *particleBufferSoAB, false);
        }
      }
    }
  }

  /**
   * Helper Method for computeRemainderInteractions.
   * This method calculates all interactions between buffers and halo buffers.
   *
   * @note This function never uses Newton3 because we do not need to calculate the effect on the halo particles.
   *
   * @tparam PairwiseFunctor
   * @param f
   * @param particleBuffers Vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   * @param useSoA Use SoA based interactions instead of AoS.
   */
  template <class PairwiseFunctor>
  void remainderHelperBufferHaloBuffer(PairwiseFunctor *f, std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                       std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers, bool useSoA) {
    if (useSoA)
      remainderHelperBufferHaloBufferSoA(f, particleBuffers, haloParticleBuffers);
    else
      remainderHelperBufferHaloBufferAoS(f, particleBuffers, haloParticleBuffers);
  }

  /**
   * AoS implementation for remainderHelperBufferHaloBuffer()
   * @tparam PairwiseFunctor
   * @param f
   * @param particleBuffers
   * @param haloParticleBuffers
   */
  template <class PairwiseFunctor>
  void remainderHelperBufferHaloBufferAoS(PairwiseFunctor *f,
                                          std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                          std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers) {
    // Here, phase / color based parallelism turned out to be more efficient than tasks
    AUTOPAS_OPENMP(parallel)
    for (int interactionOffset = 0; interactionOffset < haloParticleBuffers.size(); ++interactionOffset) {
      AUTOPAS_OPENMP(for)
      for (size_t i = 0; i < particleBuffers.size(); ++i) {
        auto &particleBuffer = particleBuffers[i];
        auto &haloBuffer = haloParticleBuffers[(i + interactionOffset) % haloParticleBuffers.size()];

        for (auto &p1 : particleBuffer) {
          for (auto &p2 : haloBuffer) {
            f->AoSFunctor(p1, p2, false);
          }
        }
      }
    }
  }

  /**
   * SoA implementation for remainderHelperBufferHaloBuffer()
   * @tparam PairwiseFunctor
   * @param f
   * @param particleBuffers
   * @param haloParticleBuffers
   */
  template <class PairwiseFunctor>
  void remainderHelperBufferHaloBufferSoA(PairwiseFunctor *f,
                                          std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                          std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers) {
    // Here, phase / color based parallelism turned out to be more efficient than tasks
    AUTOPAS_OPENMP(parallel)
    for (int interactionOffset = 0; interactionOffset < haloParticleBuffers.size(); ++interactionOffset) {
      AUTOPAS_OPENMP(for)
      for (size_t i = 0; i < particleBuffers.size(); ++i) {
        auto &particleBufferSoA = particleBuffers[i]._particleSoABuffer;
        auto &haloBufferSoA =
            haloParticleBuffers[(i + interactionOffset) % haloParticleBuffers.size()]._particleSoABuffer;
        f->SoAFunctorPair(particleBufferSoA, haloBufferSoA, false);
      }
    }
  }
};
}  // namespace autopas
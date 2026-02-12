/**
 * @file RemainderHandler.h
 * @author muehlhaeusser
 * @date 12.02.1997
 */

#pragma once

#include <array>
#include <cmath>
#include <memory>
#include <mutex>
#include <numeric>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/checkFunctorType.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

/**
 * Handles interactions involving particle buffers (particles not yet inserted into the main container).
 * This includes Buffer <-> Container and Buffer <-> Buffer interactions.
 *
 * @tparam Particle_T
 */
template <typename Particle_T>
class RemainderHandler {
 public:
  /**
   * Constructor.
   * Initializes buffer locks based on thread count.
   */
  RemainderHandler() : _bufferLocks(std::max(2, autopas_get_max_threads())) {
    for (auto &lockPtr : _bufferLocks) {
      lockPtr = std::make_unique<std::mutex>();
    }
  }

  /**
   * Initialize or update the spatial locks used during the remainder traversal.
   * If the locks are already initialized but the container size changed, surplus locks will
   * be deleted, new locks are allocated and locks that are still necessary are reused.
   *
   * @note Should the generated number of locks exceed 1e6, the number of locks is reduced so that
   * the number of locks per dimensions is proportional to the domain side lengths and smaller than the limit.
   *
   * @param boxLength Size per dimension for the box that should be covered by locks (should include halo).
   * @param interactionLengthInv Inverse of the side length of the virtual boxes one lock is responsible for.
   *
   */
  void initSpatialLocks(const std::array<double, 3> &boxLength, double interactionLengthInv) {
    using namespace autopas::utils::ArrayMath::literals;
    using utils::ArrayMath::ceil;
    using utils::ArrayUtils::static_cast_copy_array;

    // The maximum number of spatial locks is capped at 1e6.
    // This limit is chosen more or less arbitrary. It is big enough so that our regular MD simulations
    // fall well within it and small enough so that no memory issues arise.
    // There were no rigorous tests for an optimal number of locks.
    // Without this cap, very large domains (or tiny cutoffs) would generate an insane number of locks,
    // that could blow up the memory.
    constexpr size_t maxNumSpacialLocks{1000000};

    // One lock per interaction length or less if this would generate too many.
    const std::array<size_t, 3> locksPerDim = [&]() {
      // First naively calculate the number of locks if we simply take the desired cell length.
      // Ceil because both decisions are possible, and we are generous gods.
      const std::array<size_t, 3> locksPerDimNaive =
          static_cast_copy_array<size_t>(ceil(boxLength * interactionLengthInv));
      const auto totalLocksNaive =
          std::accumulate(locksPerDimNaive.begin(), locksPerDimNaive.end(), 1ul, std::multiplies<>());
      // If the number of locks is within the limits everything is fine and we can return.
      if (totalLocksNaive <= maxNumSpacialLocks) {
        return locksPerDimNaive;
      } else {
        // If the number of locks grows too large, calculate the locks per dimension proportionally to the side lengths.
        // Calculate side length relative to dimension 0.
        const std::array<double, 3> boxSideProportions = {
            1.,
            boxLength[0] / boxLength[1],
            boxLength[0] / boxLength[2],
        };
        // With this, calculate the number of locks the first dimension should receive.
        const auto prodProportions =
            std::accumulate(boxSideProportions.begin(), boxSideProportions.end(), 1., std::multiplies<>());
        // Needs floor, otherwise we exceed the limit.
        const auto locksInFirstDimFloat = std::floor(std::cbrt(maxNumSpacialLocks * prodProportions));
        // From this and the proportions relative to the first dimension, we can calculate the remaining number of locks
        return std::array<size_t, 3>{
            static_cast<size_t>(locksInFirstDimFloat),
            static_cast<size_t>(locksInFirstDimFloat / boxSideProportions[1]),
            static_cast<size_t>(locksInFirstDimFloat / boxSideProportions[2]),
        };
      }
    }();

    _spatialLocks.resize(locksPerDim[0]);
    for (auto &lockVecVec : _spatialLocks) {
      lockVecVec.resize(locksPerDim[1]);
      for (auto &lockVec : lockVecVec) {
        lockVec.resize(locksPerDim[2]);
        for (auto &lockPtr : lockVec) {
          if (not lockPtr) {
            lockPtr = std::make_unique<std::mutex>();
          }
        }
      }
    }
  }

  /**
   * Select the right Remainder function depending on the interaction type and newton3 setting.
   *
   * @tparam Functor
   * @tparam ContainerType
   * @param functor
   * @param container
   * @param particleBuffers
   * @param haloParticleBuffers
   * @param newton3
   * @param useSoA Use SoA functor calls where it is feasible. This still leaves some interactions with the AoS functor
   * but the switch is useful to disable SoA calls if a functor doesn't support them at all.
   * @return
   */
  template <class Functor, class ContainerType>
  void computeRemainderInteractions(Functor &functor, ContainerType &container,
                                    std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                    std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers, bool newton3,
                                    bool useSoA) {
    if constexpr (utils::isPairwiseFunctor<Functor>()) {
      if (newton3) {
        computeRemainderInteractions2B<true>(&functor, container, particleBuffers, haloParticleBuffers, useSoA);
      } else {
        computeRemainderInteractions2B<false>(&functor, container, particleBuffers, haloParticleBuffers, useSoA);
      }
    } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
      if (newton3) {
        computeRemainderInteractions3B<true>(&functor, container, particleBuffers, haloParticleBuffers);
      } else {
        computeRemainderInteractions3B<false>(&functor, container, particleBuffers, haloParticleBuffers);
      }
    }
  }

 private:
  /// Locks for regions in the domain.
  std::vector<std::vector<std::vector<std::unique_ptr<std::mutex>>>> _spatialLocks;
  /// Locks for the particle buffers.
  std::vector<std::unique_ptr<std::mutex>> _bufferLocks;

  // ------------------------- 2-Body Interactions -------------------------

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
  void computeRemainderInteractions2B(PairwiseFunctor *f, ContainerType &container,
                                      std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                                      std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers, bool useSoA) {
    // Sanity check. If this is violated feel free to add some logic here that adapts the number of locks.
    if (_bufferLocks.size() < particleBuffers.size()) {
      utils::ExceptionHandler::exception("RemainderHandler: Not enough locks for non-halo buffers!");
    }

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

  /**
   * Helper Method for computeRemainderInteractions2B.
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
    using autopas::utils::ArrayUtils::static_cast_copy_array;
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
   * Helper Method for computeRemainderInteractions2B.
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
   * Helper Method for computeRemainderInteractions2B.
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

  // ------------------------- 3-Body Interactions -------------------------

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
   * @tparam TriwiseFunctor Functor type to use for the interactions.
   * @param f
   * @param container Reference to the container. Preferably pass the container with the actual static type.
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class ContainerType, class TriwiseFunctor>
  void computeRemainderInteractions3B(TriwiseFunctor *f, ContainerType &container,
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
    remainderHelper3bBufferContainerContainerAoS<newton3>(bufferParticles, numOwnedBufferParticles, container, f);
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
    timerBufferContainerContainer.stop();
#endif

    AutoPasLog(TRACE, "Timer Buffer <-> Buffer <-> Buffer       : {}", timerBufferBufferBuffer.getTotalTime());
    AutoPasLog(TRACE, "Timer Buffer <-> Buffer <-> Container    : {}", timerBufferBufferContainer.getTotalTime());
    AutoPasLog(TRACE, "Timer Buffer <-> Container <-> Container : {}", timerBufferContainerContainer.getTotalTime());
  }

  /**
   * Helper Method for computeRemainderInteractions3B. This method collects pointers to all owned halo particle buffers.
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
   * Helper Method for computeRemainderInteractions3B. This method calculates all interactions between all triplets
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
   * Helper Method for computeRemainderInteractions3B. This method calculates all interactions between two buffer
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
   * Helper Method for computeRemainderInteractions3B. This method calculates all interactions between one buffer
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
                                                    size_t numOwnedBufferParticles, ContainerType &container,
                                                    TriwiseFunctor *f) {
    // Todo: parallelize without race conditions - https://github.com/AutoPas/AutoPas/issues/904
    using utils::ArrayUtils::static_cast_copy_array;
    using namespace autopas::utils::ArrayMath::literals;

    const double cutoff = container.getCutoff();
    for (auto i = 0; i < bufferParticles.size(); ++i) {
      Particle_T &p1 = *bufferParticles[i];
      const auto boxMin = p1.getR() - cutoff;
      const auto boxMax = p1.getR() + cutoff;

      auto p2Iter = container.getRegionIterator(
          boxMin, boxMax, IteratorBehavior::ownedOrHalo | IteratorBehavior::forceSequential, nullptr);
      for (; p2Iter.isValid(); ++p2Iter) {
        Particle_T &p2 = *p2Iter;

        auto p3Iter = p2Iter;
        ++p3Iter;

        for (; p3Iter.isValid(); ++p3Iter) {
          Particle_T &p3 = *p3Iter;

          if constexpr (newton3) {
            f->AoSFunctor(p1, p2, p3, true);
          } else {
            if (i < numOwnedBufferParticles) {
              f->AoSFunctor(p1, p2, p3, false);
            }
            if (p2.isOwned()) {
              f->AoSFunctor(p2, p1, p3, false);
            }
            if (p3.isOwned()) {
              f->AoSFunctor(p3, p1, p2, false);
            }
          }
        }
      }
    }
  }
};
}  // namespace autopas
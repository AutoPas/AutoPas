/**
 * @file VerletListHelpers.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <atomic>
#include <unordered_map>
#include <vector>

#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
namespace autopas {

/**
 * Class of helpers for the VerletLists class.
 * @tparam Particle_T
 */
template <class Particle_T>
class VerletListHelpers {
 public:
  /**
   * Flat Compressed-Row-Storage (CRS) neighbor list.
   *
   * For particle i:
   *   - neighbor count: offsets[i+1] - offsets[i]
   *   - neighbor slice: indices.data() + offsets[i]
   */
  struct NeighborListCRS {
    std::vector<size_t> offsets;                            ///< size N+1
    std::vector<size_t, AlignedAllocator<size_t>> indices;  ///< flat neighbor indices

    /**
     * Number of particles tracked by this list.
     * @return number of particles
     */
    [[nodiscard]] size_t size() const { return offsets.empty() ? 0u : offsets.size() - 1u; }
    /** Number of neighbors of particle i.
     * @param i particle index
     * @return number of neighbors
     */
    [[nodiscard]] size_t count(const size_t i) const { return offsets[i + 1] - offsets[i]; }
    /** Pointer to the first neighbor index of particle i (const).
     * @param i particle index
     * @return pointer to the first neighbor index
     */
    [[nodiscard]] const size_t *begin(const size_t i) const { return indices.data() + offsets[i]; }
    /** Pointer to the first neighbor index of particle i (mutable).
     * @param i particle index
     * @return pointer to the first neighbor index
     */
    [[nodiscard]] size_t *begin(const size_t i) { return indices.data() + offsets[i]; }
    /** Returns a span of the neighbors of particle i.
     * @param i particle index
     * @return span of neighbor indices
     */
    [[nodiscard]] std::span<const size_t> getNeighbors(const size_t i) const { return {begin(i), count(i)}; }
  };

  /**
   * Policy for managing neighbor lists using the Compressed-Row-Storage (CRS) format.
   *
   * The neighbor lists are "temporary" neighbor lists that get compressed to CRS neighbor lists afterward.
   * This policy is required by the InteractionListGeneratorFunctor.
   */
  class CRSNeighborListPolicy {
   public:
    /**
     * Constructor.
     * @param neighborLists temporary neighbor list vectors
     * @param particleToIndex map to connect particle pointers to their dense SoA indices
     */
    CRSNeighborListPolicy(std::vector<std::vector<size_t>> &neighborLists,
                          std::unordered_map<const Particle_T *, size_t> &particleToIndex)
        : _neighborLists(neighborLists), _particleToIndex(particleToIndex) {}

    /**
     * Adds particle j to the neighbor list of particle i.
     * @param i the first particle
     * @param j the neighbor particle
     */
    void add(Particle_T *i, Particle_T *j) { _neighborLists[_particleToIndex[i]].push_back(_particleToIndex[j]); }

   private:
    std::vector<std::vector<size_t>> &_neighborLists;
    std::unordered_map<const Particle_T *, size_t> &_particleToIndex;
  };

  /**
   * Pass-1 functor for the two-pass CRS build.
   *
   * For every interacting pair (i, j) it increments a per-particle neighbor counter. The counters are then
   * prefix-summed to produce the CRS offsets before pass 2 runs.
   *
   * Thread safety: each counter is a std::atomic<size_t> padded to a full cache line (64 bytes) so that adjacent
   * particles never share a cache line. This eliminates false-sharing stalls that would otherwise dominate on
   * multi-socket / many-core systems.
   */
  class VerletListCounterFunctor : public PairwiseFunctor<Particle_T, VerletListCounterFunctor> {
   public:
    /**
     * Structure of the SoAs defined by the particle.
     */
    using SoAArraysType = typename Particle_T::SoAArraysType;

    /**
     * One atomic counter per particle, each padded to its own cache line.
     */
    struct alignas(64) PaddedAtomic {
      /**
       * The actual counter value.
       */
      std::atomic<size_t> value{0};
    };

    /**
     * @param counts        Per-particle neighbor counters (size N).
     * @param particleToIndex Map from particle pointer to its dense SoA index.
     * @param interactionLength cutoff + skin.
     */
    VerletListCounterFunctor(std::vector<PaddedAtomic> &counts,
                             const std::unordered_map<const Particle_T *, size_t> &particleToIndex,
                             double interactionLength)
        : PairwiseFunctor<Particle_T, VerletListCounterFunctor>(interactionLength),
          _counts(counts),
          _particleToIndex(particleToIndex),
          _interactionLengthSquared(interactionLength * interactionLength) {}

    std::string getName() override { return "VerletListCounterFunctor"; }
    bool isRelevantForTuning() override { return false; }
    bool allowsNewton3() override { return true; }
    bool allowsNonNewton3() override { return true; }

    /**
     * @copydoc autopas::PairwiseFunctor::AoSFunctor()
     */
    void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) override {
      using namespace autopas::utils::ArrayMath::literals;
      if (i.isDummy() or j.isDummy()) return;
      auto dist = i.getR() - j.getR();
      if (utils::ArrayMath::dot(dist, dist) < _interactionLengthSquared) {
        _counts[_particleToIndex.at(&i)].value.fetch_add(1, std::memory_order_relaxed);
        // newton3 ignored: AoSFunctor(j, i) is also called for newton3=false.
      }
    }

    /**
     * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
     */
    void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
      if (soa.size() == 0) return;
      auto **const __restrict ptrptr = soa.template begin<Particle_T::AttributeNames::ptr>();
      const double *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
      const double *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
      const double *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
      const size_t n = soa.size();
      for (size_t i = 0; i < n; ++i) {
        const size_t iIdx = _particleToIndex.at(ptrptr[i]);
        size_t localCount = 0;
        for (size_t j = i + 1; j < n; ++j) {
          const double dx = xptr[i] - xptr[j], dy = yptr[i] - yptr[j], dz = zptr[i] - zptr[j];
          if (dx * dx + dy * dy + dz * dz < _interactionLengthSquared) {
            ++localCount;
            if (not newton3) {
              _counts[_particleToIndex.at(ptrptr[j])].value.fetch_add(1, std::memory_order_relaxed);
            }
          }
        }
        _counts[iIdx].value.fetch_add(localCount, std::memory_order_relaxed);
      }
    }

    /**
     * @param soa1 SoA of first cell
     * @param soa2 SoA of second cell
     */
    void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool /*newton3*/) override {
      if (soa1.size() == 0 || soa2.size() == 0) return;
      auto **const __restrict ptr1ptr = soa1.template begin<Particle_T::AttributeNames::ptr>();
      const double *const __restrict x1ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
      const double *const __restrict y1ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
      const double *const __restrict z1ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
      auto **const __restrict ptr2ptr = soa2.template begin<Particle_T::AttributeNames::ptr>();
      const double *const __restrict x2ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
      const double *const __restrict y2ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
      const double *const __restrict z2ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
      // newton3 ignored: SoAFunctorPair(soa2,soa1) is also called for newton3=false.
      const size_t n1 = soa1.size(), n2 = soa2.size();
      for (size_t i = 0; i < n1; ++i) {
        const size_t iIdx = _particleToIndex.at(ptr1ptr[i]);
        size_t localCount = 0;
        for (size_t j = 0; j < n2; ++j) {
          const double dx = x1ptr[i] - x2ptr[j], dy = y1ptr[i] - y2ptr[j], dz = z1ptr[i] - z2ptr[j];
          if (dx * dx + dy * dy + dz * dz < _interactionLengthSquared) ++localCount;
        }
        _counts[iIdx].value.fetch_add(localCount, std::memory_order_relaxed);
      }
    }

    /**
     * @copydoc autopas::PairwiseFunctor::getNeededAttr()
     */
    constexpr static std::array<typename Particle_T::AttributeNames, 4> getNeededAttr() {
      return {Particle_T::AttributeNames::ptr, Particle_T::AttributeNames::posX, Particle_T::AttributeNames::posY,
              Particle_T::AttributeNames::posZ};
    }

    /**
     * @copydoc autopas::PairwiseFunctor::getComputedAttr()
     */
    constexpr static std::array<typename Particle_T::AttributeNames, 0> getComputedAttr() { return {}; }

   private:
    std::vector<PaddedAtomic> &_counts;
    const std::unordered_map<const Particle_T *, size_t> &_particleToIndex;
    double _interactionLengthSquared;
  };

  /**
   * Pass-2 functor for the two-pass CRS build.
   *
   * The CRS offsets are already set. For each interacting pair (i, j) this functor writes j's index directly into the
   * pre-allocated _neighborList.indices slice for particle i, advancing a per-particle fill cursor with a relaxed
   * atomic fetch-add.
   *
   * The fill cursors reuse the same cache-line-padded PaddedAtomic array as the counter pass (reset to the CRS start
   * offset before this pass runs).
   */
  class VerletListFillerFunctor : public PairwiseFunctor<Particle_T, VerletListFillerFunctor> {
   public:
    /**
     * Structure of the SoAs defined by the particle.
     */
    using SoAArraysType = typename Particle_T::SoAArraysType;
    /**
     * Atomic counter
     */
    using PaddedAtomic = typename VerletListCounterFunctor::PaddedAtomic;

    /**
     * @param neighborList    The CRS structure with offsets already filled and
     *                        indices pre-allocated.
     * @param fillPos         Per-particle fill cursors, initialized to
     *                        neighborList.offsets[i] before this pass starts.
     * @param particleToIndex Map from particle pointer to its dense SoA index.
     * @param interactionLength cutoff + skin.
     */
    VerletListFillerFunctor(NeighborListCRS &neighborList, std::vector<PaddedAtomic> &fillPos,
                            const std::unordered_map<const Particle_T *, size_t> &particleToIndex,
                            double interactionLength)
        : PairwiseFunctor<Particle_T, VerletListFillerFunctor>(interactionLength),
          _neighborList(neighborList),
          _fillPos(fillPos),
          _particleToIndex(particleToIndex),
          _interactionLengthSquared(interactionLength * interactionLength) {}

    std::string getName() override { return "VerletListFillerFunctor"; }
    bool isRelevantForTuning() override { return false; }
    bool allowsNewton3() override { return true; }
    bool allowsNonNewton3() override { return true; }

    /**
     * @copydoc autopas::PairwiseFunctor::AoSFunctor()
     */
    void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) override {
      using namespace autopas::utils::ArrayMath::literals;
      if (i.isDummy() or j.isDummy()) return;
      auto dist = i.getR() - j.getR();
      if (utils::ArrayMath::dot(dist, dist) < _interactionLengthSquared) {
        const size_t iIdx = _particleToIndex.at(&i);
        const size_t jIdx = _particleToIndex.at(&j);
        _neighborList.indices[_fillPos[iIdx].value.fetch_add(1, std::memory_order_relaxed)] = jIdx;
        // newton3 ignored: AoSFunctor(j, i) is also called for newton3=false.
      }
    }

    /**
     * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
     */
    void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
      if (soa.size() == 0) return;
      auto **const __restrict ptrptr = soa.template begin<Particle_T::AttributeNames::ptr>();
      const double *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
      const double *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
      const double *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
      const size_t n = soa.size();
      for (size_t i = 0; i < n; ++i) {
        const size_t iIdx = _particleToIndex.at(ptrptr[i]);
        for (size_t j = i + 1; j < n; ++j) {
          const double dx = xptr[i] - xptr[j], dy = yptr[i] - yptr[j], dz = zptr[i] - zptr[j];
          if (dx * dx + dy * dy + dz * dz < _interactionLengthSquared) {
            const size_t jIdx = _particleToIndex.at(ptrptr[j]);
            _neighborList.indices[_fillPos[iIdx].value.fetch_add(1, std::memory_order_relaxed)] = jIdx;
            if (not newton3) {
              _neighborList.indices[_fillPos[jIdx].value.fetch_add(1, std::memory_order_relaxed)] = iIdx;
            }
          }
        }
      }
    }

    /**
     * @param soa1 SoA of first cell
     * @param soa2 SoA of second cell
     */
    void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool /*newton3*/) override {
      if (soa1.size() == 0 || soa2.size() == 0) return;
      auto **const __restrict ptr1ptr = soa1.template begin<Particle_T::AttributeNames::ptr>();
      const double *const __restrict x1ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
      const double *const __restrict y1ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
      const double *const __restrict z1ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
      auto **const __restrict ptr2ptr = soa2.template begin<Particle_T::AttributeNames::ptr>();
      const double *const __restrict x2ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
      const double *const __restrict y2ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
      const double *const __restrict z2ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
      // newton3 ignored: SoAFunctorPair(soa2,soa1) is also called for newton3=false.
      const size_t n1 = soa1.size(), n2 = soa2.size();
      for (size_t i = 0; i < n1; ++i) {
        const size_t iIdx = _particleToIndex.at(ptr1ptr[i]);
        for (size_t j = 0; j < n2; ++j) {
          const double dx = x1ptr[i] - x2ptr[j], dy = y1ptr[i] - y2ptr[j], dz = z1ptr[i] - z2ptr[j];
          if (dx * dx + dy * dy + dz * dz < _interactionLengthSquared) {
            _neighborList.indices[_fillPos[iIdx].value.fetch_add(1, std::memory_order_relaxed)] =
                _particleToIndex.at(ptr2ptr[j]);
          }
        }
      }
    }

    /**
     * @copydoc autopas::PairwiseFunctor::getNeededAttr()
     */
    constexpr static std::array<typename Particle_T::AttributeNames, 4> getNeededAttr() {
      return {Particle_T::AttributeNames::ptr, Particle_T::AttributeNames::posX, Particle_T::AttributeNames::posY,
              Particle_T::AttributeNames::posZ};
    }

    /**
     * @copydoc autopas::PairwiseFunctor::getComputedAttr()
     */
    constexpr static std::array<typename Particle_T::AttributeNames, 0> getComputedAttr() { return {}; }

   private:
    NeighborListCRS &_neighborList;
    std::vector<PaddedAtomic> &_fillPos;
    const std::unordered_map<const Particle_T *, size_t> &_particleToIndex;
    double _interactionLengthSquared;
  };

  /**
   * This functor checks the validity of neighborhood lists.
   * If a pair of particles has a distance of less than the cutoff radius it
   * checks whether the pair is represented in the CRS neighbor list.
   * If the pair is not present in the list the neighborhood lists are invalid
   * and neighborlistsAreValid() will return false.
   */
  class VerletListValidityCheckerFunctor : public PairwiseFunctor<Particle_T, VerletListValidityCheckerFunctor> {
   public:
    /**
     * Structure of the SoAs defined by the particle.
     */
    using SoAArraysType = typename Particle_T::SoAArraysType;

    /**
     * Constructor
     * @param neighborList  The CRS neighbor list to validate.
     * @param particleIndex Map from particle pointer to its dense SoA index.
     * @param cutoff        The cutoff radius (pairs within this are expected to be listed).
     */
    VerletListValidityCheckerFunctor(const NeighborListCRS &neighborList,
                                     const std::unordered_map<const Particle_T *, size_t> &particleIndex, double cutoff)
        : PairwiseFunctor<Particle_T, VerletListValidityCheckerFunctor>(cutoff),
          _neighborList(neighborList),
          _particleIndex(particleIndex),
          _cutoffsquared(cutoff * cutoff),
          _valid(true) {}

    std::string getName() override { return "VerletListValidityCheckerFunctor"; }

    bool isRelevantForTuning() override { return false; }

    bool allowsNewton3() override {
      utils::ExceptionHandler::exception(
          "VLCAllCellsGeneratorFunctor::allowsNewton3() is not implemented, because it should not be called.");
      return true;
    }

    bool allowsNonNewton3() override {
      utils::ExceptionHandler::exception(
          "VLCAllCellsGeneratorFunctor::allowsNonNewton3() is not implemented, because it should not be called.");
      return true;
    }

    void AoSFunctor(Particle_T &i, Particle_T &j, bool /*newton3*/) override {
      using namespace autopas::utils::ArrayMath::literals;

      auto dist = i.getR() - j.getR();
      double distSquared = utils::ArrayMath::dot(dist, dist);
      if (distSquared < _cutoffsquared) {
        // Thread-safe: reads only from the immutable CRS structure and stack variables.
        const size_t iIdx = _particleIndex.at(&i);
        const size_t jIdx = _particleIndex.at(&j);
        const size_t *beg = _neighborList.begin(iIdx);
        const size_t *end = beg + _neighborList.count(iIdx);
        if (std::find(beg, end, jIdx) == end) {
          // this is thread safe, as _valid is atomic
          _valid = false;
        }
      }
    }

    /**
     * Returns whether the neighbour list are valid.
     * Call this after performing the pairwise traversal
     * @return
     */
    bool neighborlistsAreValid() { return _valid; }

   private:
    const NeighborListCRS &_neighborList;
    const std::unordered_map<const Particle_T *, size_t> &_particleIndex;
    double _cutoffsquared;

    // needs to be thread safe
    std::atomic<bool> _valid;
  };

};  // class VerletListHelpers
}  // namespace autopas

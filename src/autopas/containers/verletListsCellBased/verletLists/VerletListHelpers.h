/**
 * @file VerletListHelpers.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <atomic>

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
   *   - neighbor count : offsets[i+1] - offsets[i]
   *   - neighbor slice : indices.data() + offsets[i]
   */
  struct NeighborListCRS {
    std::vector<size_t> offsets;                            ///< size N+1
    std::vector<size_t, AlignedAllocator<size_t>> indices;  ///< flat neighbor indices

    /// Number of particles tracked by this list.
    [[nodiscard]] size_t size() const { return offsets.empty() ? 0u : offsets.size() - 1u; }
    /// Number of neighbors of particle i.
    [[nodiscard]] size_t count(size_t i) const { return offsets[i + 1] - offsets[i]; }
    /// Pointer to the first neighbor index of particle i (const).
    [[nodiscard]] const size_t *begin(size_t i) const { return indices.data() + offsets[i]; }
    /// Pointer to the first neighbor index of particle i (mutable).
    [[nodiscard]] size_t *begin(size_t i) { return indices.data() + offsets[i]; }
  };

  /**
   * This functor can generate verlet lists using the typical pairwise traversal.
   */
  class VerletListGeneratorFunctor : public PairwiseFunctor<Particle_T, VerletListGeneratorFunctor> {
   public:
    /**
     * Structure of the SoAs defined by the particle.
     */
    using SoAArraysType = typename Particle_T::SoAArraysType;

    /**
     * Constructor
     * @param neighborLists Temporary ragged list: one inner vector per particle index.
     *                      Populated during the linked-cells traversal; used by
     *                      VerletLists::buildCRS() to assemble the final flat CRS.
     * @param particleToIndex Map from particle pointers to their dense SoA index.
     *                      Built once by VerletLists::buildParticleIndex().
     * @param interactionLength cutoff + skin radius.
     */
    VerletListGeneratorFunctor(std::vector<std::vector<size_t>> &neighborLists,
                               const std::unordered_map<const Particle_T *, size_t> &particleToIndex,
                               double interactionLength)
        : PairwiseFunctor<Particle_T, VerletListGeneratorFunctor>(interactionLength),
          _neighborLists(neighborLists),
          _particleToIndex(particleToIndex),
          _interactionLengthSquared(interactionLength * interactionLength) {}

    std::string getName() override { return "VerletListGeneratorFunctor"; }

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

      if (i.isDummy() or j.isDummy()) {
        return;
      }
      auto dist = i.getR() - j.getR();

      double distSquared = utils::ArrayMath::dot(dist, dist);
      if (distSquared < _interactionLengthSquared) {
        // Thread-safe: particle i lives in a specific cell, and each cell is
        // accessed by exactly one thread at a time (guaranteed by the traversal).
        // _neighborLists[iIdx] is therefore never accessed concurrently.
        _neighborLists[_particleToIndex.at(&i)].push_back(_particleToIndex.at(&j));
        // No newton3 here: AoSFunctor(j, i) is also called when newton3=false.
      }
    }

    /**
     * SoAFunctor for verlet list generation. (single cell version)
     * @param soa the soa
     * @param newton3 whether to use newton 3
     */
    void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
      if (soa.size() == 0) return;

      auto **const __restrict ptrptr = soa.template begin<Particle_T::AttributeNames::ptr>();
      const double *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
      const double *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
      const double *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();

      size_t numPart = soa.size();
      for (unsigned int i = 0; i < numPart; ++i) {
        const size_t iIdx = _particleToIndex.at(ptrptr[i]);
        auto &currentList = _neighborLists[iIdx];

        for (unsigned int j = i + 1; j < numPart; ++j) {
          const double drx = xptr[i] - xptr[j];
          const double dry = yptr[i] - yptr[j];
          const double drz = zptr[i] - zptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double distSquared = drx2 + dry2 + drz2;

          if (distSquared < _interactionLengthSquared) {
            const size_t jIdx = _particleToIndex.at(ptrptr[j]);
            currentList.push_back(jIdx);
            if (not newton3) {
              // SoAFunctorSingle is called once for both newton3=true and newton3=false,
              // so we must explicitly record the reverse direction here.
              _neighborLists[jIdx].push_back(iIdx);
            }
          }
        }
      }
    }

    /**
     * SoAFunctor for the verlet list generation. (two cell version)
     * @param soa1 soa of first cell
     * @param soa2 soa of second cell
     * @note newton3 is ignored here, as for newton3=false SoAFunctorPair(soa2, soa1) will also be called.
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

      // newton3 is ignored: for newton3=false, SoAFunctorPair(soa2, soa1) is also called by the traversal.
      size_t numPart1 = soa1.size();
      for (unsigned int i = 0; i < numPart1; ++i) {
        const size_t iIdx = _particleToIndex.at(ptr1ptr[i]);
        auto &currentList = _neighborLists[iIdx];

        size_t numPart2 = soa2.size();

        for (unsigned int j = 0; j < numPart2; ++j) {
          const double drx = x1ptr[i] - x2ptr[j];
          const double dry = y1ptr[i] - y2ptr[j];
          const double drz = z1ptr[i] - z2ptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double distSquared = drx2 + dry2 + drz2;

          if (distSquared < _interactionLengthSquared) {
            currentList.push_back(_particleToIndex.at(ptr2ptr[j]));
          }
        }
      }
    }

    /**
     * @copydoc autopas::Functor::getNeededAttr()
     */
    constexpr static std::array<typename Particle_T::AttributeNames, 4> getNeededAttr() {
      return std::array<typename Particle_T::AttributeNames, 4>{
          Particle_T::AttributeNames::ptr, Particle_T::AttributeNames::posX, Particle_T::AttributeNames::posY,
          Particle_T::AttributeNames::posZ};
    }

    /**
     * @copydoc autopas::Functor::getComputedAttr()
     */
    constexpr static std::array<typename Particle_T::AttributeNames, 0> getComputedAttr() {
      return std::array<typename Particle_T::AttributeNames, 0>{/*Nothing*/};
    }

   private:
    std::vector<std::vector<size_t>> &_neighborLists;
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

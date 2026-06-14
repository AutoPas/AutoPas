/**
 * @file VerletListHelpers.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <atomic>

#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
namespace autopas {
class PointerToIndexMap {
 public:
  template <typename PointerType>
  void build(const std::vector<PointerType> &keys) {
    size_t size = keys.size();
    size_t tableSize = 1;
    while (tableSize < size * 2) tableSize *= 2;
    _table.assign(tableSize, Entry{nullptr, 0});
    _mask = tableSize - 1;

    for (size_t val = 0; val < size; ++val) {
      const void *key = static_cast<const void *>(keys[val]);
      size_t h = (reinterpret_cast<size_t>(key) >> 3) & _mask;
      while (_table[h].key != nullptr) {
        h = (h + 1) & _mask;
      }
      _table[h] = {key, val};
    }
  }

  [[nodiscard]] inline size_t at(const void *key) const noexcept {
    size_t h = (reinterpret_cast<size_t>(key) >> 3) & _mask;
    while (_table[h].key != key) {
      h = (h + 1) & _mask;
    }
    return _table[h].val;
  }

 private:
  struct Entry {
    const void *key;
    size_t val;
  };
  std::vector<Entry> _table;
  size_t _mask = 0;
};
/**
 * Class of helpers for the VerletLists class.
 * @tparam Particle_T
 */
template <class Particle_T>
class VerletListHelpers {
 public:
  /**
   * Neighbor list AoS style.
   */
  using NeighborListAoSType = std::unordered_map<Particle_T *, std::vector<Particle_T *>>;

  /**
   * Neighbor pairs list AoS style.
   */
  using NeighborPairsListAoSType = std::unordered_map<Particle_T *, std::vector<std::pair<Particle_T *, Particle_T *>>>;

  class CRSNeighborList {
   public:
    void clear() {
      _offsets.clear();
      _neighbors.clear();
    }

    void resizeParticles(size_t numParticles) {
      _offsets.assign(numParticles + 1, 0);
      _neighbors.clear();
    }

    [[nodiscard]] size_t size() const noexcept { return _offsets.empty() ? 0 : _offsets.size() - 1; }

    [[nodiscard]] std::span<const size_t> neighborsOf(size_t particleIndex) const noexcept {
      return {_neighbors.data() + _offsets[particleIndex], _offsets[particleIndex + 1] - _offsets[particleIndex]};
    }

    [[nodiscard]] std::span<size_t> neighborsOf(size_t particleIndex) noexcept {
      return {_neighbors.data() + _offsets[particleIndex], _offsets[particleIndex + 1] - _offsets[particleIndex]};
    }

    std::vector<size_t> &offsets() noexcept { return _offsets; }
    std::vector<size_t> &neighbors() noexcept { return _neighbors; }

    [[nodiscard]] const std::vector<size_t> &offsets() const noexcept { return _offsets; }
    [[nodiscard]] const std::vector<size_t> &neighbors() const noexcept { return _neighbors; }

   private:
    std::vector<size_t> _offsets;
    std::vector<size_t> _neighbors;
  };

  class CRSPairNeighborList {
   public:
    void clear() {
      _offsets.clear();
      _neighborPairs.clear();
    }

    void resizeParticles(size_t numParticles) {
      _offsets.assign(numParticles + 1, 0);
      _neighborPairs.clear();
    }

    [[nodiscard]] size_t size() const noexcept { return _offsets.empty() ? 0 : _offsets.size() - 1; }

    [[nodiscard]] std::span<const std::pair<Particle_T *, Particle_T *>> neighborPairsOf(
        size_t particleIndex) const noexcept {
      return {_neighborPairs.data() + _offsets[particleIndex], _offsets[particleIndex + 1] - _offsets[particleIndex]};
    }

    [[nodiscard]] std::span<std::pair<Particle_T *, Particle_T *>> neighborPairsOf(size_t particleIndex) noexcept {
      return {_neighborPairs.data() + _offsets[particleIndex], _offsets[particleIndex + 1] - _offsets[particleIndex]};
    }

    std::vector<size_t> &offsets() noexcept { return _offsets; }
    std::vector<std::pair<Particle_T *, Particle_T *>> &neighborPairs() noexcept { return _neighborPairs; }

    [[nodiscard]] const std::vector<size_t> &offsets() const noexcept { return _offsets; }
    [[nodiscard]] const std::vector<std::pair<Particle_T *, Particle_T *>> &neighborPairs() const noexcept {
      return _neighborPairs;
    }

   private:
    std::vector<size_t> _offsets;
    std::vector<std::pair<Particle_T *, Particle_T *>> _neighborPairs;
  };

  class CRSNeighborCounterFunctor : public PairwiseFunctor<Particle_T, CRSNeighborCounterFunctor> {
   public:
    using SoAArraysType = typename Particle_T::SoAArraysType;

    CRSNeighborCounterFunctor(const PointerToIndexMap &particleToIndex, std::vector<std::size_t> &counts,
                              double interactionLength)
        : PairwiseFunctor<Particle_T, CRSNeighborCounterFunctor>(interactionLength),
          _particleToIndex(particleToIndex),
          _counts(counts),
          _interactionLengthSquared(interactionLength * interactionLength) {}

    std::string getName() override { return "CRSNeighborCounterFunctor"; }

    bool isRelevantForTuning() override { return false; }

    bool allowsNewton3() override {
      utils::ExceptionHandler::exception(
          "CRSNeighborCounterFunctor::allowsNewton3() is not implemented, because it should not be called.");
      return true;
    }

    bool allowsNonNewton3() override {
      utils::ExceptionHandler::exception(
          "CRSNeighborCounterFunctor::allowsNonNewton3() is not implemented, because it should not be called.");
      return true;
    }

    void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
      const auto dr = utils::ArrayMath::sub(i.getR(), j.getR());
      const auto dist2 = utils::ArrayMath::dot(dr, dr);

      if (dist2 <= _interactionLengthSquared) {
        const auto indexI = _particleToIndex.at(&i);
        ++_counts[indexI + 1];
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
    const PointerToIndexMap &_particleToIndex;
    std::vector<std::size_t> &_counts;
    double _interactionLengthSquared;
  };

  class CRSNeighborFillFunctor : public PairwiseFunctor<Particle_T, CRSNeighborFillFunctor> {
   public:
    using SoAArraysType = typename Particle_T::SoAArraysType;

    CRSNeighborFillFunctor(const PointerToIndexMap &particleToIndex, std::vector<std::size_t> &writeOffsets,
                           std::vector<std::size_t> &neighbors, double interactionLength)
        : PairwiseFunctor<Particle_T, CRSNeighborFillFunctor>(interactionLength),
          _particleToIndex(particleToIndex),
          _writeOffsets(writeOffsets),
          _neighbors(neighbors),
          _interactionLengthSquared(interactionLength * interactionLength) {}

    std::string getName() override { return "CRSNeighborFillFunctor"; }

    bool isRelevantForTuning() override { return false; }

    bool allowsNewton3() override {
      utils::ExceptionHandler::exception(
          "CRSNeighborFillFunctor::allowsNewton3() is not implemented, because it should not be called.");
      return true;
    }

    bool allowsNonNewton3() override {
      utils::ExceptionHandler::exception(
          "CRSNeighborFillFunctor::allowsNonNewton3() is not implemented, because it should not be called.");
      return true;
    }

    void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
      const auto dr = utils::ArrayMath::sub(i.getR(), j.getR());
      const auto dist2 = utils::ArrayMath::dot(dr, dr);

      if (dist2 <= _interactionLengthSquared) {
        const auto indexI = _particleToIndex.at(&i);
        const auto indexJ = _particleToIndex.at(&j);

        _neighbors[_writeOffsets[indexI]++] = indexJ;
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
    const PointerToIndexMap &_particleToIndex;
    std::vector<std::size_t> &_writeOffsets;
    std::vector<std::size_t> &_neighbors;
    double _interactionLengthSquared;
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
     * @param verletListsAoS
     * @param interactionLength
     */
    VerletListGeneratorFunctor(NeighborListAoSType &verletListsAoS, double interactionLength)
        : PairwiseFunctor<Particle_T, VerletListGeneratorFunctor>(interactionLength),
          _verletListsAoS(verletListsAoS),
          _interactionLengthSquared(interactionLength * interactionLength) {}

    std::string getName() override { return "VerletListGeneratorFunctor"; }

    bool isRelevantForTuning() override { return false; }

    bool allowsNewton3() override {
      utils::ExceptionHandler::exception(
          "VerletListGeneratorFunctor::allowsNewton3() is not implemented, because it should not be called.");
      return true;
    }

    bool allowsNonNewton3() override {
      utils::ExceptionHandler::exception(
          "VerletListGeneratorFunctor::allowsNonNewton3() is not implemented, because it should not be called.");
      return true;
    }

    void AoSFunctor(Particle_T &i, Particle_T &j, bool /*newton3*/) override {
      using namespace autopas::utils::ArrayMath::literals;

      if (i.isDummy() or j.isDummy()) {
        return;
      }
      auto dist = i.getR() - j.getR();

      double distsquare = utils::ArrayMath::dot(dist, dist);
      if (distsquare < _interactionLengthSquared) {
        // this is thread safe, only if particle i is accessed by only one
        // thread at a time. which is ensured, as particle i resides in a
        // specific cell and each cell is only accessed by one thread at a time
        // (ensured by traversals)
        // also the list is not allowed to be resized!

        _verletListsAoS.at(&i).push_back(&j);
        // no newton3 here, as AoSFunctor(j,i) will also be called if newton3 is disabled.
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
        auto &currentList = _verletListsAoS.at(ptrptr[i]);

        for (unsigned int j = i + 1; j < numPart; ++j) {
          const double drx = xptr[i] - xptr[j];
          const double dry = yptr[i] - yptr[j];
          const double drz = zptr[i] - zptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double dr2 = drx2 + dry2 + drz2;

          if (dr2 < _interactionLengthSquared) {
            currentList.push_back(ptrptr[j]);
            if (not newton3) {
              // we need this here, as SoAFunctorSingle will only be called once for both newton3=true and false.
              _verletListsAoS.at(ptrptr[j]).push_back(ptrptr[i]);
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

      size_t numPart1 = soa1.size();
      for (unsigned int i = 0; i < numPart1; ++i) {
        auto &currentList = _verletListsAoS.at(ptr1ptr[i]);

        size_t numPart2 = soa2.size();

        for (unsigned int j = 0; j < numPart2; ++j) {
          const double drx = x1ptr[i] - x2ptr[j];
          const double dry = y1ptr[i] - y2ptr[j];
          const double drz = z1ptr[i] - z2ptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double dr2 = drx2 + dry2 + drz2;

          if (dr2 < _interactionLengthSquared) {
            currentList.push_back(ptr2ptr[j]);
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
    NeighborListAoSType &_verletListsAoS;
    double _interactionLengthSquared;
  };

  /**
   * This functor can generate verlet lists of neighbor pairs for a triwise traversal.
   */
  class PairVerletListGeneratorFunctor : public TriwiseFunctor<Particle_T, PairVerletListGeneratorFunctor> {
   public:
    /**
     * Structure of the SoAs defined by the particle.
     */
    using SoAArraysType = typename Particle_T::SoAArraysType;

    /**
     * Constructor
     * @param pairVerletListsAoS
     * @param interactionLength
     */
    PairVerletListGeneratorFunctor(std::vector<std::vector<std::pair<Particle_T *, Particle_T *>>> &pairVerletListsAoS,
                                   const PointerToIndexMap &particlePtr2indexMap, double interactionLength)
        : TriwiseFunctor<Particle_T, PairVerletListGeneratorFunctor>(interactionLength),
          _pairVerletListsAoS(pairVerletListsAoS),
          _particlePtr2indexMap(particlePtr2indexMap),
          _interactionLengthSquared(interactionLength * interactionLength) {}

    std::string getName() override { return "PairVerletListGeneratorFunctor"; }

    bool isRelevantForTuning() override { return false; }

    bool allowsNewton3() override {
      utils::ExceptionHandler::exception(
          "PairVerletListGeneratorFunctor::allowsNewton3() is not implemented, because it should not be called.");
      return true;
    }

    bool allowsNonNewton3() override {
      utils::ExceptionHandler::exception(
          "PairVerletListGeneratorFunctor::allowsNonNewton3() is not implemented, because it should not be called.");
      return true;
    }

    void AoSFunctor(Particle_T &i, Particle_T &j, Particle_T &k, bool /*newton3*/) override {
      using namespace autopas::utils::ArrayMath::literals;
      if (i.isDummy() or j.isDummy() or k.isDummy()) return;

      const auto distIJ = j.getR() - i.getR();
      const auto distIK = k.getR() - i.getR();
      const auto distJK = k.getR() - j.getR();

      if (utils::ArrayMath::dot(distIJ, distIJ) < _interactionLengthSquared and
          utils::ArrayMath::dot(distIK, distIK) < _interactionLengthSquared and
          utils::ArrayMath::dot(distJK, distJK) < _interactionLengthSquared) {
        const auto indexI = _particlePtr2indexMap.at(&i);
        _pairVerletListsAoS[indexI].push_back(std::make_pair(&j, &k));
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
    std::vector<std::vector<std::pair<Particle_T *, Particle_T *>>> &_pairVerletListsAoS;
    const PointerToIndexMap &_particlePtr2indexMap;
    double _interactionLengthSquared;
  };

  /**
   * This functor checks the validity of neighborhood lists.
   * If a pair of particles has a distance of less than the cutoff radius it
   * checks whether the pair is represented in the verlet list.
   * If the pair is not present in the list the neigborhood lists are invalid
   * and neighborlistsAreValid()  will return false.
   * @todo: SoA?
   */
  class VerletListValidityCheckerFunctor : public PairwiseFunctor<Particle_T, VerletListValidityCheckerFunctor> {
   public:
    /**
     * Structure of the SoAs defined by the particle.
     */
    using SoAArraysType = typename Particle_T::SoAArraysType;

    /**
     * Constructor
     * @param verletListsAoS
     * @param cutoff
     */
    VerletListValidityCheckerFunctor(NeighborListAoSType &verletListsAoS, double cutoff)
        : PairwiseFunctor<Particle_T, VerletListValidityCheckerFunctor>(cutoff),
          _verletListsAoS(verletListsAoS),
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

    void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) override {
      using namespace autopas::utils::ArrayMath::literals;

      auto dist = i.getR() - j.getR();
      double distsquare = utils::ArrayMath::dot(dist, dist);
      if (distsquare < _cutoffsquared) {
        // this is thread safe, we have variables on the stack
        auto found = std::find(_verletListsAoS[&i].begin(), _verletListsAoS[&i].end(), &j);
        if (found == _verletListsAoS[&i].end()) {
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
    NeighborListAoSType &_verletListsAoS;
    double _cutoffsquared;

    // needs to be thread safe
    std::atomic<bool> _valid;
  };

};  // class VerletListHelpers
}  // namespace autopas

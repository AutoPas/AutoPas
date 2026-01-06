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
   * Neighbor list AoS style.
   */
  using NeighborListAoSType = std::unordered_map<Particle_T *, std::vector<Particle_T *>>;

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

  class VerletListGeneratorFunctorSoA : public PairwiseFunctor<Particle_T, VerletListGeneratorFunctor> {
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
    VerletListGeneratorFunctorSoA(NeighborListSoAType &verletListsSoA, std::unordered_map<const Particle_T *,
                                  size_t> particlePtr2indexMap, double interactionLength)
    : PairwiseFunctor<Particle_T, VerletListGeneratorFunctor>(interactionLength),
      _verletListsSoA(verletListsSoA),
      _particlePtr2indexMap(particlePtr2indexMap),
      _interactionLengthSquared(interactionLength * interactionLength){}

    std::string getName() override { return "VerletListGeneratorFunctorSoA"; }

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

    /**
     * SoAFunctor for verlet list generation. (single cell version)
     * @param soa the soa
     * @param newton3 whether to use newton 3
     */
    //ex soa is the SoABuffer of the current cell
    void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
      if (soa.size() == 0) return;

      auto **const __restrict ptrptr = soa.template begin<Particle_T::AttributeNames::ptr>();
      const double *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
      const double *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
      const double *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();

      size_t numPart = soa.size();
      for (unsigned int i = 0; i < numPart; ++i) {
        auto &currentListSoAI = _verletListsSoA[_particlePtr2indexMap.at(ptrptr[i])];

        for (unsigned int j = i + 1; j < numPart; ++j) {
          const double drx = xptr[i] - xptr[j];
          const double dry = yptr[i] - yptr[j];
          const double drz = zptr[i] - zptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double dr2 = drx2 + dry2 + drz2;

          if (dr2 < this->_interactionLengthSquared) {
            currentListSoAI.push_back(_particlePtr2indexMap.at(ptrptr[j]));
            if (not newton3) {
              // we need this here, as SoAFunctorSingle will only be called once for both newton3=true and false.
              _verletListsSoA[_particlePtr2indexMap.at(ptrptr[j])].push_back(_particlePtr2indexMap.at(ptrptr[i]));
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

        size_t numPart2 = soa2.size();
        auto  &currentListSoA = _verletListsSoA[_particlePtr2indexMap.at(ptr1ptr[i])];

        for (unsigned int j = 0; j < numPart2; ++j) {
          const double drx = x1ptr[i] - x2ptr[j];
          const double dry = y1ptr[i] - y2ptr[j];
          const double drz = z1ptr[i] - z2ptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double dr2 = drx2 + dry2 + drz2;

          if (dr2 < this->_interactionLengthSquared) {
          currentListSoA.push_back(_particlePtr2indexMap.at(ptr2ptr[j]));
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
    NeighborListSoAType &_verletListsSoA;
    std::unordered_map<const Particle_T *, size_t> _particlePtr2indexMap;
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

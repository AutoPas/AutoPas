/**
 * @file VerletListHelpers.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <atomic>

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
namespace autopas {

/**
 * Class of helpers for the VerletLists class.
 * @tparam Particle
 */
template <class Particle>
class VerletListHelpers {
 public:
  /**
   * Neighbor list AoS style.
   */
  using NeighborListAoSType = std::unordered_map<Particle *, std::vector<Particle *>>;

  /**
   * This functor can generate verlet lists using the typical pairwise traversal.
   */
  class VerletListGeneratorFunctor : public Functor<Particle, VerletListGeneratorFunctor> {
   public:
    /**
     * Structure of the SoAs defined by the particle.
     */
    using SoAArraysType = typename Particle::SoAArraysType;

    /**
     * Constructor
     * @param verletListsAoS
     * @param interactionLength
     */
    VerletListGeneratorFunctor(NeighborListAoSType &verletListsAoS, double interactionLength)
        : Functor<Particle, VerletListGeneratorFunctor>(interactionLength),
          _verletListsAoS(verletListsAoS),
          _interactionLengthSquared(interactionLength * interactionLength) {}

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

    void AoSFunctor(Particle &i, Particle &j, bool /*newton3*/) override {
      if (i.isDummy() or j.isDummy()) {
        return;
      }
      auto dist = utils::ArrayMath::sub(i.getR(), j.getR());

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
      if (soa.getNumberOfParticles() == 0) return;

      auto **const __restrict ptrptr = soa.template begin<Particle::AttributeNames::ptr>();
      double *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
      double *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
      double *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

      size_t numPart = soa.getNumberOfParticles();
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
      if (soa1.getNumberOfParticles() == 0 || soa2.getNumberOfParticles() == 0) return;

      auto **const __restrict ptr1ptr = soa1.template begin<Particle::AttributeNames::ptr>();
      double *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
      double *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
      double *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();

      auto **const __restrict ptr2ptr = soa2.template begin<Particle::AttributeNames::ptr>();
      double *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
      double *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
      double *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

      size_t numPart1 = soa1.getNumberOfParticles();
      for (unsigned int i = 0; i < numPart1; ++i) {
        auto &currentList = _verletListsAoS.at(ptr1ptr[i]);

        size_t numPart2 = soa2.getNumberOfParticles();

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
     * @copydoc Functor::getNeededAttr()
     */
    constexpr static std::array<typename Particle::AttributeNames, 4> getNeededAttr() {
      return std::array<typename Particle::AttributeNames, 4>{
          Particle::AttributeNames::ptr, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
          Particle::AttributeNames::posZ};
    }

    /**
     * @copydoc Functor::getComputedAttr()
     */
    constexpr static std::array<typename Particle::AttributeNames, 0> getComputedAttr() {
      return std::array<typename Particle::AttributeNames, 0>{/*Nothing*/};
    }

   private:
    NeighborListAoSType &_verletListsAoS;
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
  class VerletListValidityCheckerFunctor : public Functor<Particle, VerletListValidityCheckerFunctor> {
   public:
    /**
     * Structure of the SoAs defined by the particle.
     */
    using SoAArraysType = typename Particle::SoAArraysType;

    /**
     * Constructor
     * @param verletListsAoS
     * @param cutoff
     */
    VerletListValidityCheckerFunctor(NeighborListAoSType &verletListsAoS, double cutoff)
        : Functor<Particle, VerletListValidityCheckerFunctor>(cutoff),
          _verletListsAoS(verletListsAoS),
          _cutoffsquared(cutoff * cutoff),
          _valid(true) {}

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

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      auto dist = utils::ArrayMath::sub(i.getR(), j.getR());
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

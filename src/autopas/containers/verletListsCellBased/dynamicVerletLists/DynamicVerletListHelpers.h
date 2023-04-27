/**
 * @file DynamicVerletListHelpers
 * @author Luis Gall
 * @date 26.04.23
*/

#pragma once

#include <atomic>

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"

namespace autopas {

/**
 * Class of helpers for the DynamicVerletLists class.
 * @tparam Particle
 */
template <class Particle>
class DynamicVerletListHelpers {
 public:
  /**
   * Neighbor list AoS style
   */
   using AoSNeighborListType = std::unordered_map<Particle *,  std::vector<Particle*>>;

   class DynamicVerletListGeneratorFunctor : public Functor<Particle, DynamicVerletListGeneratorFunctor> {
    public:

     DynamicVerletListGeneratorFunctor(AoSNeighborListType &verletListsAoS, std::unordered_map<Particle*, std::array<double, 3>> &particlePtr2rebuildPosition, double interactionLength)
      : Functor<Particle, DynamicVerletListGeneratorFunctor>(interactionLength),
          _verletListsAoS(verletListsAoS),
           _particlePtr2rebuildPosition(particlePtr2rebuildPosition),
           _interactionLengthSquared(interactionLength * interactionLength) {}

    bool isRelevantForTuning() override { return false; }

    bool allowsNewton3() override {
      utils::ExceptionHandler::exception("DynamicVerletListGeneratorFunctor::allowsNewton3() is not implemented because is should not be called");
      return true;
    }

    bool allowsNonNewton3() override {
      utils::ExceptionHandler::exception("DynamicVerletListGeneratorFunctor::allowsNonNewton3() is not implemented because is should not be called");
      return true;
    }

    void AoSFunctor(Particle &i, Particle &j, bool /*newtone3*/) override {
      if (i.isDummy() or j.isDummy()) {
        return;
      }

      auto distance = utils::ArrayMath::sub(i.getR(), j.getR());

      double distanceSquare = utils::ArrayMath::dot(distance, distance);
      if (distanceSquare < _interactionLengthSquared) {
        _verletListsAoS.at(&i).push_back(&j);
        _particlePtr2rebuildPosition.at(&i) = i.getR();
      }
    }

    /**
     * @copydoc autopas::Functor::getNeededAttr()
     */
    constexpr static std::array<typename Particle::AttributeNames, 4> getNeededAttr() {
      return std::array<typename Particle::AttributeNames, 4>{
          Particle::AttributeNames::ptr, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
          Particle::AttributeNames::posZ};
    }

    /**
     * @copydoc autopas::Functor::getComputedAttr()
     */
    constexpr static std::array<typename Particle::AttributeNames, 0> getComputedAttr() {
      return std::array<typename Particle::AttributeNames, 0>{/*Nothing*/};
    }

    private:
     AoSNeighborListType & _verletListsAoS;
     double _interactionLengthSquared;
     std::unordered_map<Particle*, std::array<double, 3>> _particlePtr2rebuildPosition;
   };
};

}
/**
 * @file DynamicVerletLists.h
 * @author Luis Gall
 * @date 26.04.23
 */

#pragma once

#include "DynamicVerletListHelpers.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

/**
 * Dynamic Verlet Lists container
 * This container is an improvement for the general VerletLists container
 * @tparam Particle
 */
template <class Particle>
class DynamicVerletLists : public VerletLists<Particle> {
  using LinkedParticleCell = FullParticleCell<Particle>;

 public:

  DynamicVerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
              const double skinPerTimestep, const unsigned int rebuildFrequency,
              const typename VerletLists<Particle>::BuildVerletListType buildVerletListType = VerletLists<Particle>::BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletLists<Particle>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency,
                                        buildVerletListType, cellSizeFactor) {}
  
  bool neighborListsAreValid() {

    auto halfSkinSquare = (this->getVerletSkin() * this->getVerletSkin()) / 4;

    for (auto& [particlePtr, neighborLists] : this->getVerletListsAoS()) {

      auto distance = utils::ArrayMath::sub(particlePtr->getR(), _particlePtr2rebuildPositionMap.at(particlePtr));
      double distanceSquare = utils::ArrayMath::dot(distance, distance);

      if (distanceSquare >=  halfSkinSquare) {
        return false;
      }
    }

    return true;
  }
  
 protected:

  size_t generateRebuildPositionMap() {
    size_t numParticles = 0;
    _particlePtr2rebuildPositionMap.clear();
    
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++numParticles) {
      _particlePtr2rebuildPositionMap[&(*iter)];
    }
    return numParticles;
  }

  // TODO : think about templating this function and reuse it from VerletLists.h
  void updateVerletListsAoS(bool useNewton3) override {
    VerletLists<Particle>::generateAoSNeighborLists();
    generateRebuildPositionMap();
    typename DynamicVerletListHelpers<Particle>::DynamicVerletListGeneratorFunctor f(this->getVerletListsAoS(), _particlePtr2rebuildPositionMap,
                                                                                     this->getCutoff() + this->getVerletSkin());
    /// @todo autotune traversal
    switch (this->_buildVerletListType) {
      case VerletLists<Particle>::BuildVerletListType::VerletAoS: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal =
              LCC08Traversal<LinkedParticleCell, typename DynamicVerletListHelpers<Particle>::DynamicVerletListGeneratorFunctor,
                             DataLayoutOption::aos, theBool>(
                  this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
                  this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      case VerletLists<Particle>::BuildVerletListType::VerletSoA: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal =
              LCC08Traversal<LinkedParticleCell, typename DynamicVerletListHelpers<Particle>::DynamicVerletListGeneratorFunctor,
                             DataLayoutOption::soa, theBool>(
                  this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
                  this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      default:
        utils::ExceptionHandler::exception("DynamicVerletLists::updateVerletListsAoS(): unsupported BuildVerletListType: {}",
                                           this->_buildVerletListType);
        break;
    }
    
    this->_soaListIsValid = false;
  }
  
 private:
  
  std::unordered_map<Particle*, std::array<double, 3>> _particlePtr2rebuildPositionMap;
  
};
}
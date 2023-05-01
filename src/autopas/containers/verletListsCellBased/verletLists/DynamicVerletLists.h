/**
 * @file DynamicVerletLists.h
 * @author Luis Gall
 * @date 26.04.23
 */

#pragma once

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

 public:

  DynamicVerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
              const double skinPerTimestep, const unsigned int rebuildFrequency,
              const typename VerletLists<Particle>::BuildVerletListType buildVerletListType = VerletLists<Particle>::BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletLists<Particle>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency,
                                        buildVerletListType, cellSizeFactor) {}
  
  bool neighborListsAreValid() override {

    auto halfSkinSquare = (this->getVerletSkin() * this->getVerletSkin()) / 4;
    bool listInvalid = false;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(|| : listInvalid)
#endif
    for (auto& [particlePtr, neighborLists] : this->getVerletListsAoS()) {
      auto rebuildPosition = _particlePtr2rebuildPositionMap.at(particlePtr);
      auto distance = utils::ArrayMath::sub(particlePtr->getR(), rebuildPosition);
      double distanceSquare = utils::ArrayMath::dot(distance, distance);

      if (distanceSquare >=  halfSkinSquare) {
        listInvalid = true;
      }
    }

    return !listInvalid;
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    generateRebuildPositionMap();
    VerletLists<Particle>::rebuildNeighborLists(traversal);
  }

 private:

  void generateRebuildPositionMap() {
    _particlePtr2rebuildPositionMap.clear();
    _particlePtr2rebuildPositionMap.reserve(this->getVerletListsAoS().size());

    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      _particlePtr2rebuildPositionMap[&(*iter)] = (*iter).getR();
    }
  }
  
  std::unordered_map<Particle*, std::array<double, 3>> _particlePtr2rebuildPositionMap;
  
};
}
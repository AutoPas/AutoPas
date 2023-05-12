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
    bool partialRebuilding = false;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for reduction(|| : listInvalid, partialRebuilding) schedule(static, 50)
#endif
    for (auto& particlePositionPair : _particlePtr2rebuildPositionBuffer) {
      auto distance = utils::ArrayMath::sub(particlePositionPair.first->getR(), particlePositionPair.second);
      double distanceSquare = utils::ArrayMath::dot(distance, distance);

      if (distanceSquare >=  halfSkinSquare) {
        auto& cell = this->_linkedCells.getCellBlock().getContainingCell(particlePositionPair.first->getR());
        cell.setDirty(true);
        partialRebuilding = true;
        listInvalid = true;
      }
    }

    this->_partialRebuilding = partialRebuilding;
    return !listInvalid;
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    generateRebuildPositionMap();
    VerletLists<Particle>::rebuildNeighborLists(traversal);
  }

 private:

  void generateRebuildPositionMap() {
    _particlePtr2rebuildPositionBuffer.clear();
    _particlePtr2rebuildPositionBuffer.reserve(this->getVerletListsAoS().size());

    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      std::pair<Particle*, std::array<double, 3>> particlePositionPair = std::make_pair(&(*iter), (*iter).getR());
      _particlePtr2rebuildPositionBuffer.emplace_back(particlePositionPair);
    }
  }
  
  std::vector<std::pair<Particle*, std::array<double, 3>>> _particlePtr2rebuildPositionBuffer;
  
};
}
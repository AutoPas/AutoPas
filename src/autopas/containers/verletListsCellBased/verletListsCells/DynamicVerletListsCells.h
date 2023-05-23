/**
* @file DynamicVerletListsCells.h
* @author Luis Gall
* @date 06.05.23
*/

#pragma once

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"

namespace autopas {

template <class Particle, class NeighborList>
class DynamicVerletListsCells : public VerletListsCells<Particle, NeighborList> {

 public:

  DynamicVerletListsCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                          const double skinPerTimestep = 0, const unsigned int rebuildFrequency = 2,
                          const double cellSizeFactor = 1.0,
                          const LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell,
                          typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType =
                              VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild)
  : VerletListsCells<Particle, NeighborList>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency, cellSizeFactor, loadEstimator, buildType) {}

  bool neighborListsAreValid() override {

    auto halfSkinSquare = (this->getVerletSkin() * this->getVerletSkin()) / 4;
    bool listInvalid = false;
    auto partialRebuilding = false;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for reduction(|| : listInvalid, partialRebuilding) schedule(static, 50)
#endif
    for (auto& particlePositionPair : _particlePtr2rebuildPositionBuffer) {
      auto distance = utils::ArrayMath::sub(particlePositionPair.first->getR(), particlePositionPair.second);
      double distanceSquare = utils::ArrayMath::dot(distance, distance);

      if (distanceSquare >=  halfSkinSquare) {
        listInvalid = true;
        auto & cell = this->_linkedCells.getCellBlock().getContainingCell(particlePositionPair.first->getR());
        cell.setDirty(true);

        // TODO : set dirty flag of this neighboring cells as well
        partialRebuilding = true;
        }
    }
    this->_partialRebuilding = partialRebuilding;
    return !listInvalid;
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    VerletListsCells<Particle, NeighborList>::rebuildNeighborLists(traversal);
    generateRebuildPositionMap();
  }

  [[nodiscard]] ContainerOption getContainerType() const override {
    if (this->_neighborList.getContainerType() == ContainerOption::pairwiseVerletLists) {
      return ContainerOption::dynamicPairwiseVerletLists;
    }
    else if(this->_neighborList.getContainerType() == ContainerOption::verletListsCells) {
      return ContainerOption::dynamicVerletListsCells;
    }
    else {
      return ContainerOption::verletListsCells;
    }
  }

 private:

  void generateRebuildPositionMap() {

    if (!this->_partialRebuilding) {
      // every particle's position needs to be updated
      _particlePtr2rebuildPositionBuffer.clear();
      _particlePtr2rebuildPositionBuffer.reserve(this->_neighborList.getNumberOfParticles());

      for (auto iter = this->begin(IteratorBehavior::owned); iter.isValid(); ++iter) {
          std::pair<Particle *, std::array<double, 3>> particlePositionPair = std::make_pair(&(*iter), (*iter).getR());
          _particlePtr2rebuildPositionBuffer.emplace_back(particlePositionPair);
      }
    }
    else {

      for (FullParticleCell<Particle> &cell : this->_linkedCells.getCells()) {
          if (cell.getDirty()) {
            // TODO : update every particles position in the rebuild position map
            cell.setDirty(false);
          }
      }
    }
  }

  std::vector<std::pair<Particle*, std::array<double, 3>>> _particlePtr2rebuildPositionBuffer;


};

}
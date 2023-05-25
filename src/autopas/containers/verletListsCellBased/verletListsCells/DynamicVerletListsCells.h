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
    for (auto iter = _particlePtr2rebuildPositionBuffer.begin(); iter != _particlePtr2rebuildPositionBuffer.end(); ++iter) {
      auto particlePositionPair = (*iter);
      auto distance = utils::ArrayMath::sub(particlePositionPair.first->getR(), particlePositionPair.second);
      double distanceSquare = utils::ArrayMath::dot(distance, distance);

      if (distanceSquare >=  halfSkinSquare) {
        listInvalid = true;
        auto & cell = this->_linkedCells.getCellBlock().getContainingCell(particlePositionPair.first->getR());
        // TODO : figure out why this is not enough synchronization because of atomic bool and critical section
#pragma omp critical
        { cell.setDirty(true); }
        partialRebuilding = true;
        // delete this entry from the buffer because it is reinserted later -> avoid duplicates
        iter = _particlePtr2rebuildPositionBuffer.erase(iter);
        --iter;
        }
    }

    if (partialRebuilding) {
        for (size_t i = 0; i < this->_linkedCells.getCells().size(); ++i) {
          if (this->_linkedCells.getCells().at(i).getDirty()) {
            std::cout << "Cell with index : " << i << " is dirty \n" << std::endl;
          }
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
            // reinsert every "dirty" particle back in the rebuildPositionBuffer
            for (auto iter = cell.begin(); iter != cell.end(); ++iter) {
              std::pair<Particle *, std::array<double, 3>> particlePositionPair = std::make_pair(&(*iter), (*iter).getR());
              _particlePtr2rebuildPositionBuffer.emplace_back(particlePositionPair);
            }
            cell.setDirty(false);
          }
      }
    }
  }

  std::vector<std::pair<Particle*, std::array<double, 3>>> _particlePtr2rebuildPositionBuffer;


};

}
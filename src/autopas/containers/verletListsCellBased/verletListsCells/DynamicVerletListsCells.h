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
        cell.setDirty(true);
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
      _particlePtr2BufferPosition.clear();
      _particlePtr2BufferPosition.reserve(this->_neighborList.getNumberOfParticles());

      size_t particleBufferPosition = 0;
      for (auto iter = this->begin(IteratorBehavior::owned); iter.isValid(); ++iter) {
          std::pair<Particle *, std::array<double, 3>> particlePositionPair = std::make_pair(&(*iter), (*iter).getR());
          _particlePtr2rebuildPositionBuffer.emplace_back(particlePositionPair);
          std::pair<Particle *, size_t> particleMapping = std::make_pair(&(*iter), particleBufferPosition);
          _particlePtr2BufferPosition.emplace(particleMapping);
          ++particleBufferPosition;
      }
    }
    else {

      for (FullParticleCell<Particle> &cell : this->_linkedCells.getCells()) {
          if (cell.getDirty()) {
            // update every "dirty" particle in the rebuildPositionBuffer
            for (auto iter = cell.begin(); iter != cell.end(); ++iter) {
              std::pair<Particle *, std::array<double, 3>> particlePositionPair = std::make_pair(&(*iter), (*iter).getR());
              if (_particlePtr2BufferPosition.find(&(*iter)) != _particlePtr2BufferPosition.end()) {
                auto particleBufferPosition = _particlePtr2BufferPosition.at(&(*iter));
                _particlePtr2rebuildPositionBuffer.at(particleBufferPosition) = particlePositionPair;
              }
              else {
                _particlePtr2rebuildPositionBuffer.emplace_back(particlePositionPair);
                _particlePtr2BufferPosition.emplace(
                   std::make_pair(&(*iter), _particlePtr2rebuildPositionBuffer.size() - 1));
              }

            }
            cell.setDirty(false);
          }
      }
    }
  }

  std::vector<std::pair<Particle*, std::array<double, 3>>> _particlePtr2rebuildPositionBuffer;
  std::unordered_map<Particle*, size_t> _particlePtr2BufferPosition;

};

}
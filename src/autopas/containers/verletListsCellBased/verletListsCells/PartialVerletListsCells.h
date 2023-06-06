/**
* @file PartialVerletListsCells.h
* @author Luis Gall
* @date 31.05.23
 */

#pragma once

#include "autopas/containers/verletListsCellBased/verletListsCells/DynamicVerletListsCells.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCPartialCellPairNeighborList.h"

namespace autopas {

template <class Particle, class NeighborList>
class PartialVerletListsCells : public DynamicVerletListsCells<Particle, NeighborList> {

 public:

  PartialVerletListsCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                          const double skinPerTimestep = 0, const unsigned int rebuildFrequency = 2,
                          const double cellSizeFactor = 1.0,
                          const LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell,
                          typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType =
                              VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild)
      : DynamicVerletListsCells<Particle, NeighborList>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency, cellSizeFactor, loadEstimator, buildType) {}

  bool neighborListsAreValid() override {

    auto halfSkinSquare = (this->getVerletSkin() * this->getVerletSkin()) / 4;
    bool listInvalid = false;
    auto partialRebuilding = false;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for reduction(|| : listInvalid, partialRebuilding) schedule(static, 50)
#endif
    for (auto iter = this->_particlePtr2rebuildPositionBuffer.begin(); iter != this->_particlePtr2rebuildPositionBuffer.end(); ++iter) {
      auto particlePositionPair = (*iter);

      if (particlePositionPair.first->isDummy() || particlePositionPair.first->isHalo()) {
        continue;
      }

      auto distance = utils::ArrayMath::sub(particlePositionPair.first->getR(), particlePositionPair.second);
      double distanceSquare = utils::ArrayMath::dot(distance, distance);

      if (distanceSquare >=  halfSkinSquare) {
        listInvalid = true;
        auto & cell = this->_linkedCells.getCellBlock().getContainingCell(particlePositionPair.first->getR());
        cell.setDirty(true);
        partialRebuilding = true;
      }
    }

    size_t numHighMovementCells = std::count_if(this->_linkedCells.getCells().begin(), this->_linkedCells.getCells().end(), [](auto & cell) { return cell.getDirty(); });
    this->_cellMovingCounter.emplace_back(numHighMovementCells);

    this->_partialRebuilding = partialRebuilding;
    return !listInvalid;
  }

  [[nodiscard]] ContainerOption getContainerType() const override {
    return ContainerOption::partialPairwiseVerletLists;
  }

 protected:

  void generateRebuildPositionMap() override {

    size_t numRebuildCells {0};
    size_t numInflowCells{0};
    size_t numOutflowCells{0};

    if (!this->_partialRebuilding) {
      // every particle's position needs to be updated
      this->_particlePtr2rebuildPositionBuffer.clear();
      this->_particlePtr2rebuildPositionBuffer.reserve(this->_neighborList.getNumberOfParticles());
      _particleId2BufferPosition.clear();
      _particleId2BufferPosition.reserve(this->_neighborList.getNumberOfParticles());

      size_t particleBufferPosition = 0;
      for (auto iter = this->begin(IteratorBehavior::owned); iter.isValid(); ++iter) {
        std::pair<Particle *, std::array<double, 3>> particlePositionPair = std::make_pair(&(*iter), (*iter).getR());
        this->_particlePtr2rebuildPositionBuffer.emplace_back(particlePositionPair);
        std::pair<unsigned long, size_t> particleMapping = std::make_pair((*iter).getID(), particleBufferPosition);
        _particleId2BufferPosition.emplace(particleMapping);
        ++particleBufferPosition;
      }
    }
    else {
      // TODO : think about chunk size
#pragma omp parallel for reduction (+:numRebuildCells, numInflowCells, numOutflowCells) schedule (dynamic, 5)
      for (FullParticleCell<Particle> &cell : this->_linkedCells.getCells()) {
        if (cell.getDirty() || cell.getInflowDirty() /* || cell.getOutflowDirty()*/) {
          ++numRebuildCells;

          if (cell.getInflowDirty()) {
            ++numInflowCells;
          }

          // update every "dirty" particle in the rebuildPositionBuffer
          for (auto iter = cell.begin(); iter != cell.end(); ++iter) {

            if ((*iter).isDummy() || (*iter).isHalo()) {
              continue;
            }

            std::pair<Particle *, std::array<double, 3>> particlePositionPair = std::make_pair(&(*iter), (*iter).getR());
            if (_particleId2BufferPosition.find((*iter).getID()) != _particleId2BufferPosition.end()) {
              auto particleBufferPosition = _particleId2BufferPosition.at((*iter).getID());
              this->_particlePtr2rebuildPositionBuffer.at(particleBufferPosition) = particlePositionPair;
            }
            else {
#pragma omp critical
              {
                this->_particlePtr2rebuildPositionBuffer.emplace_back(particlePositionPair);
                _particleId2BufferPosition.emplace(
                    std::make_pair((*iter).getID(), this->_particlePtr2rebuildPositionBuffer.size() - 1));
              }
            }
          }
          cell.setDirty(false);
          cell.setInflowDirty(false);
        }
        if (cell.getOutflowDirty()) {
          ++numOutflowCells;
          cell.setOutflowDirty(false);
        }
      }
    }
    this->_cellOutflowCounter.emplace_back(numOutflowCells);
    this->_cellRebuildCounter.emplace_back(numRebuildCells);
    this->_cellInflowCounter.emplace_back(numInflowCells);
  }

 private:
  // TODO : think about how to determine particle's id type with template parameter
  std::unordered_map<unsigned long, size_t> _particleId2BufferPosition;
};

}
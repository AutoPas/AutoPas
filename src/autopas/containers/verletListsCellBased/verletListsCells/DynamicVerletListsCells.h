/**
 * @file DynamicVerletListsCells.h
 * @author Luis Gall
 * @date 06.05.23
 */

#pragma once

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"

namespace autopas {

/**
 * Dynamic variant of the VerletListsCells container
 * Differs only in the way the necessity of the neighbor lists rebuilds is determined.
 * @tparam Particle
 * @tparam NeighborList The neighbor list used by this container.
 */
template <class Particle, class NeighborList>
class DynamicVerletListsCells : public VerletListsCells<Particle, NeighborList> {
 public:
  /**
   * Constructor of the VerletListsCells class.
   * The neighbor lists are build using a search radius of cutoff + skin*rebuildfrequency.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param rebuildFrequency the rebuild Frequency
   * @param skinPerTimestep the skin radius per Timestep
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param loadEstimator load estimation algorithm for balanced traversals
   * @param buildType data layout of the particles which are used to generate the neighbor lists
   */
  DynamicVerletListsCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                          const double skinPerTimestep = 0, const unsigned int rebuildFrequency = 2,
                          const double cellSizeFactor = 1.0,
                          const LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell,
                          typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType =
                              VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild)
      : VerletListsCells<Particle, NeighborList>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency,
                                                 cellSizeFactor, loadEstimator, buildType) {}

  bool neighborListsAreValid() override {
    auto halfSkinSquare = (this->getVerletSkin() * this->getVerletSkin()) / 4;
    bool listInvalid = false;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for reduction(|| : listInvalid) schedule(static, 50)
#endif
    for (auto &particlePositionPair : _particlePtr2rebuildPositionBuffer) {
      auto distance = utils::ArrayMath::sub(particlePositionPair.first->getR(), particlePositionPair.second);
      double distanceSquare = utils::ArrayMath::dot(distance, distance);

      if (distanceSquare >= halfSkinSquare) {
        listInvalid = true;
      }
    }

    return !listInvalid;
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    generateRebuildPositionMap();
    VerletListsCells<Particle, NeighborList>::rebuildNeighborLists(traversal);
  }

  [[nodiscard]] ContainerOption getContainerType() const override {
    if (this->_neighborList.getContainerType() == ContainerOption::pairwiseVerletLists) {
      return ContainerOption::dynamicPairwiseVerletLists;
    } else if (this->_neighborList.getContainerType() == ContainerOption::verletListsCells) {
      return ContainerOption::dynamicVerletListsCells;
    } else {
      return ContainerOption::verletListsCells;
    }
  }

 private:
  void generateRebuildPositionMap() {
    _particlePtr2rebuildPositionBuffer.clear();
    _particlePtr2rebuildPositionBuffer.reserve(this->_neighborList.getNumberOfParticles());

    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      std::pair<Particle *, std::array<double, 3>> particlePositionPair = std::make_pair(&(*iter), (*iter).getR());
      _particlePtr2rebuildPositionBuffer.emplace_back(particlePositionPair);
    }
  }

  std::vector<std::pair<Particle *, std::array<double, 3>>> _particlePtr2rebuildPositionBuffer;
};

}  // namespace autopas
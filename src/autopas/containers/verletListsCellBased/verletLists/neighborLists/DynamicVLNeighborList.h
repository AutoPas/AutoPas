/**
 * @file DynamicVLNeighborList.h
 * @author Luis Gall
 * @date 04.05.2023
 *
 * oriented on
 * @file VerletListsCells.h
 * @author nguyen
 * @date 30.08.18
*/

#pragma once

#include "DynamicVLGeneratorFunctor.h"
#include "VLNeighborListsInterface.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

template <class Particle>
class DynamicVLNeighborList : public VLNeighborListInterface<Particle> {

 public:

  using listType = typename NewVerletListHelpers<Particle>::DynamicNeighborListsType;
  using LinkedParticleCell = FullParticleCell<Particle>;
  /**
   * @copydoc VLNeighborListInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::dynamicVerletLists; }

  bool neighborListsAreValid() override {

    auto halfSkinSquare = (this->_internalLinkedCells->getVerletSkin()
                           * this->_internalLinkedCells->getVerletSkin()) / 4;

    bool listInvalid = false;
    size_t buckets = _aosNeighborList.bucket_count();

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for reduction(|| : listInvalid) schedule(dynamic)
#endif
    for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
      auto endIter = _aosNeighborList.end(bucketId);
      for (auto bucketIter = _aosNeighborList.begin(bucketId); bucketIter != endIter; ++bucketIter) {
        Particle &particle = *(bucketIter->first);

        auto distance = utils::ArrayMath::sub(particle.getR(), bucketIter->second.second);
        auto distanceSquare = utils::ArrayMath::dot(distance, distance);

        if (distanceSquare >= halfSkinSquare) {
          listInvalid = true;
        }
      }
    }

    return !listInvalid;

  }

  void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double interactionLength,
                            typename NewVerletListHelpers<Particle>::VLBuildType::Value buildType) override {

    _aosNeighborList.clear();
    this->_internalLinkedCells = &linkedCells;

    for (auto iter = linkedCells.begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      _aosNeighborList[&(*iter)] = std::make_pair(std::vector<Particle*>{}, (*iter).getR());
    }

    applyBuildFunctor(linkedCells, useNewton3, interactionLength, buildType);
  }

  /**
   * @copydoc VLNeighborListsInterface::getNumberOfPartners()
   */
  size_t getNumberOfPartners(const Particle *particle) const override {
    return _aosNeighborList.at(const_cast<Particle*>(particle)).first.size();
  }

  /**
   * Returns the neighbor list in AoS layout.
   * @return Neighbor list in AoS layout.
   */
  typename NewVerletListHelpers<Particle>::DynamicNeighborListsType &getAoSNeighborList() { return _aosNeighborList; }

  /**
   * Returns the neighbor list in SoA layout.
   * @return Neighbor list in SoA layout.
   */
  auto &getSoANeighborList() { return _soaNeighborLists; }

  void generateSoAFromAoS(LinkedCells<Particle> &linkedCells) override {

    // TODO : implement

  }

  void setUpTraversal(TraversalInterface *traversal) override {

    auto vlTraversal = dynamic_cast<VLTraversalInterface<LinkedParticleCell> *>(traversal);

    if (vlTraversal) {
      vlTraversal->setCellsAndNeighborLists(this->_internalLinkedCells->getCells(), _aosNeighborList,
                                            _soaNeighborLists);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in DynamicVLNeighborList.h. TraversalID: {}",
          traversal->getTraversalType());
    }

  }

 protected:

  void applyBuildFunctor(LinkedCells<Particle> & linkedCells, bool useNewton3, double interactionLength,
                         typename NewVerletListHelpers<Particle>::VLBuildType::Value buildType) override {

    DynamicVLGeneratorFunctor<Particle> f (_aosNeighborList, interactionLength);

    /// @todo autotune traversal
    switch (buildType) {
      case NewVerletListHelpers<Particle>::VLBuildType::aosBuild: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal =
              LCC08Traversal<LinkedParticleCell, DynamicVLGeneratorFunctor<Particle>,
                             DataLayoutOption::aos, theBool>(
                  linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, interactionLength,
                  linkedCells.getCellBlock().getCellLength());
          linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      case NewVerletListHelpers<Particle>::VLBuildType::soaBuild: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal =
              LCC08Traversal<LinkedParticleCell, DynamicVLGeneratorFunctor<Particle>,
                             DataLayoutOption::soa, theBool>(
                  linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, interactionLength,
                  linkedCells.getCellBlock().getCellLength());
          linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      default:
        utils::ExceptionHandler::exception("DynamicVLNeighborList::buildAoSNeighborList(): unsupported VLBuildType: {}",
                                           buildType);
        break;
    }
  }

 private:

  typename NewVerletListHelpers<Particle>::DynamicNeighborListsType _aosNeighborList{};

  std::unordered_map<const Particle *, size_t> _particlePtr2indexMap;

  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> _soaNeighborLists;
};
}
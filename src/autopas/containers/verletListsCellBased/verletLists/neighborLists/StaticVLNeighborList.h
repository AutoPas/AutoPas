/**
 * @file StaticVLNeighborList.h
 * @author Luis Gall
 * @date 04.05.2023
 *
 * oriented on
 * @file VerletListsCells.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "StaticVLGeneratorFunctor.h"
#include "VLNeighborListsInterface.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

template <class Particle>
class StaticVLNeighborList : public VLNeighborListInterface<Particle> {

 public:

  using listType = typename NewVerletListHelpers<Particle>::StaticNeighborListsType;
  using LinkedParticleCell = FullParticleCell<Particle>;
  /**
   * @copydoc VLNeighborListInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletLists; }

  bool neighborListsAreValid() override {
    return true;
  }

  void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double interactionLength,
                            typename NewVerletListHelpers<Particle>::VLBuildType::Value buildType) override {

    _aosNeighborList.clear();
    this->_internalLinkedCells = &linkedCells;

    for (auto iter = linkedCells.begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      _aosNeighborList[&(*iter)];
    }

    applyBuildFunctor(linkedCells, useNewton3, interactionLength, buildType);
  }

  /**
   * @copydoc VLNeighborListsInterface::getNumberOfPartners()
   */
  size_t getNumberOfPartners(const Particle *particle) const override {
    return _aosNeighborList.at(const_cast<Particle *>(particle)).size();
  }

  /**
   * Returns the neighbor list in AoS layout.
   * @return Neighbor list in AoS layout.
   */
  typename NewVerletListHelpers<Particle>::StaticNeighborListsType &getAoSNeighborList() { return _aosNeighborList; }

  /**
   * Returns the neighbor list in SoA layout.
   * @return Neighbor list in SoA layout.
   */
  auto &getSoANeighborList() { return _soaNeighborLists; }

  void generateSoAFromAoS(LinkedCells<Particle> &linkedCells) override {

    // TODO : implement

  }

  void setUpTraversal(TraversalInterface *traversal) override {

    auto vlTraversal = dynamic_cast<VLTraversalInterface<LinkedParticleCell>*>(traversal);

    if (vlTraversal) {
      vlTraversal->setCellsAndNeighborLists(this->_internalLinkedCells->getCells(), _aosNeighborList,
                                            _soaNeighborLists);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in StaticVLNeighborList.h. TraversalID: {}",
          traversal->getTraversalType());
    }

  }

 protected:

  void applyBuildFunctor(LinkedCells<Particle> & linkedCells, bool useNewton3, double interactionLength,
                         typename NewVerletListHelpers<Particle>::VLBuildType::Value buildType) override {

    StaticVLGeneratorFunctor<Particle> f (_aosNeighborList, interactionLength);

    /// @todo autotune traversal
    switch (buildType) {
      case NewVerletListHelpers<Particle>::VLBuildType::aosBuild: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal =
              LCC08Traversal<LinkedParticleCell, StaticVLGeneratorFunctor<Particle>,
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
              LCC08Traversal<LinkedParticleCell, StaticVLGeneratorFunctor<Particle>,
                             DataLayoutOption::soa, theBool>(
                  linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, interactionLength,
                  linkedCells.getCellBlock().getCellLength());
          linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      default:
        utils::ExceptionHandler::exception("StaticVLNeighborList::buildAoSNeighborList(): unsupported VLBuildType: {}",
                                           buildType);
        break;
    }
  }

 private:

  typename NewVerletListHelpers<Particle>::StaticNeighborListsType _aosNeighborList{};

  std::unordered_map<const Particle *, size_t> _particlePtr2indexMap;

  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> _soaNeighborLists;
};
}
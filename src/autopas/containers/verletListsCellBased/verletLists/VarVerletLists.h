/**
 * @file VarVerletLists.h
 * @author humig
 * @date 20.05.19
 */

#pragma once

#include "VerletLists.h"
#include "traversals/VarVerletTraversalInterface.h"

namespace autopas {

/**
 * Variable Verlet Lists container with different neighbor lists.

 * @tparam Particle The particle type this container contains.
 * @tparam NeighborList The Neighbor List this Verlet Container uses.
 */
template <class Particle, class NeighborList>
class VarVerletLists
    : public VerletListsLinkedBase<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                   typename VerletListHelpers<Particle>::SoAArraysType> {
  typedef typename VerletListHelpers<Particle>::SoAArraysType SoAArraysType;
  typedef typename VerletListHelpers<Particle>::VerletListParticleCellType LinkedParticleCell;

 public:
  /**
   * Constructor of the Variable VerletLists class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius.
   * @param cellSizeFactor cell size factor ralative to cutoff
   */
  VarVerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                 const double skin, const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle, LinkedParticleCell, SoAArraysType>(
            boxMin, boxMax, cutoff, skin, compatibleTraversals::allVarVLAsBuildCompatibleTraversals(), cellSizeFactor),
        _neighborList{} {}

  ContainerOption getContainerType() override { return _neighborList.getContainerType(); }

  /**
   * @copydoc VerletLists::iteratePairwise
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor *f, Traversal *traversal) {
    if (auto *traversalInterface = dynamic_cast<VarVerletTraversalInterface<NeighborList> *>(traversal)) {
      traversalInterface->initVerletTraversal(_neighborList);
      traversalInterface->iterateVerletLists(_neighborList);
      traversalInterface->endVerletTraversal(_neighborList);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VarVerletLists::iteratePairwise");
    }
  }

  /**
   * @copydoc VerletLists::checkNeighborListsAreValid
   */
  bool checkNeighborListsAreValid(bool useNewton3 = true) {
    // if a particle was added or deleted, ... the list is definitely invalid
    if (not this->_neighborListIsValid) {
      return false;
    }
    // if a particle moved more than skin/2 outside of its cell the list is
    // invalid
    if (this->isContainerUpdateNeeded()) {
      return false;
    }

    return _neighborList.checkNeighborListValidity(useNewton3, this->getCutoff());
  }

  /**
   * Returns the number of neighbor pairs in the list.
   * @return the number of neighbor pairs in the list.
   */
  long getNumberOfNeighborPairs() const { return _neighborList.getNumberOfNeighborPairs(); }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();
    _neighborList.buildNeighborList(this->_linkedCells, traversal->getUseNewton3());
    // the neighbor list is now valid
    this->_neighborListIsValid = true;

    if (traversal->getDataLayout() == DataLayoutOption::soa and not _neighborList.isSoAListValid()) {
      _neighborList.generateSoAFromAoS();
    }
  }

 private:
  NeighborList _neighborList;
};

}  // namespace autopas

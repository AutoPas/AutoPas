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

 * @tparam Particle
 * @tparam NeighborList The Neighbor List this Verlet Container uses.
 */
template <class Particle, class NeighborList>
class VarVerletLists
    : public VerletListsLinkedBase<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                   typename VerletListHelpers<Particle>::SoAArraysType> {
  typedef FullParticleCell<Particle> ParticleCell;
  typedef typename VerletListHelpers<Particle>::SoAArraysType SoAArraysType;
  typedef typename VerletListHelpers<Particle>::VerletListParticleCellType LinkedParticleCell;

 public:
  /**
   * @copydoc VerletLists::VerletLists
   *
   * @todo TODO Decide if buildVerletListType makes sense and implement it if it does
   */
  VarVerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                 const double skin, const unsigned int rebuildFrequency = 1,
                 const typename VerletLists<Particle>::BuildVerletListType buildVerletListType =
                     VerletLists<Particle>::BuildVerletListType::VerletSoA)
      : VerletListsLinkedBase<Particle, LinkedParticleCell, SoAArraysType>(boxMin, boxMax, cutoff, skin,
                                                                           rebuildFrequency, getAllTraversals()),
        _neighborList{} {}

  std::vector<TraversalOption> getAllTraversals() override { return _neighborList.getAllTraversals(); }

  ContainerOption getContainerType() override { return _neighborList.getContainerType(); }

  /**
   * @copydoc VerletLists::iteratePairwise
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor *f, Traversal *traversal, bool useNewton3Deprecated = true) {
    if (auto *traversalInterface =
            dynamic_cast<VarVerletTraversalInterface<FullParticleCell<Particle>, NeighborList> *>(traversal)) {
      if (this->needsRebuild()) {
        // TODO: See if newton3 type of the list fits to the traversal
        this->rebuildVerletLists(traversalInterface->usesNewton3());
      }

      if (traversalInterface->getDataLayout() == DataLayoutOption::soa and not _neighborList.isSoAListValid()) {
        _neighborList.generateSoAFromAoS();
      }

      traversalInterface->initVerletTraversal(_neighborList);
      traversalInterface->iterateVerletLists(_neighborList);
      traversalInterface->endVerletTraversal(_neighborList);

      this->_traversalsSinceLastRebuild++;
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VarVerletLists::iteratePairwise");
    }
  }

  /**
   * @todo implement
   * @param useNewton3
   * @return
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

    // TODO: Implement some validity checker functor like verlet lists

    return false;
  }

  long getNumberOfNeighborPairs() const { return _neighborList.getNumberOfNeighborPairs(); }

 protected:
  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in iteratePairwiseAoS() and iteratePairwiseSoA() appropriately!
   * @param useNewton3
   */
  void rebuildVerletLists(bool useNewton3 = true) {
    this->_verletBuiltNewton3 = useNewton3;
    _neighborList.buildNeighborList(this->_linkedCells, useNewton3);
    // TODO: Call neighbor list to build with soa type as parameter
    // the neighbor list is now valid
    this->_neighborListIsValid = true;
    this->_traversalsSinceLastRebuild = 0;
  }

 private:
  NeighborList _neighborList;
};

}  // namespace autopas

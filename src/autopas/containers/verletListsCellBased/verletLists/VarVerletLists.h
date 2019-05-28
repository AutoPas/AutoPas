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
template<class Particle, class NeighborList>
class VarVerletLists : public VerletLists<Particle> {

 public:

  /**
   * @copydoc VerletLists::VerletLists
   */
  VarVerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                 const double skin, const unsigned int rebuildFrequency = 1,
                 const typename VerletLists<Particle>::BuildVerletListType buildVerletListType =
                 VerletLists<Particle>::BuildVerletListType::VerletSoA)
      : VerletLists<Particle>(
      boxMin, boxMax, cutoff, skin, rebuildFrequency),
        _neighborList{} {}
  //TODO: Add applicable traversals to constructor call to verlet lists (and in that constructor)

  //TODO: Update this
  std::vector<TraversalOption> getAllTraversals() override { return _neighborList.getAllTraversals(); }

  ContainerOption getContainerType() override { return ContainerOption::verletLists; }

  /**
   * @copydoc VerletLists::iteratePairwise
   *
   * This method implements the pairwise iteration of this container,
   * but does not override VerletLists::iteratePairwise as that is not virtual.
   * // TODO: Find solution for this.
   */
  template<class ParticleFunctor, class Traversal>
  void iteratePairwiseVar(ParticleFunctor *f, Traversal *traversal, bool useNewton3Deprecated = true) {
    if (auto *traversalInterface = dynamic_cast<VarVerletTraversalInterface<NeighborList> *>(traversal)) {
      if (this->needsRebuild()) {
        //TODO: See if newton3 type of the list fits to the traversal
        this->rebuildVerletLists(traversalInterface->usesNewton3());
      }

      this->_traversalsSinceLastRebuild++;
      //TODO: Handle SoA stuff when it is implemented
      traversalInterface->iterateVerletLists(_neighborList);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VarVerletLists::iteratePairwiseVar");
    }
  }

 protected:

  /**
   * Update the verlet lists for AoS usage
   * @param useNewton3
   */
  void updateVerletListsAoS(bool useNewton3) override {
    VerletLists<Particle>::updateVerletListsAoS(useNewton3);
    _neighborList.buildNeighborList(this->_linkedCells, useNewton3);
    // TODO: Call neighbor list to build with soa type as parameter
  }

 private:
  NeighborList _neighborList;
};

}  // namespace autopas


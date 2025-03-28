/**
 * @file VarVerletLists.h
 * @author humig
 * @date 20.05.19
 */

#pragma once

#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/traversals/VVLTraversalInterface.h"

namespace autopas {

/**
 * Variable Verlet Lists container with different neighbor lists.
 * @tparam Particle_T The particle type this container contains.
 * @tparam NeighborList The Neighbor List this Verlet Container uses.
 */
template <class Particle_T, class NeighborList>
class VarVerletLists : public VerletListsLinkedBase<Particle_T> {
 public:
  /**
   * Constructor of the Variable VerletLists class.
   * The neighbor lists are build using a search radius of cutoff + skin*rebuildfrequency.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius per Timestep.
   * @param rebuildFrequency The rebuild Frequency.
   * @param cellSizeFactor cell size factor relative to cutoff
   */
  VarVerletLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
                 const double skin, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle_T>(boxMin, boxMax, cutoff, skin, rebuildFrequency,
                                          compatibleTraversals::allVarVLAsBuildCompatibleTraversals(), cellSizeFactor),
        _neighborList{} {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return _neighborList.getContainerType(); }

  void computeInteractions(TraversalInterface *traversal) override {
    auto *traversalInterface = dynamic_cast<VVLTraversalInterface<NeighborList> *>(traversal);
    if (traversalInterface) {
      traversalInterface->setNeighborListToTraverse(_neighborList);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VarVerletLists::computeInteractions");
    }

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  /**
   * Returns the number of neighbor pairs in the list.
   * @return the number of neighbor pairs in the list.
   */
  long getNumberOfNeighborPairs() const { return _neighborList.getNumberOfNeighborPairs(); }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();
    _neighborList.buildAoSNeighborList(this->_linkedCells, traversal->getUseNewton3());
    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);

    if (traversal->getDataLayout() == DataLayoutOption::soa and not _neighborList.isSoAListValid()) {
      _neighborList.generateSoAFromAoS();
    }
  }

 private:
  NeighborList _neighborList;
};

}  // namespace autopas

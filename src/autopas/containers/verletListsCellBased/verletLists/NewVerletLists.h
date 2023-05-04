/**
* @file NewVerletLists.h
* @author Luis Gall
* @date 04.05.20223
*/

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/NewVerletListHelpers.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/neighborLists/StaticVLNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletLists/neighborLists/DynamicVLNeighborList.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

template <class Particle, class NeighborList>
class NewVerletLists : public VerletListsLinkedBase<Particle> {
 using LinkedParticleCell = FullParticleCell<Particle>;

public:

 NewVerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                  const double skinPerTimestep = 0, const unsigned int rebuildFrequency = 2,
                  typename NewVerletListHelpers<Particle>::VLBuildType::Value buildType =
                      NewVerletListHelpers<Particle>::VLBuildType::soaBuild,
                  const double cellSizeFactor = 1.0)
     : VerletListsLinkedBase<Particle>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency,
                                       compatibleTraversals::allVLCompatibleTraversals(), cellSizeFactor),
       _buildType(buildType) {}

 /**
  * @copydoc ParticleContainerInterface::getContainerType()
  */
 [[nodiscard]] ContainerOption getContainerType() const override { return _neighborList.getContainerType(); }

 void iteratePairwise(TraversalInterface *traversal) override {

   _neighborList.setUpTraversal(traversal);

   traversal->initTraversal();
   traversal->traverseParticlePairs();
   traversal->endTraversal();
 }

 /**
  * Gets the number of neighbors over all neighbor lists that belong to this particle.
  * @param particle
  * @return the size of the neighbor list(s) of this particle
  */
 size_t getNumberOfPartners(const Particle *particle) const { return _neighborList.getNumberOfPartners(particle); }

 void rebuildNeighborLists(TraversalInterface *traversal) override {

   this->_verletBuiltNewton3 = traversal->getUseNewton3();

   _neighborList.buildAoSNeighborList(this->_linkedCells, this->_verletBuiltNewton3,
                                      this->getInteractionLength(), _buildType);

   if (traversal->getDataLayout() == DataLayoutOption::soa) {
     _neighborList.generateSoAFromAoS(this->_linkedCells);
   }

   this->_neighborListIsValid.store(true, std::memory_order_relaxed);
 }

 bool neighborListsAreValid() override { return _neighborList.neighborListsAreValid(); }

 /**
  * Return the cell length of the underlying linked cells structure, normally needed only for unit tests.
  * @return
  */
 [[nodiscard]] const std::array<double, 3> &getCellLength() const {
   return this->_linkedCells.getCellBlock().getCellLength();
 }

private:
 /**
  * Neighbor list abstraction for neighbor list used in the container.
  */
 NeighborList _neighborList;

 /**
  * Data layout of the particles which are used to generate the neighbor lists.
  */
 typename NewVerletListHelpers<Particle>::VLBuildType::Value _buildType;
};
}  // namespace autopas

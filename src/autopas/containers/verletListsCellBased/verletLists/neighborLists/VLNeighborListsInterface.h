/**
* @file VLNeighborListInterface.h
* @author Luis Gall
* @date 04.05.2023
*/

#pragma once

#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/verletLists/NewVerletListHelpers.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/ContainerOption.h"

namespace autopas {
/**
* Interface of neighbor lists to be used with VerletLists container.
* @tparam Particle Type of particle to be used for the neighbor list.
*/
template <class Particle>
class VLNeighborListInterface {
public:
 /**
  * Virtual default destructor.
  */
 virtual ~VLNeighborListInterface() = default;

 virtual void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double interactionLength,
                                   typename NewVerletListHelpers<Particle>::VLBuildType::Value buildType) = 0;

 /**
  * Gets the number of neighbors over all neighbor lists that belong to this particle.
  * @param particle
  * @return the size of the neighbor list(s) of this particle
  */
 virtual size_t getNumberOfPartners(const Particle *particle) const = 0;

 /**
   * Returns the container type of this neighbor list and the container it belongs to.
   * @return ContainerOption for this neighbor list and the container it belongs to.
  */
 [[nodiscard]] virtual ContainerOption getContainerType() const = 0;

 virtual bool neighborListsAreValid() = 0;

 /**
   * Generates neighbor list in SoA layout from available neighbor list in AoS layout.
   * Copies the structure of the AoS neighbor list and replaces the particle pointers with the global indices of the
   * particles.
   * @param linkedCells Underlying linked cells structure.
  */
 virtual void generateSoAFromAoS(LinkedCells<Particle> &linkedCells) = 0;

 /**
   * Loads cells into structure of arrays.
   * @tparam TFunctor
   * @param f Functor that handles the loading.
   * @return loaded structure of arrays
  */
 template <class TFunctor>
 auto *loadSoA(TFunctor *f) {
   _soa.clear();
   size_t offset = 0;
   for (auto &cell : _internalLinkedCells->getCells()) {
     f->SoALoader(cell, _soa, offset);
     offset += cell.numParticles();
   }
   return &_soa;
 }

 /**
   * Extracts cells from structure of arrays.
   * @tparam TFunctor
   * @param f Functor that handles the extraction.
  */
 template <class TFunctor>
 void extractSoA(TFunctor *f) {
   size_t offset = 0;
   for (auto &cell : _internalLinkedCells->getCells()) {
     f->SoAExtractor(cell, _soa, offset);
     offset += cell.numParticles();
   }
 }

 /**
  * Assigns the current traversal to the correct traversal interface. The choice of traversal interface depends on the
  * type(s) of neighbor list allowed to use the current traversal. Currently VLCCellPairTraversalInterface handles the
  * traversals allowed solely for VLCCellPairNeighborList, while VLCTraversalInterface handles the traversals allowed
  * for both VLCCellPairNeighborList and VLCAllCellsNeighborList.
  * @param traversal the current traversal
  */
 virtual void setUpTraversal(TraversalInterface *traversal) = 0;

protected:

 LinkedCells<Particle>* _internalLinkedCells;

 /**
  * Structure of arrays necessary for SoA data layout.
  */
 SoA<typename Particle::SoAArraysType> _soa;

 virtual void applyBuildFunctor(LinkedCells<Particle> &linkedCells, bool useNewton3,
                                double interactionLength, typename NewVerletListHelpers<Particle>::VLBuildType::Value buildType) = 0;
};

}  // namespace autopas

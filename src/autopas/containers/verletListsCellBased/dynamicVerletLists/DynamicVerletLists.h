/**
 * @file DynamicVerletLists.h
 * @author Luis Gall
 * @date 26.04.23
 */

#pragma once

#include "DynamicVerletListHelpers.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

/**
 * Dynamic Verlet Lists container
 * This container is an improvement for the general VerletLists container
 * @tparam Particle
 */
template <class Particle>
class DynamicVerletLists : public VerletListsLinkedBase<Particle> {
    using LinkedParticleCell = FullParticleCell<Particle>;

 public:

  enum BuildVerletListType {
    VerletAoS,
    VerletSoA
  };

  DynamicVerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
              const double skinPerTimestep, const unsigned int rebuildFrequency,
              const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency,
                                        compatibleTraversals::allVLCompatibleTraversals(), cellSizeFactor),
        _buildVerletListType(buildVerletListType) {}

 [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::dynamicVerletLists; };

 void iteratePairwise(TraversalInterface *traversal) override {
   // Check if traversal is allowed for this container and give it the data it needs.
   auto *verletTraversalInterface = dynamic_cast<VLTraversalInterface<LinkedParticleCell> *>(traversal);
   if (verletTraversalInterface) {
     typename VerletListHelpers<Particle>::NeighborListAoSType converted {};
     verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), converted,
                                                        _soaNeighborLists);
   } else {
     autopas::utils::ExceptionHandler::exception(
         "trying to use a traversal of wrong type in DynamicVerletLists::iteratePairwise");
   }

   traversal->initTraversal();
   traversal->traverseParticlePairs();
   traversal->endTraversal();
 }

 typename DynamicVerletListHelpers<Particle>::AoSNeighborListType &getDynamicVerletListsAoS() { return _aosNeighborLists; };

 void rebuildNeighborLists(TraversalInterface *traversal) override {
   this->_verletBuiltNewton3 = traversal->getUseNewton3();
   this->updateVerletListsAoS(traversal->getUseNewton3());
   this->_neighborListIsValid.store(true, std::memory_order_relaxed);

   if (not _soaListsIsValid and traversal->getDataLayout() == DataLayoutOption::soa) {
     // generateSoAListFromAoSVerletLists();
   }
 }

protected:

 void updateVerletListsAoS(bool useNewton3) {
   generateAoSNeighborLists();
   typename DynamicVerletListHelpers<Particle>::DynamicVerletListGeneratorFunctor f(_aosNeighborLists,
                                                                      this->getCutoff() + this->getVerletSkin());
   /// @todo autotune traversal
   switch (_buildVerletListType) {
     case BuildVerletListType::VerletAoS: {
       utils::withStaticBool(useNewton3, [&](auto theBool) {
         auto traversal =
             LCC08Traversal<LinkedParticleCell, typename DynamicVerletListHelpers<Particle>::DynamicVerletListGeneratorFunctor,
                            DataLayoutOption::aos, theBool>(
                 this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
                 this->_linkedCells.getCellBlock().getCellLength());
         this->_linkedCells.iteratePairwise(&traversal);
       });
       break;
     }
     case BuildVerletListType::VerletSoA: {
       utils::withStaticBool(useNewton3, [&](auto theBool) {
         auto traversal =
             LCC08Traversal<LinkedParticleCell, typename DynamicVerletListHelpers<Particle>::DynamicVerletListGeneratorFunctor,
                            DataLayoutOption::soa, theBool>(
                 this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
                 this->_linkedCells.getCellBlock().getCellLength());
         this->_linkedCells.iteratePairwise(&traversal);
       });
       break;
     }
     default:
       utils::ExceptionHandler::exception("DynamicVerletLists::updateVerletListsAoS(): unsupported BuildVerletListType: {}",
                                          _buildVerletListType);
       break;
   }

   _soaListsIsValid = false;
 }

 size_t generateAoSNeighborLists() {
   size_t numParticles = 0;
   _aosNeighborLists.clear();

   for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++numParticles) {
     _aosNeighborLists[&(*iter)];
   }

   return numParticles;
 }

 private:

  typename DynamicVerletListHelpers<Particle>::AoSNeighborListType _aosNeighborLists;

  std::unordered_map<const Particle *, size_t> _particlePtr2indexMap;

  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> _soaNeighborLists;

  bool _soaListsIsValid{false};

  BuildVerletListType _buildVerletListType;
};
}
/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LCC01Traversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08NeighborListBuilding3B.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise or triwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @tparam Particle_T
 */
template <class Particle_T>
class VerletLists : public VerletListsLinkedBase<Particle_T> {
  /**
   * Type of the Particle.
   */
  using ParticleType = Particle_T;
  /**
   * Type of the ParticleCell used by the underlying linked cells.
   */
  using ParticleCellType = FullParticleCell<Particle_T>;

 public:
  /**
   * Enum that specifies how the verlet lists should be build
   */
  enum BuildVerletListType {
    /**
     * Build it using AoS
     */
    VerletAoS,
    /**
     * Build it using AoS
     */
    VerletSoA,
  };

  /**
   * Constructor of the VerletLists class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius per timestep.
   * @param buildVerletListType Specifies how the verlet list should be build, see BuildVerletListType
   * @param cellSizeFactor cell size factor ralative to cutoff
   */
  VerletLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
              const double skin, const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle_T>(boxMin, boxMax, cutoff, skin, cellSizeFactor),
        _buildVerletListType(buildVerletListType) {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletLists; }

  void computeInteractions(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *verletTraversalInterface = dynamic_cast<VLTraversalInterface<ParticleCellType> *>(traversal);
    if (verletTraversalInterface) {
      verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), _aosNeighborLists,
                                                         _soaNeighborLists, _aosNeighborPairsLists);
    } else {
      utils::ExceptionHandler::exception("trying to use a traversal of wrong type in VerletLists::computeInteractions");
    }

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  /**
   * get the actual neighbor list
   * @return the neighbor list
   */
  typename VerletListHelpers<Particle_T>::NeighborListAoSType &getVerletListsAoS() { return _aosNeighborLists; }

  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in computeInteractions()
   * @param traversal
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();

    // Check for triwise traversals
    switch (traversal->getTraversalType()) {
      // Standard pairwise traversal
      case TraversalOption::vl_list_iteration: {
        this->updateVerletListsAoS<InteractionTypeOption::pairwise>(traversal->getUseNewton3());
        break;
      }
      case TraversalOption::vl_list_intersection_sorted_3b: {
        this->updateVerletListsAoS<InteractionTypeOption::triwise>(traversal->getUseNewton3());

        // sort neighborLists for efficient intersecting
        /// @todo paralelize sorting
        size_t buckets = _aosNeighborLists.bucket_count();

        AUTOPAS_OPENMP(parallel for schedule(dynamic))
        for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
          auto endIter = _aosNeighborLists.end(bucketId);
          for (auto bucketIter = _aosNeighborLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
            auto &neighborPtrList = bucketIter->second;
            std::sort(neighborPtrList.begin(), neighborPtrList.end());
          }
        }
        break;
      }
      case TraversalOption::vl_pair_list_iteration_3b: {
        // build 3Body verletLists through VLIteration traversal
        this->updatePairVerletListsAoS3B(traversal->getUseNewton3());
        break;
      }
      // Default builds normal neighbor lists including halo particles.
      default: {
        this->updateVerletListsAoS<InteractionTypeOption::triwise>(traversal->getUseNewton3());
      }
    }

    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);

    if (not _soaListIsValid and traversal->getDataLayout() == DataLayoutOption::soa) {
      // only do this if we need it, i.e., if we are using soa!
      generateSoAListFromAoSVerletLists();
    }
  }

 protected:
  /**
   * Update the verlet lists for AoS usage
   * @param useNewton3
   */
  template <InteractionTypeOption::Value interactionType>
  void updateVerletListsAoS(bool useNewton3) {
    generateAoSNeighborLists();
    typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctor f(_aosNeighborLists,
                                                                         this->getCutoff() + this->getVerletSkin());

    /// @todo autotune traversal
    DataLayoutOption dataLayout;
    if (_buildVerletListType == BuildVerletListType::VerletAoS) {
      dataLayout = DataLayoutOption::aos;
    } else if (_buildVerletListType == BuildVerletListType::VerletSoA) {
      dataLayout = DataLayoutOption::soa;
    } else {
      utils::ExceptionHandler::exception("VerletLists::updateVerletListsAoS(): unsupported BuildVerletListType: {}",
                                         _buildVerletListType);
    }

    if constexpr (interactionType == InteractionTypeOption::pairwise) {
      auto pairwiseTraversal =
          LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctor>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);
      this->_linkedCells.computeInteractions(&pairwiseTraversal);
    } else {
      auto triwiseTraversal =
          LCC08NeighborListBuilding3B<ParticleCellType,
                                              typename VerletListHelpers<Particle_T>::PairVerletListGeneratorFunctor>(
        &f, dataLayout, useNewton3);
    this->computeInteractions(&traversal);
  }

  /**
   * Clears and then generates the AoS neighbor lists.
   * The Id Map is used to map the id of a particle to the actual particle.
   * @return Number of particles in the container
   */
  size_t generateAoSNeighborLists() {
    size_t numParticles = 0;
    _aosNeighborLists.clear();
    // DON'T simply parallelize this loop!!! this needs modifications if you want to parallelize it!
    // We have to iterate also over dummy particles here to ensure a correct size of the arrays.
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++numParticles) {
      // create the verlet list entries for all particles
      _aosNeighborLists[&(*iter)];
    }

    return numParticles;
  }

  /**
   * Clears and then generates the AoS neighbor pairs lists.
   */
  void generateAoSNeighborPairsLists() {
    _aosNeighborPairsLists.clear();
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      // create the pair verlet list entries for all particles
      _aosNeighborPairsLists[&(*iter)];
    }
  }

  /**
   * Fills SoA neighbor list with particle indices.
   */
  void generateSoAListFromAoSVerletLists() {
    // resize the list to the size of the aos neighborlist
    _soaNeighborLists.resize(_aosNeighborLists.size());
    // clear the aos 2 soa map
    _particlePtr2indexMap.clear();

    _particlePtr2indexMap.reserve(_aosNeighborLists.size());
    size_t index = 0;

    // Here we have to iterate over all particles, as particles might be later on marked for deletion, and we cannot
    // differentiate them from particles already marked for deletion.
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++index) {
      // set the map
      _particlePtr2indexMap[&(*iter)] = index;
    }
    size_t accumulatedListSize = 0;
    for (const auto &[particlePtr, neighborPtrVector] : _aosNeighborLists) {
      accumulatedListSize += neighborPtrVector.size();
      size_t i_id = _particlePtr2indexMap[particlePtr];
      // each soa neighbor list should be of the same size as for aos
      _soaNeighborLists[i_id].resize(neighborPtrVector.size());
      size_t j = 0;
      for (auto &neighborPtr : neighborPtrVector) {
        _soaNeighborLists[i_id][j] = _particlePtr2indexMap[neighborPtr];
        j++;
      }
    }

    AutoPasLog(DEBUG,
               "VerletLists::generateSoAListFromAoSVerletLists: average verlet list "
               "size is {}",
               static_cast<double>(accumulatedListSize) / _aosNeighborLists.size());
    _soaListIsValid = true;
  }

 private:
  /**
   * Neighbor Lists: Map of particle pointers to vector of particle pointers.
   */
  typename VerletListHelpers<Particle_T>::NeighborListAoSType _aosNeighborLists;

  /**
   * Neighbor Pairs Lists: Map of particle pointers to vector of pairs of particle pointers. (To find triplets.)
   */
  typename VerletListHelpers<Particle_T>::NeighborPairsListAoSType _aosNeighborPairsLists;

  /**
   * Mapping of every particle, represented by its pointer, to an index.
   * The index indexes all particles in the container.
   */
  std::unordered_map<const Particle_T *, size_t> _particlePtr2indexMap;

  /**
   * verlet list for SoA:
   * For every Particle, identified via the _particlePtr2indexMap, a vector of its neighbor indices is stored.
   */
  std::vector<std::vector<size_t, AlignedAllocator<size_t>>> _soaNeighborLists;

  /**
   * Shows if the SoA neighbor list is currently valid.
   */
  bool _soaListIsValid{false};

  /**
   * Specifies for what data layout the verlet lists are build.
   */
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas

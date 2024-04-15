/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LCC01Traversal3B.h"
#include "autopas/containers/linkedCells/traversals/LCC08NeighborListBuilding3B.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal3B.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise
 * interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @tparam Particle
 */
template <class Particle>
class VerletLists : public VerletListsLinkedBase<Particle> {
  using LinkedParticleCell = FullParticleCell<Particle>;

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
   * The neighbor lists are build using a search radius of cutoff + skinPerTimestep*rebuildFrequency.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skinPerTimestep The skin radius per timestep.
   * @param rebuildFrequency rebuild fequency.
   * @param buildVerletListType Specifies how the verlet list should be build, see BuildVerletListType
   * @param cellSizeFactor cell size factor ralative to cutoff
   */
  VerletLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
              const double skinPerTimestep, const unsigned int rebuildFrequency,
              const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency,
                                        compatibleTraversals::allVLCompatibleTraversals(), cellSizeFactor),
        _buildVerletListType(buildVerletListType) {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletLists; }

  void iteratePairwise(TraversalInterface<InteractionTypeOption::pairwise> *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *verletTraversalInterface = dynamic_cast<VLTraversalInterface<LinkedParticleCell> *>(traversal);
    if (verletTraversalInterface) {
      verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), _aosNeighborLists,
                                                         _soaNeighborLists, _pairwiseAosNeighborLists);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VerletLists::computeInteractions");
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  void iterateTriwise(TraversalInterface<InteractionTypeOption::threeBody> *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *verletTraversalInterface = dynamic_cast<VLTraversalInterface<LinkedParticleCell> *>(traversal);
    if (verletTraversalInterface) {
      verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), _aosNeighborLists,
                                                         _soaNeighborLists, _pairwiseAosNeighborLists);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VerletLists::computeInteractions");
    }

    traversal->initTraversal();
    traversal->traverseParticleTriplets();
    traversal->endTraversal();
  }

  /**
   * get the actual neighbour list
   * @return the neighbour list
   */
  typename VerletListHelpers<Particle>::NeighborListAoSType &getVerletListsAoS() { return _aosNeighborLists; }

  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in iteratePairwise()!
   * @param traversal
   */
  void rebuildNeighborLists(TraversalInterface<InteractionTypeOption::pairwise> *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();
    this->updateVerletListsAoS(traversal->getUseNewton3());
    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);

    if (not _soaListIsValid and traversal->getDataLayout() == DataLayoutOption::soa) {
      // only do this if we need it, i.e., if we are using soa!
      generateSoAListFromAoSVerletLists();
    }
  }

  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in iterateTriwise()!
   * @param traversal
   */
  void rebuildNeighborLists(TraversalInterface<InteractionTypeOption::threeBody> *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();

    switch (traversal->getTraversalType()) {
      case TraversalOption::vl_list_intersection_sorted_3b: {
        this->updateVerletListsAoS3B(traversal->getUseNewton3());

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
        this->updatePairwiseVerletListsAoS3B(traversal->getUseNewton3());
        break;
      }
      default: {
        this->updateVerletListsAoS3B(traversal->getUseNewton3());
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
  virtual void updateVerletListsAoS(bool useNewton3) {
    generateAoSNeighborLists();
    typename VerletListHelpers<Particle>::VerletListGeneratorFunctor f(_aosNeighborLists,
                                                                       this->getCutoff() + this->getVerletSkin());

    /// @todo autotune traversal
    switch (_buildVerletListType) {
      case BuildVerletListType::VerletAoS: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal =
              LCC08Traversal<LinkedParticleCell, typename VerletListHelpers<Particle>::VerletListGeneratorFunctor,
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
              LCC08Traversal<LinkedParticleCell, typename VerletListHelpers<Particle>::VerletListGeneratorFunctor,
                             DataLayoutOption::soa, theBool>(
                  this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
                  this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      default:
        utils::ExceptionHandler::exception("VerletLists::updateVerletListsAoS(): unsupported BuildVerletListType: {}",
                                           _buildVerletListType);
        break;
    }

    _soaListIsValid = false;
  }

  /**
   * Update the verlet lists for AoS usage, also builds Halo neighborlists for 3Body-Traversals that require those
   * @param useNewton3
   */
  virtual void updateVerletListsAoS3B(bool useNewton3) {
    generateAoSNeighborLists();
    typename VerletListHelpers<Particle>::VerletListGeneratorFunctor f(_aosNeighborLists,
                                                                       this->getCutoff() + this->getVerletSkin());

    /// @todo autotune traversal
    switch (_buildVerletListType) {
      case BuildVerletListType::VerletAoS: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal = LCC08NeighborListBuilding3B<LinkedParticleCell,
                                                       typename VerletListHelpers<Particle>::VerletListGeneratorFunctor,
                                                       DataLayoutOption::aos, theBool>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      case BuildVerletListType::VerletSoA: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal = LCC08NeighborListBuilding3B<LinkedParticleCell,
                                                       typename VerletListHelpers<Particle>::VerletListGeneratorFunctor,
                                                       DataLayoutOption::soa, theBool>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      default:
        utils::ExceptionHandler::exception("VerletLists::updateVerletListsAoS3B(): unsupported BuildVerletListType: {}",
                                           _buildVerletListType);
        break;
    }

    _soaListIsValid = false;
  }

  /**
   * Update the pairwise verlet lists for AoS usage
   * @param useNewton3
   */
  virtual void updatePairwiseVerletListsAoS3B(bool useNewton3) {
    updateVerletListsAoS3B(false);
    generatePairwiseAoSNeighborLists();
    typename VerletListHelpers<Particle>::PairwiseVerletListGeneratorFunctor f(
        _pairwiseAosNeighborLists, this->getCutoff() + this->getVerletSkin());

    /// @todo autotune traversal
    switch (_buildVerletListType) {
      case BuildVerletListType::VerletAoS: {
        auto traversal =
            VLListIterationTraversal3B<LinkedParticleCell,
                                       typename VerletListHelpers<Particle>::PairwiseVerletListGeneratorFunctor,
                                       DataLayoutOption::aos, false>(&f);
        this->iterateTriwise(&traversal);
        break;
      }
      case BuildVerletListType::VerletSoA: {
        // there are no SoA 3-body traversals, so we print out a warning and use AoS Layout instead
        AutoPasLog(WARN, "Pairwise Verlet-Lists can currently only be built with AoS DataLayout, using that instead!");
        auto traversal =
            VLListIterationTraversal3B<LinkedParticleCell,
                                       typename VerletListHelpers<Particle>::PairwiseVerletListGeneratorFunctor,
                                       DataLayoutOption::aos, false>(&f);
        this->iterateTriwise(&traversal);
        break;
      }
      default:
        utils::ExceptionHandler::exception(
            "VerletLists::updatePairwiseVerletListsAoS3B(): unsupported BuildVerletListType: {}", _buildVerletListType);
        break;
    }
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
   * Clears and then generates the pairwise AoS neighbor lists.
   */
  void generatePairwiseAoSNeighborLists() {
    _pairwiseAosNeighborLists.clear();
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      // create the pairwise verlet list entries for all particles
      _pairwiseAosNeighborLists[&(*iter)];
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
  typename VerletListHelpers<Particle>::NeighborListAoSType _aosNeighborLists;

  /**
   * Pairwise Neighbor Lists: Map of particle pointers to vector of pairs of particle pointers.
   */
  typename VerletListHelpers<Particle>::PairwiseNeighborListAoSType _pairwiseAosNeighborLists;

  /**
   * Mapping of every particle, represented by its pointer, to an index.
   * The index indexes all particles in the container.
   */
  std::unordered_map<const Particle *, size_t> _particlePtr2indexMap;

  /**
   * verlet list for SoA:
   * For every Particle, identified via the _particlePtr2indexMap, a vector of its neighbor indices is stored.
   */
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> _soaNeighborLists;

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

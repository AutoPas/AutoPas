/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/C08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VerletTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticSelectorMacros.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise
 * interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @note This class does NOT work with RMM cells and is not intended to!
 * @tparam Particle
 * @todo deleting particles should also invalidate the verlet lists - should be
 * implemented somehow
 */
template <class Particle>
class VerletLists
    : public VerletListsLinkedBase<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                   typename VerletListHelpers<Particle>::SoAArraysType> {
  typedef VerletListHelpers<Particle> verlet_internal;
  typedef FullParticleCell<Particle> ParticleCell;
  typedef typename VerletListHelpers<Particle>::SoAArraysType SoAArraysType;
  typedef typename VerletListHelpers<Particle>::VerletListParticleCellType LinkedParticleCell;

 public:
  /**
   * Enum that specifies how the verlet lists should be build
   */
  enum BuildVerletListType {
    VerletAoS,  /// build it using AoS
    VerletSoA   /// build it using SoA
  };

  /**
   * Constructor of the VerletLists class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * The rebuildFrequency should be chosen, s.t. the particles do not move more
   * than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius.
   * @param rebuildFrequency Specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals.
   * @param buildVerletListType Specifies how the verlet list should be build, see BuildVerletListType
   * @param cellSizeFactor cell size factor ralative to cutoff
   */
  VerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
              const double skin, const unsigned int rebuildFrequency = 1,
              const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle, LinkedParticleCell, SoAArraysType>(
            boxMin, boxMax, cutoff, skin, rebuildFrequency, compatibleTraversals::allVLCompatibleTraversals(),
            cellSizeFactor),
        _soaListIsValid(false),
        _buildVerletListType(buildVerletListType) {}

  ContainerOption getContainerType() override { return ContainerOption::verletLists; }

  /**
   * @copydoc LinkedCells::iteratePairwise
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor *f, Traversal *traversal, bool useNewton3 = true) {
    if (this->needsRebuild()) {
      rebuildVerletLists(useNewton3);
    }

    if (auto *traversalInterface = dynamic_cast<VerletTraversalInterface<LinkedParticleCell> *>(traversal)) {
      if (not _soaListIsValid and traversalInterface->getDataLayout() == DataLayoutOption::soa) {
        // only do this if we need it, i.e., if we are using soa!
        generateSoAListFromAoSVerletLists();
      }
      traversalInterface->initTraversal(this->_linkedCells.getCells());
      traversalInterface->iterateVerletLists(_aosNeighborLists, _soaNeighborLists);
      traversalInterface->endTraversal(this->_linkedCells.getCells());
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VerletLists::iteratePairwise");
    }
    this->_traversalsSinceLastRebuild++;
  }

  /**
   * get the actual neighbour list
   * @return the neighbour list
   */
  typename verlet_internal::AoS_verletlist_storage_type &getVerletListsAoS() { return _aosNeighborLists; }

  /**
   * Checks whether the neighbor lists are valid.
   * A neighbor list is valid if all pairs of particles whose interaction should
   * be calculated are represented in the neighbor lists.
   * @param useNewton3 specified whether newton 3 should be used
   * @return whether the list is valid
   * @note This check involves pair-wise interaction checks and is thus
   * relatively costly.
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

    // particles can also simply be very close already:
    typename verlet_internal::template VerletListValidityCheckerFunctor<LinkedParticleCell> validityCheckerFunctor(
        _aosNeighborLists, ((this->getCutoff() - this->_skin) * (this->getCutoff() - this->_skin)));

    auto traversal =
        C08Traversal<LinkedParticleCell,
                     typename verlet_internal::template VerletListValidityCheckerFunctor<LinkedParticleCell>,
                     DataLayoutOption::aos, true>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(),
                                                  &validityCheckerFunctor);
    this->_linkedCells.iteratePairwise(&validityCheckerFunctor, &traversal, useNewton3);

    return validityCheckerFunctor.neighborlistsAreValid();
  }

 protected:
  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in iteratePairwiseAoS() and iteratePairwiseSoA() appropriately!
   * @param useNewton3
   */
  void rebuildVerletLists(bool useNewton3 = true) {
    this->_verletBuiltNewton3 = useNewton3;
    this->updateVerletListsAoS(useNewton3);
    // the neighbor list is now valid
    this->_neighborListIsValid = true;
    this->_traversalsSinceLastRebuild = 0;
  }

  /**
   * Update the verlet lists for AoS usage
   * @param useNewton3
   */
  virtual void updateVerletListsAoS(bool useNewton3) {
    updateIdMapAoS();
    typename verlet_internal::VerletListGeneratorFunctor f(_aosNeighborLists, this->getCutoff());

    /// @todo autotune traversal
    switch (_buildVerletListType) {
      case BuildVerletListType::VerletAoS: {
        AUTOPAS_WITH_STATIC_BOOL(useNewton3, {
          auto traversal = C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor,
                                        DataLayoutOption::aos, c_useNewton3>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          this->_linkedCells.iteratePairwise(&f, &traversal);
        })
        break;
      }
      case BuildVerletListType::VerletSoA: {
        AUTOPAS_WITH_STATIC_BOOL(useNewton3, {
          auto traversal = C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor,
                                        DataLayoutOption::soa, c_useNewton3>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          this->_linkedCells.iteratePairwise(&f, &traversal);
        })
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
   * update the AoS id maps.
   * The Id Map is used to map the id of a particle to the actual particle
   * @return
   */
  size_t updateIdMapAoS() {
    size_t i = 0;
    _aosNeighborLists.clear();
    // DON'T simply parallelize this loop!!! this needs modifications if you
    // want to parallelize it!
    for (auto iter = this->begin(); iter.isValid(); ++iter, ++i) {
      // create the verlet list entries for all particles
      _aosNeighborLists[&(*iter)];
    }

    return i;
  }

  /**
   * Converts the verlet list stored for AoS usage into one for SoA usage
   */
  void generateSoAListFromAoSVerletLists() {
    // resize the list to the size of the aos neighborlist
    _soaNeighborLists.resize(_aosNeighborLists.size());
    // clear the aos 2 soa map
    _aos2soaMap.clear();

    _aos2soaMap.reserve(_aosNeighborLists.size());
    size_t i = 0;
    for (auto iter = this->begin(); iter.isValid(); ++iter, ++i) {
      // set the map
      _aos2soaMap[&(*iter)] = i;
    }
    size_t accumulatedListSize = 0;
    for (auto &aosList : _aosNeighborLists) {
      accumulatedListSize += aosList.second.size();
      size_t i_id = _aos2soaMap[aosList.first];
      // each soa neighbor list should be of the same size as for aos
      _soaNeighborLists[i_id].resize(aosList.second.size());
      size_t j = 0;
      for (auto neighbor : aosList.second) {
        _soaNeighborLists[i_id][j] = _aos2soaMap.at(neighbor);
        j++;
      }
    }
    AutoPasLog(debug,
               "VerletLists::generateSoAListFromAoSVerletLists: average verlet list "
               "size is {}",
               static_cast<double>(accumulatedListSize) / _aosNeighborLists.size());
    _soaListIsValid = true;
  }

 private:
  /// verlet lists.
  typename verlet_internal::AoS_verletlist_storage_type _aosNeighborLists;

  /// map converting from the aos type index (Particle *) to the soa type index
  /// (continuous, size_t)
  std::unordered_map<Particle *, size_t> _aos2soaMap;

  /// verlet list for SoA:
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> _soaNeighborLists;

  // specifies if the SoA neighbor list is currently valid
  bool _soaListIsValid;

  /// specifies how the verlet lists are build
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas

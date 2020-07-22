/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
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
 * @todo deleting particles should also invalidate the verlet lists - should be implemented somehow
 */
template <class Particle>
class VerletLists
    : public VerletListsLinkedBase<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                   typename VerletListHelpers<Particle>::SoAArraysType> {
  using verlet_internal = VerletListHelpers<Particle>;
  using ParticleCell = FullParticleCell<Particle>;
  using SoAArraysType = typename VerletListHelpers<Particle>::SoAArraysType;
  using LinkedParticleCell = typename VerletListHelpers<Particle>::VerletListParticleCellType;

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
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius.
   * @param buildVerletListType Specifies how the verlet list should be build, see BuildVerletListType
   * @param cellSizeFactor cell size factor ralative to cutoff
   */
  VerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
              const double skin, const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle, LinkedParticleCell, SoAArraysType>(
            boxMin, boxMax, cutoff, skin, compatibleTraversals::allVLCompatibleTraversals(), cellSizeFactor),
        _soaListIsValid(false),
        _buildVerletListType(buildVerletListType) {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletLists; }

  void iteratePairwise(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *verletTraversalInterface = dynamic_cast<VLTraversalInterface<LinkedParticleCell> *>(traversal);
    if (verletTraversalInterface) {
      verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), _aosNeighborLists,
                                                         _soaNeighborLists);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VerletLists::iteratePairwise");
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  /**
   * get the actual neighbour list
   * @return the neighbour list
   */
  typename verlet_internal::AoS_verletlist_storage_type &getVerletListsAoS() { return _aosNeighborLists; }

  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in iteratePairwiseAoS() and iteratePairwiseSoA() appropriately!
   * @param traversal
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();
    this->updateVerletListsAoS(traversal->getUseNewton3());
    // the neighbor list is now valid
    this->_neighborListIsValid = true;

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
    updateIdMapAoS();
    typename verlet_internal::VerletListGeneratorFunctor f(_aosNeighborLists, this->getCutoff() + this->getSkin());

    /// @todo autotune traversal
    switch (_buildVerletListType) {
      case BuildVerletListType::VerletAoS: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal = LCC08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor,
                                          DataLayoutOption::aos, theBool>(
              this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
              this->_linkedCells.getCellBlock().getCellLength());
          this->_linkedCells.iteratePairwise(&traversal);
        });
        break;
      }
      case BuildVerletListType::VerletSoA: {
        utils::withStaticBool(useNewton3, [&](auto theBool) {
          auto traversal = LCC08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor,
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

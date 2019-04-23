/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

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
  using floatPrecision = typename Particle::ParticleFloatingPointType;

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
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals.
   * @param buildVerletListType specifies how the verlet list should be build, see BuildVerletListType
   */
  VerletLists(const std::array<floatPrecision, 3> boxMin, const std::array<floatPrecision, 3> boxMax,
              const floatPrecision cutoff, const floatPrecision skin, const unsigned int rebuildFrequency = 1,
              const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA)
      : VerletListsLinkedBase<Particle, LinkedParticleCell, SoAArraysType>(
            boxMin, boxMax, cutoff, skin, rebuildFrequency, allVLApplicableTraversals()),
        _soaListIsValid(false),
        _buildVerletListType(buildVerletListType) {}

  /**
   * Lists all traversal options applicable for the Verlet Lists container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOption>& allVLApplicableTraversals() {
    static const std::vector<TraversalOption> v{TraversalOption::verletTraversal};
    return v;
  }

  std::vector<TraversalOption> getAllTraversals() override { return allVLApplicableTraversals(); }

  ContainerOption getContainerType() override { return ContainerOption::verletLists; }

  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in iteratePairwiseAoS() and iteratePairwiseSoA() appropriately!
   * @param useNewton3
   */
  void rebuild(bool useNewton3 = true) {
    this->updateVerletListsAoS(useNewton3);
    // the neighbor list is now valid
    this->_neighborListIsValid = true;
    this->_traversalsSinceLastRebuild = 0;
  }

  /**
   * @copydoc LinkedCells::iteratePairwise
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    if (needsRebuild()) {
      rebuild(useNewton3);
    }

    if (auto* traversalInterface = dynamic_cast<VerletTraversalInterface<LinkedParticleCell>*>(traversal)) {
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
  typename verlet_internal::AoS_verletlist_storage_type& getVerletListsAoS() { return _aosNeighborLists; }

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

  TraversalSelector<ParticleCell> generateTraversalSelector() override {
    // @FIXME dummyTraversal is a workaround because this container does not yet use traversals like it should
    return TraversalSelector<ParticleCell>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo());
  }

  /**
   * specifies whether the neighbor lists need to be rebuild
   * @return true if the neighbor lists need to be rebuild, false otherwise
   */
  bool needsRebuild() {
    AutoPasLog(debug, "Neighborlist is valid: {}", this->_neighborListIsValid);
    // if the neighbor list is NOT valid or we have not rebuilt since this->_rebuildFrequency steps
    return (not this->_neighborListIsValid) or (this->_traversalsSinceLastRebuild >= this->_rebuildFrequency);
  }

 protected:
  /**
   * update the verlet lists for AoS usage
   * @param useNewton3 CURRENTLY NOT USED!
   * @todo Build verlet lists according to newton 3.
   */
  virtual void updateVerletListsAoS(bool useNewton3) {
    updateIdMapAoS();
    typename verlet_internal::VerletListGeneratorFunctor f(_aosNeighborLists, this->getCutoff());

    /// @todo autotune traversal
    switch (_buildVerletListType) {
      case BuildVerletListType::VerletAoS: {
        auto traversal = C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor,
                                      DataLayoutOption::aos, true>(
            this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
        this->_linkedCells.iteratePairwise(&f, &traversal);
        break;
      }
      case BuildVerletListType::VerletSoA: {
        auto traversal = C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor,
                                      DataLayoutOption::soa, true>(
            this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
        this->_linkedCells.iteratePairwise(&f, &traversal);
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
    i = 0;
    size_t accumulatedListSize = 0;
    for (auto& aosList : _aosNeighborLists) {
      accumulatedListSize += aosList.second.size();
      size_t i_id = _aos2soaMap[aosList.first];
      // each soa neighbor list should be of the same size as for aos
      _soaNeighborLists[i_id].resize(aosList.second.size());
      size_t j = 0;
      for (auto neighbor : aosList.second) {
        _soaNeighborLists[i_id][j] = _aos2soaMap.at(neighbor);
        j++;
      }
      i++;
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
  std::unordered_map<Particle*, size_t> _aos2soaMap;

  /// verlet list for SoA:
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> _soaNeighborLists;

  // specifies if the SoA neighbor list is currently valid
  bool _soaListIsValid;

  /// specifies how the verlet lists are build
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas

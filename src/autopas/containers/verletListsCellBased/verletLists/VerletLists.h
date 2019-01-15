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

 private:
  static const std::vector<TraversalOptions>& VLApplicableTraversals() {
    // @todo: implement some traversals for this
    static const std::vector<TraversalOptions> v{};
    return v;
  }

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
  VerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
              const double skin, const unsigned int rebuildFrequency = 1,
              const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA)
      : VerletListsLinkedBase<Particle, LinkedParticleCell, SoAArraysType>(
            boxMin, boxMax, cutoff, skin, rebuildFrequency, allVLApplicableTraversals()),
        _soaListIsValid(false),
        _soa(),
        _buildVerletListType(buildVerletListType) {}

  /**
   * Lists all traversal options applicable for the Verlet Lists container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOptions>& allVLApplicableTraversals() {
    // @FIXME This is a workaround because this container does not yet use traversals like it should
    static const std::vector<TraversalOptions> v{TraversalOptions::dummyTraversal};
    return v;
  }

  ContainerOptions getContainerType() override { return ContainerOptions::verletLists; }

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
   * @copydoc LinkedCells::iteratePairwiseAoS
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseAoS(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    if (needsRebuild()) {  // if we need to rebuild the list, we should rebuild it!
      rebuild(useNewton3);
    }
    this->iterateVerletListsAoS(f, useNewton3);
    // we iterated, so increase traversal counter
    this->_traversalsSinceLastRebuild++;
  }

  /**
   * @copydoc LinkedCells::iteratePairwiseSoA
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseSoA(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    if (needsRebuild()) {
      rebuild(useNewton3);
      generateSoAListFromAoSVerletLists();
    } else if (not _soaListIsValid) {
      generateSoAListFromAoSVerletLists();
    }
    iterateVerletListsSoA(f, useNewton3);
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
                     typename verlet_internal::template VerletListValidityCheckerFunctor<LinkedParticleCell>, false,
                     true>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &validityCheckerFunctor);
    this->_linkedCells.iteratePairwiseAoS(&validityCheckerFunctor, &traversal, useNewton3);

    return validityCheckerFunctor.neighborlistsAreValid();
  }

  TraversalSelector<ParticleCell> generateTraversalSelector(std::vector<TraversalOptions> traversalOptions) override {
    //    std::vector<TraversalOptions> allowedAndApplicable;
    //
    //    std::sort(traversalOptions.begin(), traversalOptions.end());
    //    std::set_intersection(this->_applicableTraversals.begin(), this->_applicableTraversals.end(),
    //    traversalOptions.begin(),
    //                          traversalOptions.end(), std::back_inserter(allowedAndApplicable));
    // @FIXME dummyTraversal is a workaround because this container does not yet use traversals like it should
    return TraversalSelector<ParticleCell>({0, 0, 0}, {dummyTraversal});
  }

  /**
   * specifies whether the neighbor lists need to be rebuild
   * @return true if the neighbor lists need to be rebuild, false otherwise
   */
  bool needsRebuild() {
    AutoPasLog(debug, "Neighborlist is valid: {}", this->_neighborListIsValid);
    // if the neighbor list is NOT validor we have not rebuild for this->_rebuildFrequency steps
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
        auto traversal =
            C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, false, true>(
                this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
        this->_linkedCells.iteratePairwiseAoS(&f, &traversal);
        break;
      }
      case BuildVerletListType::VerletSoA: {
        auto traversal =
            C08Traversal<LinkedParticleCell, typename verlet_internal::VerletListGeneratorFunctor, true, true>(
                this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
        this->_linkedCells.iteratePairwiseSoA(&f, &traversal);
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
   * iterate over the verlet lists using the AoS traversal
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3
   */
  template <class ParticleFunctor>
  void iterateVerletListsAoS(ParticleFunctor* f, const bool useNewton3) {
    // @todo optimize iterateVerletListsAoS, e.g. by using openmp-capable
    /// traversals

#if defined(AUTOPAS_OPENMP)
    if (not useNewton3) {
      size_t buckets = _aosNeighborLists.bucket_count();
      // @todo find a sensible chunk size
#pragma omp parallel for schedule(dynamic)
      for (size_t b = 0; b < buckets; b++) {
        auto endIter = _aosNeighborLists.end(b);
        for (auto it = _aosNeighborLists.begin(b); it != endIter; ++it) {
          Particle& i = *(it->first);
          for (auto j_ptr : it->second) {
            Particle& j = *j_ptr;
            f->AoSFunctor(i, j, false);
          }
        }
      }
    } else
#endif
    {
      for (auto& list : _aosNeighborLists) {
        Particle& i = *list.first;
        for (auto j_ptr : list.second) {
          Particle& j = *j_ptr;
          f->AoSFunctor(i, j, useNewton3);
        }
      }
    }
  }

  /**
   * iterate over the verlet lists using the SoA traversal
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3
   */
  template <class ParticleFunctor>
  void iterateVerletListsSoA(ParticleFunctor* f, const bool useNewton3) {
    // @todo optimize iterateVerletListsSoA, e.g. by using traversals with
    /// openmp possibilities

    // load data from cells into soa
    loadVerletSoA(f);

    // @todo here you can (sort of) use traversals, by modifying iFrom and iTo.
    const size_t iFrom = 0;
    const size_t iTo = _soaNeighborLists.size();

#if defined(AUTOPAS_OPENMP)
    if (not useNewton3) {
      // @todo find a sensible chunk size
      const size_t chunkSize = std::max((iTo - iFrom) / (omp_get_max_threads() * 10), 1ul);
#pragma omp parallel for schedule(dynamic, chunkSize)
      for (size_t i = iFrom; i < iTo; i++) {
        f->SoAFunctor(_soa, _soaNeighborLists, i, i + 1, useNewton3);
      }
    } else
#endif
    {
      // iterate over SoA
      f->SoAFunctor(_soa, _soaNeighborLists, iFrom, iTo, useNewton3);
    }

    // extract SoA
    extractVerletSoA(f);
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
   * Load the particle information from the cell and store it in the global SoA
   * using functor.SoALoader(...)
   * @tparam ParticleFunctor the type of the functor
   * @param functor the SoAExtractor method of this functor is used. use the
   * actual
   */
  template <class ParticleFunctor>
  void loadVerletSoA(ParticleFunctor* functor) {
    size_t offset = 0;
    for (auto& cell : this->_linkedCells.getCells()) {
      functor->SoALoader(cell, _soa, offset);
      offset += cell.numParticles();
    }
  }

  /**
   * Extracts the particle information from the global SoA using
   * functor.SoAExtractor(...)
   * @tparam ParticleFunctor the type of the functor
   * @param functor the SoAExtractor method of this functor is used. use the
   * actual
   */
  template <class ParticleFunctor>
  void extractVerletSoA(ParticleFunctor* functor) {
    size_t offset = 0;
    for (auto& cell : this->_linkedCells.getCells()) {
      functor->SoAExtractor(cell, _soa, offset);
      offset += cell.numParticles();
    }
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

  /// global SoA of verlet lists
  SoA<typename Particle::SoAArraysType> _soa;

  /// specifies how the verlet lists are build
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas

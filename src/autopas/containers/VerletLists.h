/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "autopas/containers/LinkedCells.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/VerletListHelpers.h"
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
class VerletLists : public ParticleContainer<Particle, autopas::FullParticleCell<Particle>> {
  typedef VerletListHelpers<Particle> verlet_internal;
  typedef FullParticleCell<Particle> ParticleCell;

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
  VerletLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin,
              unsigned int rebuildFrequency = 1,
              BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA)
      : ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff + skin, allVLApplicableTraversals()),
        _linkedCells(boxMin, boxMax, cutoff + skin),
        _skin(skin),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _neighborListIsValid(false),
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
    _neighborListIsValid = true;
    _traversalsSinceLastRebuild = 0;
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
    _traversalsSinceLastRebuild++;
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
    _traversalsSinceLastRebuild++;
  }

  /**
   * get the actual neighbour list
   * @return the neighbour list
   */
  typename verlet_internal::AoS_verletlist_storage_type& getVerletListsAoS() { return _aosNeighborLists; }

  /**
   * @copydoc LinkedCells::addParticle()
   * @note This function invalidates the neighbor lists.
   */
  void addParticle(Particle& p) override {
    _neighborListIsValid = false;
    _linkedCells.addParticle(p);
  }

  /**
   * @copydoc LinkedCells::addHaloParticle()
   * @note This function invalidates the neighbor lists.
   */
  void addHaloParticle(Particle& haloParticle) override {
    _neighborListIsValid = false;
    _linkedCells.addHaloParticle(haloParticle);
  }

  unsigned long getNumParticles() override { return _linkedCells.getNumParticles(); }

  /**
   * @copydoc LinkedCells::deleteHaloParticles
   * @note This function invalidates the neighbor lists.
   */
  void deleteHaloParticles() override {
    _neighborListIsValid = false;
    _linkedCells.deleteHaloParticles();
  }

  /**
   * @copydoc ParticleContainerInterface::deleteAllParticles
   * @note This function invalidates the neighbor lists.
   */
  void deleteAllParticles() override {
    _neighborListIsValid = false;
    _linkedCells.deleteAllParticles();
  }

  /**
   * @copydoc LinkedCells::updateContainer()
   * @note This function invalidates the neighbor lists.
   */
  void updateContainer() override {
    AutoPasLog(debug, "Updating container");
    _neighborListIsValid = false;
    _linkedCells.updateContainer();
  }

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
    if (not _neighborListIsValid) {
      return false;
    }
    // if a particle moved more than skin/2 outside of its cell the list is
    // invalid
    if (this->isContainerUpdateNeeded()) {
      return false;
    }

    // particles can also simply be very close already:
    typename verlet_internal::template VerletListValidityCheckerFunctor<ParticleCell> validityCheckerFunctor(
        _aosNeighborLists, ((this->getCutoff() - _skin) * (this->getCutoff() - _skin)));

    auto traversal =
        C08Traversal<FullParticleCell<Particle, typename verlet_internal::SoAArraysType>,
                     typename verlet_internal::template VerletListValidityCheckerFunctor<ParticleCell>, false, true>(
            _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &validityCheckerFunctor);
    _linkedCells.iteratePairwiseAoS(&validityCheckerFunctor, &traversal, useNewton3);

    return validityCheckerFunctor.neighborlistsAreValid();
  }

  bool isContainerUpdateNeeded() override {
    std::atomic<bool> outlierFound(false);
#ifdef AUTOPAS_OPENMP
    // @todo: find a sensible value for ???
#pragma omp parallel for shared(outlierFound)  // if (this->_cells.size() / omp_get_max_threads() > ???)
#endif
    for (size_t cellIndex1d = 0; cellIndex1d < _linkedCells.getCells().size(); ++cellIndex1d) {
      std::array<double, 3> boxmin{0., 0., 0.};
      std::array<double, 3> boxmax{0., 0., 0.};
      _linkedCells.getCellBlock().getCellBoundingBox(cellIndex1d, boxmin, boxmax);
      boxmin = ArrayMath::addScalar(boxmin, -_skin / 2.);
      boxmax = ArrayMath::addScalar(boxmax, +_skin / 2.);
      for (auto iter = _linkedCells.getCells()[cellIndex1d].begin(); iter.isValid(); ++iter) {
        if (utils::notInBox(iter->getR(), boxmin, boxmax)) {
          outlierFound = true;  // we need an update
          break;
        }
      }
      if (outlierFound) cellIndex1d = _linkedCells.getCells().size();
    }
    if (outlierFound) {
      AutoPasLog(debug,
                 "VerletLists: containerUpdate needed! Particles are fast. You "
                 "might want to increase the skin radius or decrease the rebuild "
                 "frequency.");
    } else {
      AutoPasLog(debug,
                 "VerletLists: containerUpdate not yet needed. Particles are slow "
                 "enough.");
    }
    return outlierFound;
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
    AutoPasLog(debug, "Neighborlist is valid: {}", _neighborListIsValid);
    return (not _neighborListIsValid)                              // if the neighborlist is NOT valid a
                                                                   // rebuild is needed
           or (_traversalsSinceLastRebuild >= _rebuildFrequency);  // rebuild with frequency
  }

  /**
   * Searches the provided halo particle and updates the found particle.
   * Searches for the provided particle within the halo cells of the container
   * and overwrites the found particle with the provided particle.
   * @param particle
   */
  void updateHaloParticle(Particle& particle) {
    auto cells = _linkedCells.getCellBlock().getNearbyHaloCells(particle.getR(), _skin);
    bool updated = false;
    for (auto cellptr : cells) {
      updated |= checkParticleInCellAndUpdate(*cellptr, particle);
      if (updated) {
        continue;
      }
    }
    if (not updated) {
      AutoPasLog(error,
                 "VerletLists: updateHaloParticle was not able to update particle at "
                 "[{}, {}, {}]",
                 particle.getR()[0], particle.getR()[1], particle.getR()[2]);
      utils::ExceptionHandler::exception("VerletLists: updateHaloParticle could not find any particle");
    }
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return _linkedCells.begin(behavior);
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(std::array<double, 3> lowerCorner,
                                                      std::array<double, 3> higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    return _linkedCells.getRegionIterator(lowerCorner, higherCorner, behavior, true);
  }

 protected:
  /**
   * Updates a found particle within cellI to the values of particleI.
   * Checks whether a particle with the same id as particleI is within the cell
   * cellI and overwrites the particle with particleI, if it is found.
   * @param cellI
   * @param particleI
   * @return
   */
  bool checkParticleInCellAndUpdate(FullParticleCell<Particle, typename verlet_internal::SoAArraysType>& cellI,
                                    Particle& particleI) {
    for (auto iterator = cellI.begin(); iterator.isValid(); ++iterator) {
      if (iterator->getID() == particleI.getID()) {
        *iterator = particleI;
        return true;
      }
    }
    return false;
  }

  /**
   * update the verlet lists for AoS usage
   * @param useNewton3 CURRENTLY NOT USED!
   * @todo Build verlet lists according to newton 3.
   */
  virtual void updateVerletListsAoS(bool useNewton3) {
    updateIdMapAoS();
    typename verlet_internal::VerletListGeneratorFunctor f(_aosNeighborLists, this->getCutoff());

    // @todo autotune traversal
    switch (_buildVerletListType) {
      case BuildVerletListType::VerletAoS: {
        auto traversal = C08Traversal<FullParticleCell<Particle, typename verlet_internal::SoAArraysType>,
                                      typename verlet_internal::VerletListGeneratorFunctor, false, true>(
            _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
        _linkedCells.iteratePairwiseAoS(&f, &traversal);
        break;
      }
      case BuildVerletListType::VerletSoA: {
        auto traversal = C08Traversal<FullParticleCell<Particle, typename verlet_internal::SoAArraysType>,
                                      typename verlet_internal::VerletListGeneratorFunctor, true, true>(
            _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
        _linkedCells.iteratePairwiseSoA(&f, &traversal);
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
#pragma omp parallel for schedule(dynamic)
      for (size_t b = 0; b < buckets; b++) {
        for (auto it = _aosNeighborLists.begin(b); it != _aosNeighborLists.end(b); it++) {
          Particle& i = *it->first;
          for (auto j_ptr : it->second) {
            Particle& j = *j_ptr;
            f->AoSFunctor(i, j, useNewton3);
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
    size_t iFrom = 0;
    size_t iTo = _soaNeighborLists.size();

#if defined(AUTOPAS_OPENMP)
    if (not useNewton3) {
#pragma omp parallel for schedule(dynamic)
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
    for (auto& cell : _linkedCells.getCells()) {
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
    for (auto& cell : _linkedCells.getCells()) {
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

  /// internal linked cells storage, handles Particle storage and used to build verlet lists
  LinkedCells<Particle, FullParticleCell<Particle, typename verlet_internal::SoAArraysType>,
              typename verlet_internal::SoAArraysType>
      _linkedCells;

  /// map converting from the aos type index (Particle *) to the soa type index
  /// (continuous, size_t)
  std::unordered_map<Particle*, size_t> _aos2soaMap;

  /// verlet list for SoA:
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> _soaNeighborLists;

  /// skin radius
  double _skin;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  // specifies if the neighbor list is currently valid
  bool _neighborListIsValid;

  // specifies if the SoA neighbor list is currently valid
  bool _soaListIsValid;

  /// global SoA of verlet lists
  SoA<typename Particle::SoAArraysType> _soa;

  /// specifies how the verlet lists are build
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas

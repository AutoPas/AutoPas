/**
 * @file VerletListsCells.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "autopas/containers/LinkedCells.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/VerletListsCellsHelpers.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * Linked Cells with Verlet Lists container.
 * The VerletListsCells class uses neighborhood lists for each cell
 * to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @tparam Particle
 */
template <class Particle>
class VerletListsCells : public ParticleContainer<Particle, FullParticleCell<Particle>> {
  typedef VerletListsCellsHelpers<Particle> verlet_internal;
  typedef typename verlet_internal::VerletListParticleCellType ParticleCell;

 private:
  const std::vector<TraversalOptions>& VLApplicableTraversals() {
    switch (_buildTraversal) {
      case c08: {
        static const std::vector<TraversalOptions> v{TraversalOptions::sliced, TraversalOptions::c01};
        return v;
      }
      case c18: {
        static const std::vector<TraversalOptions> v{TraversalOptions::c18, TraversalOptions::c01};
        return v;
      }
      case c01: {
        static const std::vector<TraversalOptions> v{TraversalOptions::c01};
        return v;
      }
      default: {
        static const std::vector<TraversalOptions> v{};
        return v;
      }
    }
  }

 public:
  /**
   * Constructor of the VerletListsCells class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * The rebuildFrequency should be chosen, s.t. the particles do not move more
   * than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals
   * @param buildTraversal the traversal used to build the verletlists
   */
  VerletListsCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff,
                   TraversalOptions buildTraversal, double skin = 0, unsigned int rebuildFrequency = 1)
      : ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff + skin),
        _linkedCells(boxMin, boxMax, cutoff + skin),
        _buildTraversal(buildTraversal),
        _skin(skin),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _neighborListIsValid(false),
        _verletBuiltNewton3(false) {}

  ContainerOptions getContainerType() override { return ContainerOptions::verletListsCells; }

  /**
   * Function to iterate over all pairs of particles.
   * This function only handles short-range interactions.
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal the traversal that will be used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseAoS(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    if (needsRebuild(useNewton3)) {  // if rebuild needed
      this->updateVerletLists(useNewton3);
    }
    traversal->traverseCellVerlet(_neighborLists);

    // we iterated, so increase traversal counter
    _traversalsSinceLastRebuild++;
  }

  /**
   * Dummy function. (Uses AoS instead)
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal the traversal that will be used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseSoA(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    iteratePairwiseAoS(f, traversal, useNewton3);
  }

  /**
   * @copydoc VerletLists::addParticle()
   */
  void addParticle(Particle& p) override {
    _neighborListIsValid = false;
    _linkedCells.addParticle(p);
  }

  /**
   * @copydoc VerletLists::addHaloParticle()
   */
  void addHaloParticle(Particle& haloParticle) override {
    _neighborListIsValid = false;
    _linkedCells.addHaloParticle(haloParticle);
  }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override {
    _neighborListIsValid = false;
    _linkedCells.deleteHaloParticles();
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  void updateContainer() override {
    AutoPasLog(debug, "updating container");
    _neighborListIsValid = false;
    _linkedCells.updateContainer();
  }

  bool isContainerUpdateNeeded() override {
    std::atomic<bool> outlierFound(false);
#ifdef AUTOPAS_OPENMP
    // TODO: find a sensible value for ???
#pragma omp parallel for shared(outlierFound)  // if (this->_cells.size() / omp_get_max_threads() > ???)
#endif
    for (size_t cellIndex1d = 0; cellIndex1d < _linkedCells.getCells().size(); ++cellIndex1d) {
      std::array<double, 3> boxmin{0., 0., 0.};
      std::array<double, 3> boxmax{0., 0., 0.};
      _linkedCells.getCellBlock().getCellBoundingBox(cellIndex1d, boxmin, boxmax);
      boxmin = ArrayMath::addScalar(boxmin, -_skin / 2.);
      boxmax = ArrayMath::addScalar(boxmax, +_skin / 2.);
      for (auto iter = _linkedCells.getCells()[cellIndex1d].begin(); iter.isValid(); ++iter) {
        if (not iter->inBox(boxmin, boxmax)) {
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
    // at the moment this is just a dummy
    return TraversalSelector<ParticleCell>({0, 0, 0}, traversalOptions);
  }

  /**
   * specifies whether the neighbor lists need to be rebuild
   * @param useNewton3 if newton3 is gonna be used to traverse
   * @return true if the neighbor lists need to be rebuild, false otherwise
   */
  bool needsRebuild(bool useNewton3) {
    AutoPasLog(debug, "VerletLists: neighborlist is valid: {}", _neighborListIsValid);
    return (not _neighborListIsValid)                             // if the neighborlist is NOT valid a
                                                                  // rebuild is needed
           or (_traversalsSinceLastRebuild >= _rebuildFrequency)  // rebuild with frequency
           or (useNewton3 != _verletBuiltNewton3);                // useNewton3 changed
  }

  /**
   * @copydoc VerletLists::updateHaloParticle()
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

  ParticleIteratorWrapper<Particle> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return _linkedCells.getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * get the dimension of the used cellblock including the haloboxes
   * @return the dimensions of the used cellblock
   */
  const std::array<std::size_t, 3>& getCellsPerDimension() {
    return _linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
  }

  /**
   * get the neighbors list of a particle
   * @param particle
   * @param useNewton3
   * @return the neighbor list of the particle
   */
  std::vector<Particle*>& getVerletList(Particle* particle, bool useNewton3 = true) {
    if (needsRebuild(useNewton3)) {  // if rebuild needed
      this->updateVerletLists(useNewton3);
    }
    auto indices = _cellMap[particle];
    return _neighborLists[indices.first][indices.second].second;
  }

 protected:
  /**
   * @copydoc VerletLists::checkParticleInCellAndUpdate()
   */
  bool checkParticleInCellAndUpdate(typename verlet_internal::VerletListParticleCellType& cellI, Particle& particleI) {
    for (auto iterator = cellI.begin(); iterator.isValid(); ++iterator) {
      if (iterator->getID() == particleI.getID()) {
        *iterator = particleI;
        return true;
      }
    }
    return false;
  }

  /**
   * update the verlet lists
   */
  void updateVerletLists(bool useNewton3) {
    // the neighbor list is now valid
    _neighborListIsValid = true;
    _traversalsSinceLastRebuild = 0;
    _verletBuiltNewton3 = useNewton3;

    // create a Verlet Lists for each cell
    _neighborLists.clear();
    auto& cells = _linkedCells.getCells();
    size_t cellsSize = cells.size();
    _neighborLists.resize(cellsSize);
    for (size_t cellIndex = 0; cellIndex < cellsSize; ++cellIndex) {
      size_t i = 0;
      for (auto iter = cells[cellIndex].begin(); iter.isValid(); ++iter, ++i) {
        Particle* particle = &*iter;
        _neighborLists[cellIndex].push_back(
            std::pair<Particle*, std::vector<Particle*>>(particle, std::vector<Particle*>()));
        _cellMap[particle] = std::pair<size_t, size_t>(cellIndex, i);
      }
    }

    typename verlet_internal::VerletListGeneratorFunctor f(_neighborLists, _cellMap, this->getCutoff());

    switch (_buildTraversal) {
      case c08: {
        if (useNewton3) {
          auto traversal = C08Traversal<typename verlet_internal::VerletListParticleCellType,
                                        typename verlet_internal::VerletListGeneratorFunctor, false, true>(
              _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          _linkedCells.iteratePairwiseAoS(&f, &traversal);
        } else {
          auto traversal = C08Traversal<typename verlet_internal::VerletListParticleCellType,
                                        typename verlet_internal::VerletListGeneratorFunctor, false, false>(
              _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          _linkedCells.iteratePairwiseAoS(&f, &traversal);
        }
        break;
      }
      case c18: {
        if (useNewton3) {
          auto traversal = C18Traversal<typename verlet_internal::VerletListParticleCellType,
                                        typename verlet_internal::VerletListGeneratorFunctor, false, true>(
              _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          _linkedCells.iteratePairwiseAoS(&f, &traversal);
        } else {
          auto traversal = C18Traversal<typename verlet_internal::VerletListParticleCellType,
                                        typename verlet_internal::VerletListGeneratorFunctor, false, false>(
              _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          _linkedCells.iteratePairwiseAoS(&f, &traversal);
        }
        break;
      }
      case c01: {
        if (not useNewton3) {
          auto traversal = C01Traversal<typename verlet_internal::VerletListParticleCellType,
                                        typename verlet_internal::VerletListGeneratorFunctor, false>(
              _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
          _linkedCells.iteratePairwiseAoS(&f, &traversal);
        } else {
          utils::ExceptionHandler::exception("VerletListsCells::updateVerletLists(): c01 does not support newton3");
        }
        break;
      }
      default:
        utils::ExceptionHandler::exception("VerletListsCells::updateVerletLists(): unsupported Traversal: {}",
                                           _buildTraversal);
        break;
    }
  }

 private:
  /// verlet lists for each particle for each cell
  typename verlet_internal::VerletList_storage_type _neighborLists;

  /// internal linked cells storage, handles Particle storage and used to build verlet lists
  LinkedCells<Particle, typename verlet_internal::VerletListParticleCellType> _linkedCells;

  /// mapping each particle to its corresponding cell and position in this cell
  std::unordered_map<Particle*, std::pair<size_t, size_t>> _cellMap;

  // the traversal used to build the verletlists
  TraversalOptions _buildTraversal;

  /// skin radius
  double _skin;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  // specifies if the neighbor list is currently valid
  bool _neighborListIsValid;

  // specifies if the current verlet list was built for newton3
  bool _verletBuiltNewton3;
};

}  // namespace autopas

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
 * @tparam useNewton3
 */
template <class Particle, bool useNewton3 = true>
class VerletListsCells : public ParticleContainer<Particle, FullParticleCell<Particle>> {
  typedef VerletListsCellsHelpers<Particle> verlet_internal;
  typedef typename verlet_internal::VerletListParticleCellType ParticleCell;

 private:
  static const std::vector<TraversalOptions>& VLApplicableTraversals() {
    static const std::vector<TraversalOptions> v{TraversalOptions::c18};
    return v;
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
   */
  VerletListsCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff,
                   double skin = 0, unsigned int rebuildFrequency = 1)
      : ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff + skin),
        _linkedCells(boxMin, boxMax, cutoff + skin),
        _skin(skin),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _neighborListIsValid(false) {}

  ContainerOptions getContainerType() override { return ContainerOptions::verletListsCells; }

  /**
   * Function to iterate over all pairs of particles.
   * This function only handles short-range interactions.
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal the traversal that will be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor* f, Traversal* traversal) {
    if (needsRebuild()) {  // if rebuild needed
      this->updateVerletLists();
    }
    this->iterateVerletLists(f);
    // we iterated, so increase traversal counter
    _traversalsSinceLastRebuild++;
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
    AutoPasLogger->debug("updating container");
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
    if (outlierFound)
      AutoPasLogger->debug(
          "VerletLists: containerUpdate needed! Particles are fast. You "
          "might want to increase the skin radius or decrease the rebuild "
          "frequency.");
    else
      AutoPasLogger->debug(
          "VerletLists: containerUpdate not yet needed. Particles are slow "
          "enough.");
    return outlierFound;
  }

  TraversalSelector<ParticleCell> generateTraversalSelector(std::vector<TraversalOptions> traversalOptions) override {
    // at the moment this is just a dummy
    return TraversalSelector<ParticleCell>({0, 0, 0}, traversalOptions);
  }

  /**
   * @copydoc VerletLists::needsRebuild()
   */
  bool needsRebuild() {
    AutoPasLogger->debug("VerletLists: neighborlist is valid: {}", _neighborListIsValid);
    return (not _neighborListIsValid)                              // if the neighborlist is NOT valid a
                                                                   // rebuild is needed
           or (_traversalsSinceLastRebuild >= _rebuildFrequency);  // rebuild with frequency
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
      AutoPasLogger->error(
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
   * get the neighbors list of a particle
   * @param particle
   * @return the neighbor list of the particle
   */
  std::vector<Particle*>& getVerletList(Particle* particle) {
    if (needsRebuild()) {  // if rebuild needed
      this->updateVerletLists();
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
  void updateVerletLists() {
    // the neighbor list is now valid
    _neighborListIsValid = true;
    _traversalsSinceLastRebuild = 0;

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

    auto traversal = C18Traversal<typename verlet_internal::VerletListParticleCellType,
                                  typename verlet_internal::VerletListGeneratorFunctor, false, true>(
        _linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f);
    _linkedCells.iteratePairwiseAoS(&f, &traversal);
  }

  /**
   * iterate over the verlet list of a given cell
   * @tparam ParticleFunctor
   * @param f
   * @param cellIndex
   */
  template <class ParticleFunctor>
  inline void iterateVerletListsCell(ParticleFunctor* f, unsigned long cellIndex) {
    for (auto& list : _neighborLists[cellIndex]) {
      Particle& i = *list.first;
      for (auto j_ptr : list.second) {
        Particle& j = *j_ptr;
        f->AoSFunctor(i, j, useNewton3);
        if (not useNewton3) {
          f->AoSFunctor(j, i, false);
        }
      }
    }
  }

  /**
   * iterate over the verlet list
   * @tparam ParticleFunctor
   * @param f
   */
  template <class ParticleFunctor>
  void iterateVerletLists(ParticleFunctor* f) {
    using std::array;

    const auto cellsPerDimension = _linkedCells.getCellBlock().getCellsPerDimensionWithHalo();
    const array<unsigned long, 3> stride = {3, 3, 2};
    array<unsigned long, 3> end;
    end[0] = cellsPerDimension[0];
    end[1] = cellsPerDimension[1];
    end[2] = cellsPerDimension[2] - 1;

#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
    {
      for (unsigned long col = 0; col < 18; ++col) {
        std::array<unsigned long, 3> start = ThreeDimensionalMapping::oneToThreeD(col, stride);

        // intel compiler demands following:
        const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
        const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];
        const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(3)
#endif
        for (unsigned long z = start_z; z < end_z; z += stride_z) {
          for (unsigned long y = start_y; y < end_y; y += stride_y) {
            for (unsigned long x = start_x; x < end_x; x += stride_x) {
              unsigned long cellIndex = ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDimension);
              iterateVerletListsCell(f, cellIndex);
            }
          }
        }
      }
    }
  }

 private:
  /// verlet lists for each particle for each cell
  typename verlet_internal::VerletList_storage_type _neighborLists;

  /// internal linked cells storage, handles Particle storage and used to build verlet lists
  LinkedCells<Particle, typename verlet_internal::VerletListParticleCellType> _linkedCells;

  /// mapping each particle to its corresponding cell and position in this cell
  std::unordered_map<Particle*, std::pair<size_t, size_t>> _cellMap;

  /// skin radius
  double _skin;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  // specifies if the neighbor list is currently valid
  bool _neighborListIsValid;
};

}  // namespace autopas

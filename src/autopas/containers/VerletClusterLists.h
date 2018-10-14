/**
 * @file VerletClusterLists.h
 * @author nguyen
 * @date 14.10.18
 */

#pragma once

#include <cmath>
#include "autopas/containers/ParticleContainer.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * Particles are divided into clusters
 * The VerletClusterLists class uses neighborhood lists for each cluster
 * to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * @tparam Particle
 */
template <class Particle>
class VerletClusterLists : public ParticleContainer<Particle, std::vector<Particle>> {
  /**
   * the index type to access the particle cells
   */
  typedef std::size_t index_t;

 public:
  /**
   * Constructor of the VerletClusterLists class.
   * The neighbor lists are build using a estimated density.
   * The box is divided into cuboids with roughly with approximately the
   * same side length. The rebuildFrequency should be chosen, s.t. the particles do
   * not move more than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param density the estimated particles per volume
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals
   * @param clusterSize size of clusters
   */
  VerletClusterLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double density,
                     double cutoff, double skin = 0, unsigned int rebuildFrequency = 1, int clusterSize = 4)
      : ParticleContainer<Particle, std::vector<Particle>>(boxMin, boxMax, cutoff + skin),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _skin(skin),
        _cutoff(cutoff),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _neighborListIsValid(false),
        _clusterSize(clusterSize) {
    _gridSideLength = cbrt(clusterSize / density);
    _gridSideLengthReciprocal = 1 / _gridSideLength;

    index_t numCells = 1;
    for (int d = 0; d < 2; d++) {
      double diff = _boxMax[d] - _boxMin[d];
      _cellsPerDim[d] = static_cast<index_t>(std::floor(diff / _gridSideLength));
      // at least one central cell
      _cellsPerDim[d] = std::max(_cellsPerDim[d], 1ul);

      numCells *= _cellsPerDim[d];
    }
    _cellsPerDim[2] = 1ul;

    _clusters.resize(numCells);
  }

  ContainerOptions getContainerType() override { return ContainerOptions::verletListsCells; }

  /**
   * Function to iterate over all pairs of particles.
   * This function only handles short-range interactions.
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal not used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseAoS(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = true) {
    if (needsRebuild()) {  // if rebuild needed
      this->updateVerletLists();
    }
    traverseVerletLists(f, useNewton3);

    // we iterated, so increase traversal counter
    _traversalsSinceLastRebuild++;
  }

  /**
   * Dummy function. (Uses AoS instead)
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal not used
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
    // add particle somewhere, because lists will be rebuild anyways
    _clusters[0].push_back(&p);
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  void updateContainer() override {
    AutoPasLogger->debug("updating container");
    _neighborListIsValid = false;
  }

  /**
   * specifies whether the neighbor lists need to be rebuild
   * @return true if the neighbor lists need to be rebuild, false otherwise
   */
  bool needsRebuild() {
    AutoPasLogger->debug("VerletLists: neighborlist is valid: {}", _neighborListIsValid);
    return (not _neighborListIsValid)                              // if the neighborlist is NOT valid a
                                                                   // rebuild is needed
           or (_traversalsSinceLastRebuild >= _rebuildFrequency);  // rebuild with frequency
  }

 protected:
  /*
   * recalculate clusters
   */
  void rebuild() {
    std::vector<Particle> invalidParticles;

    for (auto& cluster : _clusters) {
      for (auto& particle : cluster) {
        invalidParticles.push_back(particle);
      }

      cluster.clear();
    }

    // put particles into grids
    for (auto& particle : invalidParticles) {
      if (inBox(particle.getR(), _boxMin, _boxMax)) {
        auto index = getIndexOfPosition(particle.getR());
        _clusters[index].push_back(particle);
      }
    }

    // sort by last dimension
    for (auto& cluster : _clusters) {
      std::sort(cluster.begin(), cluster.end(),
                [](const Particle& a, const Particle& b) -> bool { return a.GetR()[2] < b.GetR()[2]; });
    }
  }

  /**
   * update the verlet lists
   */
  void updateVerletLists() {
    // the neighbor list is now valid
    _neighborListIsValid = true;
    _traversalsSinceLastRebuild = 0;

    rebuild();

    const int boxRange = static_cast<int>(ceil((_cutoff + _skin) / _gridSideLength));

    const int gridMaxX = _cellsPerDim[0] - 1;
    const int gridMaxY = _cellsPerDim[1] - 1;
    for (int yi = 0; yi <= gridMaxY; yi++) {
      for (int xi = 0; xi <= gridMaxX; xi++) {
        auto ci = _clusters[index1D(xi, yi)];
        const int minX = std::max(xi - boxRange, 0);
        const int minY = std::max(yi - boxRange, 0);
        const int maxX = std::min(xi + boxRange, gridMaxX);
        const int maxY = std::min(yi + boxRange, gridMaxY);

        for (int yj = minY; yj <= maxY; yj++) {
          for (int xj = minX; xj <= maxX; xj++) {
            auto cj = _clusters[index1D(xj, yj)];

            // TODO
          }
        }
      }
    }
  }

  template <class ParticleFunctor>
  void traverseVerletLists(ParticleFunctor* functor, bool useNewton3) {
    // TODO
  }

  inline index_t getIndexOfPosition(const std::array<double, 3>& pos) const {
    std::array<index_t, 2> cellIndex{};

    for (int dim = 0; dim < 2; dim++) {
      const long int value = (static_cast<long int>(floor((pos[dim] - _boxMin[dim]) * _gridSideLengthReciprocal))) + 1l;
      const index_t nonnegativeValue = static_cast<index_t>(std::max(value, 0l));
      const index_t nonLargerValue = std::min(nonnegativeValue, _cellsPerDim[dim] - 1);
      cellIndex[dim] = nonLargerValue;
      /// @todo this is a sanity check to prevent doubling of particles, but
      /// could be done better!
      if (pos[dim] >= _boxMax[dim]) {
        cellIndex[dim] = _cellsPerDim[dim] - 1;
      } else if (pos[dim] < _boxMin[dim]) {
        cellIndex[dim] = 0;
      }
    }

    return cellIndex[0] + cellIndex[1] * _cellsPerDim[0];
    // in very rare cases rounding is stupid, thus we need a check...
    /// @todo when the border and flag manager is there
  }

  inline index_t index1D(const int x, const int y) { return x + y * _cellsPerDim[0]; }

 private:
  // TODO
  // std::vector<std::vector<Particle*>> _neighborLists;

  /// internal storage, particles are split into a grid in xy-dimension
  std::vector<std::vector<Particle>> _clusters;
  int _clusterSize;

  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;

  // side length of xy-grid and reciprocal
  double _gridSideLength;
  double _gridSideLengthReciprocal;

  // dimensions of grid
  std::array<index_t, 3> _cellsPerDim;

  /// skin radius
  double _skin;

  /// cutoff
  double _cutoff;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  // specifies if the neighbor list is currently valid
  bool _neighborListIsValid;
};

}  // namespace autopas

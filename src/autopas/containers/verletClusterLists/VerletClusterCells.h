/**
 * @file VerletClusterCluster.h
 * @author jspahl
 * @date 25.3.19
 */

#pragma once

#include <algorithm>
#include <cmath>
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/cellPairTraversals/VerletClusterTraversalInterface.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

/**
 * Particles are divided into clusters.
 * The VerletClusterLists class uses neighborhood lists for each cluster pair
 * to calculate pairwise interactions.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * @tparam Particle
 */
template <class Particle>
class VerletClusterCells : public ParticleContainer<Particle, FullParticleCell<Particle>> {
  /**
   * the index type to access the particle cells
   */
  typedef std::size_t index_t;

 public:
  /**
   * Constructor of the VerletClusterLists class.
   * The neighbor lists are build using a estimated density.
   * The box is divided into cuboids with roughly the
   * same side length. The rebuildFrequency should be chosen, s.t. the particles do
   * not move more than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals.
   * @param clusterSize size of clusters
   */
  VerletClusterCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff,
                     double skin = 0, unsigned int rebuildFrequency = 1, int clusterSize = 32)
      : ParticleContainer<Particle, FullParticleCell<Particle>>(boxMin, boxMax, cutoff + skin,
                                                                allVCLApplicableTraversals()),
        _firstHaloClusterId(1),
        _clusterSize(clusterSize),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _skin(skin),
        _cutoff(cutoff),
        _cutoffSqr(cutoff * cutoff),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _needsRebuild(false) {
    _clusters.resize(1);
    _dummyStarts.resize(1, 0);
  }

  /**
   * Lists all traversal options applicable for the Verlet Lists container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOption>& allVCLApplicableTraversals() {
    // traversal not used but prevents usage of newton3
    static const std::vector<TraversalOption> v{TraversalOption::verletClusterCellsTraversal};
    return v;
  }

  std::vector<TraversalOption> getAllTraversals() override { return allVCLApplicableTraversals(); }

  ContainerOption getContainerType() override { return ContainerOption::verletClusterCells; }

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
  void iteratePairwise(ParticleFunctor* f, Traversal* traversal, bool useNewton3 = false) {
    auto* traversalInterface = dynamic_cast<VerletClusterTraversalInterface<FullParticleCell<Particle>>*>(traversal);
    if (!traversalInterface) {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in VerletClusterCells::iteratePairwise");
    }
    if (needsRebuild()) {
      this->rebuild();
      traversalInterface->rebuild(_cellsPerDim, _clusterSize, _clusters, _boundingBoxes);
    }

    traversal->initTraversal(_clusters);
    traversalInterface->traverseCellPairs(_clusters);
    traversal->endTraversal(_clusters);

    // we iterated, so increase traversal counter
    _traversalsSinceLastRebuild++;
  }

  /**
   * @copydoc VerletLists::addParticle()
   */
  void addParticle(Particle& p) override {
    if (autopas::utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
      _needsRebuild = false;
      _clusters[0].resize(_dummyStarts[0]);
      // add particle somewhere, because lists will be rebuild anyways
      _clusters[0].addParticle(p);
      ++_dummyStarts[0];
    } else {
      utils::ExceptionHandler::exception(
          "VerletCluster: trying to add particle that is not inside the bounding box.\n" + p.toString());
    }
  }

  /**
   * @copydoc VerletLists::addHaloParticle()
   */
  void addHaloParticle(Particle& haloParticle) override {
    _needsRebuild = false;
    _haloInsertQueue.push_back(haloParticle);
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.addHaloParticle not yet implemented.");
  }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override {
    _clusters.resize(_firstHaloClusterId);
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.deleteHaloParticles not yet implemented.");
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  void updateContainer() override {
    AutoPasLog(debug, "updating container");
    _needsRebuild = false;
  }

  bool isContainerUpdateNeeded() override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.isContainerUpdateNeeded not yet implemented");
    return false;
  }

  TraversalSelector<FullParticleCell<Particle>> generateTraversalSelector() override {
    return TraversalSelector<FullParticleCell<Particle>>(_cellsPerDim);
  }

  /**
   * Specifies whether the neighbor lists need to be rebuild.
   * @return true if the neighbor lists need to be rebuild, false otherwise
   */
  bool needsRebuild() {
    AutoPasLog(debug, "VerletLists: neighborlist is valid: {}", _needsRebuild);
    // if the neighbor list is NOT valid or we have not rebuild for _rebuildFrequency steps
    return (not _needsRebuild) or (_traversalsSinceLastRebuild >= _rebuildFrequency);
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, FullParticleCell<Particle>>(&this->_clusters));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(std::array<double, 3> lowerCorner,
                                                      std::array<double, 3> higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    // @todo implement this if bounding boxes are here
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.getRegionIterator not yet implemented.");
    return ParticleIteratorWrapper<Particle>();
  }

  /**
   * Deletes all Dummy Particles in the container
   */
  void deleteDummyParticles() {
    for (size_t i = 0; i < _dummyStarts.size(); ++i) {
      _clusters[i].resize(_dummyStarts[i]);
    }
    _needsRebuild = false;
  }

 protected:
  /**
   * Recalculate grids and clusters,
   * build verlet lists and
   * pad clusters.
   * @param useNewton3
   */
  void rebuild() {
    deleteDummyParticles();
    // get the dimensions and volumes of the box
    std::array<double, 3> boxSize{};
    double volume = 1.0;

    for (int d = 0; d < 3; ++d) {
      boxSize[d] = _boxMax[d] - _boxMin[d];
      volume *= boxSize[d];
    }

    // get all particles and clear clusters
    std::vector<Particle> invalidParticles;
    for (auto& cluster : _clusters) {
      for (index_t i = 0; i < cluster._particles.size(); ++i) {
        if (utils::inBox(cluster._particles[i].getR(), _boxMin, _boxMax))
          invalidParticles.push_back(cluster._particles[i]);
      }
      cluster.clear();
    }

    // estimate particle density
    double density = (double)invalidParticles.size() / volume;

    // guess optimal grid side length
    _gridSideLength = std::cbrt(_clusterSize / density);
    _gridSideLengthReciprocal = 1 / _gridSideLength;

    // get cells per dimension
    index_t sizeGrid = 1;
    for (int d = 0; d < 2; d++) {
      _cellsPerDim[d] = static_cast<index_t>(std::ceil(boxSize[d] * _gridSideLengthReciprocal));

      sizeGrid *= _cellsPerDim[d];
    }
    _cellsPerDim[2] = static_cast<index_t>(1);

    // resize to number of grids
    _clusters.resize(sizeGrid);
    _dummyStarts.resize(sizeGrid);

    // put particles into grid cells
    for (auto& particle : invalidParticles) {
      int index = (int)((particle.getR()[0] - _boxMin[0]) * _gridSideLengthReciprocal) +
                  (int)((particle.getR()[1] - _boxMin[1]) * _gridSideLengthReciprocal) * _cellsPerDim[0];
      _clusters[index].addParticle(particle);
    }

    index_t numClusters = 0;
    // sort by last dimension
    for (auto& cluster : _clusters) {
      cluster.sortByDim(2);
      numClusters += ((cluster.numParticles()) / _clusterSize) + 1;
    }
    numClusters = std::max((index_t)1, numClusters);

    _clusters.resize(numClusters);
    _boundingBoxes.resize(numClusters, {_boxMax[0] + 8 * _cutoff, _boxMax[1] + 8 * _cutoff, _boxMax[2] + 8 * _cutoff,
                                        _boxMin[0] - 8 * _cutoff, _boxMin[1] - 8 * _cutoff, _boxMin[2] - 8 * _cutoff});

    index_t cid = sizeGrid;
    for (index_t i = 0; i < sizeGrid; ++i) {
      index_t dummyStart = _clusters[i].numParticles() % _clusterSize;
      // no cell with only dummy particles except when grid already empty
      if (dummyStart == 0 && !(_clusters[i].numParticles() == 0)) {
        dummyStart = _clusterSize;
      }
      for (index_t clusterStart = dummyStart; clusterStart < _clusters[i].numParticles();
           clusterStart += _clusterSize) {
        for (index_t pid = 0; pid < _clusterSize; ++pid) {
          Particle& p = _clusters[i]._particles[clusterStart + pid];
          expandBoundingBox(_boundingBoxes[cid], p);
          _clusters[cid].addParticle(p);
        }
        ++cid;
      }
      _clusters[i].resize(dummyStart);
      for (auto& p : _clusters[i]._particles) {
        expandBoundingBox(_boundingBoxes[i], p);
      }
      _dummyStarts[i] = dummyStart;
      for (index_t pid = dummyStart; pid < _clusterSize; ++pid) {
        Particle dummyParticle =
            Particle({_boxMax[0] + 8 * _cutoff + static_cast<typename Particle::ParticleFloatingPointType>(i),
                      _boxMax[1] + 8 * _cutoff + static_cast<typename Particle::ParticleFloatingPointType>(pid),
                      _boxMax[2] + 8 * _cutoff},
                     {0., 0., 0.}, ULONG_MAX);
        _clusters[i].addParticle(dummyParticle);
      }
    }

    _traversalsSinceLastRebuild = 0;
  }

 private:
  /**
   * Expands a bounding Box such the Particle is in it
   * @param box
   * @param p
   */
  void expandBoundingBox(std::array<typename Particle::ParticleFloatingPointType, 6>& box, Particle& p) {
    for (int i = 0; i < 3; ++i) {
      box[i] = std::min(box[i], p.getR()[i] - _skin - _cutoff);
      box[3 + i] = std::max(box[3 + i], p.getR()[i] + _skin + _cutoff);
    }
  }

  /// internal storage, particles are split into a grid in xy-dimension and then cut in z dimension
  std::vector<FullParticleCell<Particle>> _clusters;

  /// first Cluster ID to be a Halo Cluster
  index_t _firstHaloClusterId;

  /// indeces where dummy particles in the cells start
  std::vector<index_t> _dummyStarts;

  /// Stores Halo Particles to be inserted
  std::vector<Particle> _haloInsertQueue;

  // number of particles in a cluster
  unsigned int _clusterSize;

  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;

  // bounding boxes of all clusters (xmin,ymin,zmin,xmax,ymax,zmax) including skin and cutoff
  std::vector<std::array<typename Particle::ParticleFloatingPointType, 6>> _boundingBoxes;

  // side length of xy-grid and reciprocal
  double _gridSideLength;
  double _gridSideLengthReciprocal;

  // dimensions of grid
  std::array<index_t, 3> _cellsPerDim;

  /// skin radius
  double _skin;

  /// cutoff
  double _cutoff;
  double _cutoffSqr;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  // specifies if the neighbor list is currently valid
  bool _needsRebuild;
};

}  // namespace autopas

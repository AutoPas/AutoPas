/**
 * @file VerletClusterLists.h
 * @author nguyen
 * @date 14.10.18
 */

#pragma once

#include <cmath>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/verletClusterLists/ClusterTower.h"
#include "autopas/containers/verletClusterLists/VerletClusterListsRebuilder.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Timer.h"

namespace autopas {

/**
 * Particles are divided into clusters.
 * The VerletClusterLists class uses neighborhood lists for each cluster
 * to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * @tparam Particle
 */
template <class Particle>
class VerletClusterLists : public ParticleContainer<FullParticleCell<Particle>> {
 public:
  /**
   * The number of particles in a full cluster. Currently, constexpr is necessary so it can be passed to ClusterTower as
   * a template parameter.
   */
  static constexpr size_t clusterSize = 4;

  /**
   * Constructor of the VerletClusterLists class.
   * The neighbor lists are build using a estimated density.
   * The box is divided into cuboids with roughly the
   * same side length.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   */
  VerletClusterLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff,
                     double skin = 0)
      : ParticleContainer<FullParticleCell<Particle>>(boxMin, boxMax, cutoff, skin),
        _numClusters(0),
        _numTowersPerInteractionLength(0),
        _neighborListIsNewton3(false) {}

  ContainerOption getContainerType() const override { return ContainerOption::verletClusterLists; }

  void iteratePairwise(TraversalInterface *traversal) override {
    AutoPasLog(debug, "Using traversal {}.", traversal->getTraversalType().to_string());

    auto *traversalInterface = dynamic_cast<VerletClustersTraversalInterface<Particle> *>(traversal);
    if (traversalInterface) {
      traversalInterface->setClusterLists(*this);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in VerletClusterLists::iteratePairwise. TraversalID: {}",
          traversal->getTraversalType());
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  /// @todo: Somehow make the iterator also iterate over the _particlesToAdd. Otherwise, e.g. ContainerSelectorTest
  /// will work but have all particles removed after it switches from this container to another.
  /// see: https://github.com/AutoPas/AutoPas/issues/155
  /**
   * Adds the given particle to the container. rebuildVerletLists() has to be called to have it actually sorted in.
   * @param p The particle to add.
   */
  void addParticle(const Particle &p) override { _particlesToAdd.push_back(p); }

  /**
   * @copydoc VerletLists::addHaloParticle()
   */
  void addHaloParticle(const Particle &haloParticle) override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.addHaloParticle not yet implemented.");
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const Particle &haloParticle) override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.updateHaloParticle not yet implemented.");
    return false;
  }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override {
    // quick and dirty: iterate over all particles and delete halo particles
    /// @todo: make this proper
    for (auto iter = this->begin(IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
      if (not iter->isOwned()) {
        internal::deleteParticle(iter);
      }
    }
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  AUTOPAS_WARN_UNUSED_RESULT
  std::vector<Particle> updateContainer() override {
    /// @todo What happens when some particles are just deleted here?
    AutoPasLog(debug, "updating container");
    // first delete all particles
    this->deleteHaloParticles();

    // next find invalid particles
    std::vector<Particle> invalidParticles;
    /// @todo: parallelize
    for (auto iter = this->begin(IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      if (not utils::inBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
        invalidParticles.push_back(*iter);
        internal::deleteParticle(iter);
      }
    }

    return invalidParticles;
  }

  bool isContainerUpdateNeeded() const override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.isContainerUpdateNeeded not yet implemented");
    return false;
  }

  TraversalSelectorInfo getTraversalSelectorInfo() const override {
    std::array<double, 3> towerSize = {_towerSideLength, _towerSideLength, this->getBoxMax()[2] - this->getBoxMin()[2]};
    std::array<unsigned long, 3> towerDimensions = {_towersPerDim[0], _towersPerDim[1], 1};
    return TraversalSelectorInfo(towerDimensions, this->getInteractionLength(), towerSize, clusterSize);
  }

  ParticleIteratorWrapper<Particle, true> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle, true>(
        new internal::ParticleIterator<Particle, internal::ClusterTower<Particle, clusterSize>, true>(
            &(this->_towers)));
  }

  ParticleIteratorWrapper<Particle, false> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    return ParticleIteratorWrapper<Particle, false>(
        new internal::ParticleIterator<Particle, internal::ClusterTower<Particle, clusterSize>, false>(
            &(this->_towers)));
  }

  ParticleIteratorWrapper<Particle, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    /// @todo implement this if bounding boxes are here
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.getRegionIterator not yet implemented.");
    return ParticleIteratorWrapper<Particle, true>();
  }

  ParticleIteratorWrapper<Particle, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    /// @todo implement this if bounding boxes are here
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.getRegionIterator not yet implemented.");
    return ParticleIteratorWrapper<Particle, false>();
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    internal::VerletClusterListsRebuilder<Particle> builder{*this, _particlesToAdd, traversal->getUseNewton3()};
    auto res = builder.rebuild();
    std::tie(_towerSideLength, _numTowersPerInteractionLength, _towersPerDim, _numClusters, _neighborListIsNewton3) = {
        res._towerSideLength, res._interactionLengthInTowers, res._towersPerDim, res._numClusters,
        res._neighborListIsNewton3};
  }

  /**
   * Helper method to iterate over all clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @tparam inParallel If the iteration should be executed in parallel or sequential.  See traverseClustersParallel()
   * for thread safety.
   * @param loopBody The lambda to execute for all clusters. Parameters given is internal::Cluster& cluster.
   */
  template <bool inParallel, class LoopBody>
  void traverseClusters(LoopBody &&loopBody) {
    if (inParallel) {
      traverseClustersParallel<LoopBody>(std::forward<LoopBody>(loopBody));
    } else {
      traverseClustersSequential<LoopBody>(std::forward<LoopBody>(loopBody));
    }
  }

  unsigned long getNumParticles() const override {
    unsigned long sum = 0;
    for (size_t index = 0; index < _towers.size(); index++) {
      sum += _towers[index].getNumActualParticles();
    }
    return sum;
  }

  /**
   * Returns the number of clusters in this container.
   * @return The number of clusters in this container.
   */
  auto getNumClusters() const { return _numClusters; }

  /**
   * Returns the grid side length of the grids in the container.
   * @return the grid side length of the grids in the container.
   */
  auto getTowerSideLength() const { return _towerSideLength; }

  /**
   * Returns the number of grids per dimension on the container.
   * @return the number of grids per dimension on the container.
   */
  auto getTowersPerDimension() const { return _towersPerDim; }

  /**
   * Returns the 2D grid for the XY-plane of this container that defines the cluster towers.
   * @return the grids of this container for usage in traversals.
   */
  auto &getTowers() { return _towers; }

  /**
   * Returns the number of particles in each cluster.
   * @return the number of particles in each cluster.
   */
  constexpr auto getClusterSize() const { return clusterSize; }

  /**
   * Returns the towers per interaction length. That is how many towers fit into one interaction length rounded up.
   * @return the number of towers per interaction length.
   */
  auto getNumTowersPerInteractionLength() const { return _numTowersPerInteractionLength; }

  /**
   * Loads all particles of the container in their correct SoA and generates the SoAViews for the clusters.
   * @tparam Functor The type of the functor to use.
   * @param functor The functor to use for loading the particles into the SoA.
   */
  template <class Functor>
  void loadParticlesIntoSoAs(Functor *functor) {
    const auto numTowers = _towers.size();
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t index = 0; index < numTowers; index++) {
      _towers[index].loadSoA(functor);
    }
  }

  /**
   * Extracts all SoAs of the container into the particles.
   * @tparam Functor The type of the functor to use.
   * @param functor The functor to use for extracting the SoAs into the particles..
   */
  template <class Functor>
  void extractParticlesFromSoAs(Functor *functor) {
    const auto numTowers = _towers.size();
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t index = 0; index < numTowers; index++) {
      _towers[index].extractSoA(functor);
    }
  }

  /**
   * Returns the tower for the given tower grid coordinates.
   * @param x The x-th tower in x direction.
   * @param y The y-th tower in y direction.
   * @return The tower for the given tower grid coordinates.
   */
  auto &getTowerAtCoordinates(const size_t x, const size_t y) { return _towers[towerIndex2DTo1D(x, y)]; }

  /**
   * Returns the 1D index for the given tower grid coordinates of a tower.
   *
   * @param x The x-coordinate of the tower.
   * @param y The y-coordinate of the tower.
   * @param towersPerDim The number of towers in each dimension.
   * @return the 1D index for the given tower grid coordinates of a tower.
   */
  static auto towerIndex2DTo1D(const size_t x, const size_t y, const std::array<size_t, 2> towersPerDim) {
    return x + y * towersPerDim[0];
  }

  /**
   * Returns the 1D index for the given 2D-coordinates of a tower.
   *
   * @param x The x-coordinate of the tower.
   * @param y The y-coordinate of the tower.
   * @return the 1D index for the given 2D-coordinates of a tower.
   */
  [[nodiscard]] size_t towerIndex2DTo1D(const size_t x, const size_t y) const {
    return towerIndex2DTo1D(x, y, _towersPerDim);
  }

 protected:
  /**
   * Helper method to sequentially iterate over all clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @param loopBody The lambda to execute for all clusters. Parameters given is internal::Cluster& cluster.
   */
  template <class LoopBody>
  void traverseClustersSequential(LoopBody &&loopBody) {
    for (size_t x = 0; x < _towersPerDim[0]; x++) {
      for (size_t y = 0; y < _towersPerDim[1]; y++) {
        auto &tower = getTowerAtCoordinates(x, y);
        for (auto &cluster : tower.getClusters()) {
          loopBody(cluster);
        }
      }
    }
  }

  /**
   * Helper method to iterate over all clusters in parallel.
   *
   * It is always safe to modify the particles in the cluster that is passed to the given loop body. However, when
   * modifying particles from other clusters, the caller has to make sure that no data races occur. Particles must not
   * be added or removed during the traversal.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @param loopBody The lambda to execute for all clusters. Parameters given is internal::Cluster& cluster.
   */
  template <class LoopBody>
  void traverseClustersParallel(LoopBody &&loopBody) {
    const auto towersPerDimX = _towersPerDim[0];
    const auto towersPerDimY = _towersPerDim[1];
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (size_t x = 0; x < towersPerDimX; x++) {
      for (size_t y = 0; y < towersPerDimY; y++) {
        auto &tower = getTowerAtCoordinates(x, y);

        for (auto &cluster : tower.getClusters()) {
          loopBody(cluster);
        }
      }
    }
  }

 private:
  /**
   * internal storage, particles are split into a grid in xy-dimension
   */
  std::vector<internal::ClusterTower<Particle, clusterSize>> _towers{1};

  /**
   * Dimensions of the 2D xy-grid.
   */
  std::array<size_t, 2> _towersPerDim{};

  /**
   * Side length of xy-grid.
   */
  double _towerSideLength{0.};

  /**
   * The number of clusters in the container.
   */
  size_t _numClusters;

  /**
   * The interaction length in number of towers it reaches.
   * static_cast<int>(std::ceil((this->getInteractionLength()) * _towerSideLengthReciprocal))
   */
  int _numTowersPerInteractionLength;

  /**
   * Specifies if the saved neighbors use newton 3 or not.
   */
  bool _neighborListIsNewton3;

  /**
   * Contains all particles that should be added to the container during the next rebuild.
   */
  std::vector<Particle> _particlesToAdd;
};

}  // namespace autopas

/**
 * @file ClusterTower.h
 * @author humig
 * @date 27.07.19
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/verletClusterLists/Cluster.h"

namespace autopas::internal {

template <class Particle, size_t clusterSize>
class ClusterTower : public ParticleCell<Particle> {
  static const Particle dummy;

 public:
  /**
   * Adds a particle to the cluster tower. If generateClusters() has already been called on this ClusterTower, clear()
   * must be called first, or dummies are messed up!
   *
   * Is allowed to be called in parallel since a lock is used on the internal cell.
   * @param particle The particle to add.
   */
  void addParticle(const Particle &particle) override { _particles.addParticle(particle); }

  void clear() override {
    _clusters.clear();
    _particles.clear();
    _numDummyParticles = 0;
  }

  size_t generateClusters() {
    if (getNumActualParticles() > 0) {
      _particles.sortByDim(2);

      auto sizeLastCluster = (_particles.numParticles() % clusterSize);
      _numDummyParticles = sizeLastCluster != 0 ? clusterSize - sizeLastCluster : 0;

      const auto &lastParticle = _particles[_particles.numParticles() - 1];
      for (size_t i = 0; i < _numDummyParticles; i++) {
        _particles.addParticle(lastParticle);
      }

      size_t numClusters = _particles.numParticles() / clusterSize;
      _clusters.reserve(numClusters);
      for (size_t index = 0; index < numClusters; index++) {
        _clusters.emplace_back(&(_particles[clusterSize * index]));
      }
    }

    return getNumClusters();
  }

  void fillUpWithDummyParticles(double dummyStartX, double dummyDistZ) {
    auto &lastCluster = getCluster(getNumClusters() - 1);
    for (size_t index = 1; index <= _numDummyParticles; index++) {
      lastCluster.getParticle(clusterSize - index) = dummy;
      lastCluster.getParticle(clusterSize - index).setR({dummyStartX, 0, dummyDistZ * index});
    }
  }

  template <class Functor>
  void loadSoA(Functor *functor) {
    functor->SoALoader(_particles, _particles._particleSoABuffer);
    for (size_t index = 0; index < getNumClusters(); index++) {
      auto &cluster = getCluster(index);
      cluster.getSoAView() = {&(_particles._particleSoABuffer), index * clusterSize, (index + 1) * clusterSize};
    }
  }

  template <class Functor>
  void extractSoA(Functor *functor) {
    functor->SoAExtractor(_particles, _particles._particleSoABuffer);
  }

  [[nodiscard]] size_t getNumDummyParticles() const { return _numDummyParticles; }

      [[nodiscard]] size_t getNumActualParticles() const {
    return _particles.numParticles() - _numDummyParticles;
  }

  [[nodiscard]] size_t getNumClusters() const { return _clusters.size(); }

      [[nodiscard]] auto &getClusters() {
    return _clusters;
  }

  [[nodiscard]] auto &getCluster(size_t index) { return _clusters[index]; }

  [[nodiscard]] auto &getCluster(size_t index) const { return _clusters[index]; }

  [[nodiscard]] unsigned long numParticles() const override { return getNumActualParticles(); }

      [[nodiscard]] SingleCellIteratorWrapper<Particle> begin() override {
    return SingleCellIteratorWrapper<Particle>{
        new SingleCellIterator<Particle, ClusterTower<Particle, clusterSize>>(this)};
  }

  /**
   * Returns the particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  decltype(auto) getParticle(size_t index) { return _particles._particles.at(index); }

  // Methods from here on: Only to comply with ParticleCell interface. SingleCellIterators work on ParticleCells, and
  // while those methods would not be needed, still complying to the whole interface should be helpful, if
  // maybe someday new necessary pure virtual methods are introduced there.

  [[nodiscard]] bool isNotEmpty() const override { return getNumActualParticles() > 0; }

  void deleteByIndex(size_t index) override {
    // TODO support deletion of particles somehow
    autopas::utils::ExceptionHandler::exception("Not supported!");
  }

  void setCellLength(std::array<double, 3> &) override {
    autopas::utils::ExceptionHandler::exception("Not supported!");
  }

  [[nodiscard]] std::array<double, 3> getCellLength() const override {
    autopas::utils::ExceptionHandler::exception("Not supported!");
    return {0, 0, 0};
  }

  private : std::vector<Cluster<Particle, clusterSize>> _clusters;
  FullParticleCell<Particle> _particles;
  size_t _numDummyParticles{};
};

template <class Particle, size_t clusterSize>
const Particle ClusterTower<Particle, clusterSize>::dummy{
    {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()},
    {0, 0, 0},
    0};

}  // namespace autopas::internal

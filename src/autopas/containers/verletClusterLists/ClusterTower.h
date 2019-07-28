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
class ClusterTower {
  static const Particle dummy;

 public:
  void addParticle(const Particle &particle) { _particles.addParticle(particle); }

  size_t generateClusters() {
    if (getNumParticles() > 0) {
      _particles.sortByDim(2);
      _numDummyParticles = clusterSize - getNumParticles() % clusterSize;
      const auto &lastParticle = _particles[getNumParticles() - 1];
      for (size_t i = 0; i < _numDummyParticles; i++) {
        _particles.addParticle(lastParticle);
      }

      size_t numClusters = getNumParticles() / clusterSize;
      _clusters.reserve(numClusters);
      for (size_t index = 0; index < numClusters; index++) {
        _clusters.emplace_back(&(_particles[clusterSize * index]));
      }
    }

    return getNumClusters();
  }

  void clear() {
    _clusters.clear();
    _particles.clear();
    _numDummyParticles = 0;
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

      [[nodiscard]] size_t getNumParticles() const {
    return _particles.numParticles();
  }

  [[nodiscard]] size_t getNumClusters() const { return _clusters.size(); }

  decltype(auto) begin() {
    return _particles.begin();
  }

  [[nodiscard]] auto &getClusters() { return _clusters; }

  [[nodiscard]] auto &getCluster(size_t index) { return _clusters[index]; }

  [[nodiscard]] bool isNotEmpty() const { return getNumParticles() > 0; }

  private : std::vector<Cluster<Particle, clusterSize>> _clusters;
  FullParticleCell<Particle> _particles;
  size_t _numDummyParticles;
};

template <class Particle, size_t clusterSize>
const Particle ClusterTower<Particle, clusterSize>::dummy{
    {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()},
    {0, 0, 0},
    0};

}  // namespace autopas::internal

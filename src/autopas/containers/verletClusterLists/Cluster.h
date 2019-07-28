/**
 * @file Cluster.h
 * @author humig
 * @date 27.07.19
 */

#pragma once

#include <vector>
#include "autopas/utils/SoAView.h"

namespace autopas::internal {

template <class Particle, size_t clusterSize>
class Cluster {
 public:
  explicit Cluster(Particle *firstParticle) : _firstParticle(firstParticle) {}

  auto &getParticle(size_t index) { return *(_firstParticle + index); }

  const auto &getParticle(size_t index) const { return *(_firstParticle + index); }

  auto &getSoAView() { return _soaView; }

  const auto &getNeighbors() const { return _neighborClusters; }

  void addNeighbor(Cluster<Particle, clusterSize> &neighbor) { _neighborClusters.push_back(&neighbor); }

 private:
  Particle *_firstParticle;
  SoAView<typename Particle::SoAArraysType> _soaView;
  std::vector<Cluster *> _neighborClusters;
};

}  // namespace autopas::internal

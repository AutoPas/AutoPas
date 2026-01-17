#pragma once

#include "autopas/containers/ParticleContainerInterface.h"
#include "traversals/KokkosVCLTraversal.h"

namespace autopas {

template <typename Particle_T>
class KokkosVerletClusterLists : public autopas::ParticleContainerInterface<Particle_T> {
 public:
  KokkosVerletClusterLists(const std::array<double, 3> &box_min, const std::array<double, 3> &box_max,
                           const double cutoff, const double verlet_skin, const int rebuild_frequency)
      : _boxMin(box_min),
        _boxMax(box_max),
        _cutoff(cutoff),
        _verletSkin(verlet_skin),
        _rebuildFrequency(rebuild_frequency) {}

  autopas::CellType getParticleCellTypeEnum() const override { return autopas::CellType::IsNoCell; }

  [[nodiscard]] autopas::ContainerOption getContainerType() const override {
    return autopas::ContainerOption::kokkosVerletClusterLists;
  }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    std::lock_guard<AutoPasLock> lock(_particlesLock);
    particles.reserve(numParticles + numParticlesHaloEstimate);
  }

 protected:
  void addParticleImpl(const Particle_T &p) override {
    std::lock_guard<AutoPasLock> lock(_particlesLock);
    particles.push_back(p);
  }
  void addHaloParticleImpl(const Particle_T &haloParticle) override {
    std::lock_guard<AutoPasLock> lock(_particlesLock);
    particles.push_back(haloParticle);
  }

 public:
  bool updateHaloParticle(const Particle_T &haloParticle) override {
    for (auto &p : particles) {
      if (p.getID() == haloParticle.getID()) {
        p = haloParticle;
        return true;
      }
    }
    return false;
  }
  void rebuildNeighborLists(autopas::TraversalInterface *traversal) override {
    auto trav = dynamic_cast<containers::kokkos::traversal::KokkosTraversalInterface<Particle_T> *>(traversal);
    trav->rebuild(this->particles);
  }

  void deleteHaloParticles() override {
    std::lock_guard<AutoPasLock> lock(_particlesLock);
    std::remove_if(particles.begin(), particles.end(),
                   [](auto &particle) { return !particle.isDummy() || !particle.isOwned(); });
  }

  void deleteAllParticles() override {
    std::lock_guard<AutoPasLock> lock(_particlesLock);
    particles.clear();
  }

  [[nodiscard]] size_t getNumberOfParticles(autopas::IteratorBehavior behavior) const override {
    return std::count_if(particles.begin(), particles.end(),
                         [&](auto &particle) { return behavior.contains(particle); });
  }
  [[nodiscard]] size_t size() const override { return particles.size(); }

  void computeInteractions(autopas::TraversalInterface *traversal) override {
    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override { return _boxMax; }
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override { return _boxMin; }
  [[nodiscard]] double getCutoff() const override { return _cutoff; }
  void setCutoff(double cutoff) override { _cutoff = cutoff; }
  [[nodiscard]] double getVerletSkin() const override { return _verletSkin; }

  [[nodiscard]] size_t getStepsSinceLastRebuild() const override { return _stepsSinceLastRebuild; }

  void setStepsSinceLastRebuild(size_t stepsSinceLastRebuild) override {
    _stepsSinceLastRebuild = stepsSinceLastRebuild;
  }

  [[nodiscard]] double getInteractionLength() const override { return _cutoff + this->getVerletSkin(); }

  [[nodiscard]] std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override { return {}; }
  [[nodiscard]] autopas::TraversalSelectorInfo getTraversalSelectorInfo() const override {
    std::array<double, 3> cellSize = {0, 0, 0};
    for (int axis = 0; axis < 3; ++axis) {
      cellSize[axis] = (_boxMax[axis] - _boxMin[axis]) / (1 << 21);
    }
    return TraversalSelectorInfo(
        // dummy value
        {1, 1, 1}, this->_cutoff , cellSize, 0);
  }

  bool deleteParticle(Particle_T &particle) override {
    std::remove_if(particles.begin(), particles.end(), [&](auto &p) { return particle == p; });
    return true;
  }

  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    for (auto &particle : particles) {
      if (behavior.contains(particle)) {
        if (utils::inBox(particle.getR(), lowerCorner, higherCorner)) {
          forEachLambda(particle);
        }
      }
    }

    if (not(behavior & IteratorBehavior::ownedOrHalo)) {
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
    }
  }

  template <typename Lambda>
  void forEach(Lambda lambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
               typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) {
    for (auto &particle : particles) {
      if (behavior.contains(particle)) {
        lambda(particle);
      }
    }
    if (additionalVectors != nullptr) {
      for (auto &vector : *additionalVectors) {
        for (auto &particle : *vector) {
          if (behavior.contains(particle)) {
            lambda(particle);
          }
        }
      }
    }
  };

  template <typename Lambda>
  void forEach(
      Lambda lambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors = nullptr) const {
    for (const auto &particle : particles) {
      if (behavior.contains(particle)) {
        lambda(particle);
      }
    }
    if (additionalVectors != nullptr) {
      for (const auto &vector : *additionalVectors) {
        for (auto &particle : *vector) {
          if (behavior.contains(particle)) {
            lambda(particle);
          }
        }
      }
    }
  };

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior,
              typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) {
    for (auto &particle : particles) {
      if (behavior.contains(particle)) {
        reduceLambda(particle, result);
      }
    }

    if (additionalVectors != nullptr) {
      for (const auto &vector : *additionalVectors) {
        for (auto &particle : *vector) {
          if (behavior.contains(particle)) {
            reduceLambda(particle, result);
          }
        }
      }
    }
  }

 private:
  AutoPasLock _particlesLock;
  std::vector<Particle_T> particles;

  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
  double _verletSkin;
  int _rebuildFrequency;
  size_t _stepsSinceLastRebuild;
};
}  // namespace autopas

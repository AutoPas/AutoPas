#pragma once

#include "autopas/utils/kokkos/KokkosSoA.h"
#include "autopas/utils/kokkos/KokkosSoAType.h"
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

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {}

 protected:
  void addParticleImpl(const Particle_T &p) override { particles.push_back(p); }
  void addHaloParticleImpl(const Particle_T &haloParticle) override { particles.push_back(haloParticle); }

 public:
  bool updateHaloParticle(const Particle_T &haloParticle) override {
    // TODO

    return false;
  }
  void rebuildNeighborLists(autopas::TraversalInterface *traversal) override {
    auto trav = dynamic_cast<containers::kokkos::traversal::KokkosTraversalInterface<Particle_T> *>(traversal);
    trav.rebuild();
  }

  void deleteHaloParticles() override {
    std::remove_if(particles.begin(), particles.end(), [](auto &particle) {
      return particle._ownershipState == OwnershipState::dummy || particle._ownershipState == OwnershipState::owned;
    });
  }

  void deleteAllParticles() override { particles.clear(); }

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

  [[nodiscard]] size_t getStepsSinceLastRebuild() const override {
    // TODO
    return 0;
  }
  void setStepsSinceLastRebuild(size_t stepsSinceLastRebuild) override {
    // TODO
  }

  [[nodiscard]] double getInteractionLength() const override { return _cutoff + this->getVerletSkin(); }

  [[nodiscard]] std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override { return {}; }
  [[nodiscard]] autopas::TraversalSelectorInfo getTraversalSelectorInfo() const override { return {}; }

  bool deleteParticle(Particle_T &particle) override {
    std::remove_if(particles.begin(), particles.end(), [&](auto &p) { return particle == p; });
    return true;
  }

  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {}

  template <typename Lambda>
  void forEach(Lambda lambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
               typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) {};

  template <typename Lambda>
  void forEach(
      Lambda lambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors = nullptr) const {};

 private:
  std::vector<Particle_T> particles;

  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
  double _verletSkin;
  int _rebuildFrequency;
};
}  // namespace autopas

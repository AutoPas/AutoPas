#pragma once
#include "autopas/utils/kokkos/KokkosSoA.h"
#include "autopas/utils/kokkos/KokkosSoAType.h"

template <typename Space, typename Layout, typename Particle_T>
class KokkosVerletClusterLists : public autopas::ParticleContainerInterface<Particle_T> {
 public:
  KokkosVerletClusterLists(const std::array<double, 3> &box_min, const std::array<double, 3> &box_max,
                           const double cutoff, const double verlet_skin, const int rebuild_frequency)
      : _boxMin(box_min),
        _boxMax(box_max),
        _cutoff(cutoff),
        _verletSkin(verlet_skin),
        _rebuildFrequency(rebuild_frequency) {}

  autopas::CellType getParticleCellTypeEnum() const override {}
  [[nodiscard]] autopas::ContainerOption getContainerType() const override {}
  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {}

 protected:
  void addParticleImpl(const Particle_T &p) override {}
  void addHaloParticleImpl(const Particle_T &haloParticle) override {}

 public:
  bool updateHaloParticle(const Particle_T &haloParticle) override {}
  void rebuildNeighborLists(autopas::TraversalInterface *traversal) override {}
  void deleteHaloParticles() override {}
  void deleteAllParticles() override {}
  [[nodiscard]] size_t getNumberOfParticles(autopas::IteratorBehavior behavior) const override {}
  [[nodiscard]] size_t size() const override {}
  [[nodiscard]] autopas::ContainerIterator<Particle_T, true, false> begin(
      autopas::IteratorBehavior behavior,
      typename autopas::ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors) override {}
  [[nodiscard]] autopas::ContainerIterator<Particle_T, false, false> begin(
      autopas::IteratorBehavior behavior,
      typename autopas::ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors)
      const override {}
  [[nodiscard]] autopas::ContainerIterator<Particle_T, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      autopas::IteratorBehavior behavior,
      typename autopas::ContainerIterator<Particle_T, true, true>::ParticleVecType *additionalVectors) override {}
  [[nodiscard]] autopas::ContainerIterator<Particle_T, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      autopas::IteratorBehavior behavior,
      typename autopas::ContainerIterator<Particle_T, false, true>::ParticleVecType *additionalVectors) const override {
  }
  void computeInteractions(autopas::TraversalInterface *traversal) override {}
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override {}
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override {}
  [[nodiscard]] double getCutoff() const override {}
  void setCutoff(double cutoff) override {}
  [[nodiscard]] double getVerletSkin() const override {}
  [[nodiscard]] size_t getStepsSinceLastRebuild() const override {}
  void setStepsSinceLastRebuild(size_t stepsSinceLastRebuild) override {}
  [[nodiscard]] double getInteractionLength() const override {}
  [[nodiscard]] std::vector<typename autopas::ParticleContainerInterface<Particle_T>::ParticleType> updateContainer(
      bool keepNeighborListsValid) override {}
  [[nodiscard]] autopas::TraversalSelectorInfo getTraversalSelectorInfo() const override {}
  std::tuple<const Particle_T *, size_t, size_t> getParticle(
      size_t cellIndex, size_t particleIndex, autopas::IteratorBehavior iteratorBehavior) const override {}
  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                             autopas::IteratorBehavior iteratorBehavior,
                                                             const std::array<double, 3> &boxMin,
                                                             const std::array<double, 3> &boxMax) const override {}
  bool deleteParticle(Particle_T &particle) override {}
  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {}


 protected:
  using Tuple = typename autopas::utils::kokkos::KokkosSoAType<Space, Layout, typename Particle_T::SoAArraysType>::Type;
  autopas::utils::kokkos::KokkosSoA<Tuple> soa;

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
  double _verletSkin;
  int _rebuildFrequency;
};
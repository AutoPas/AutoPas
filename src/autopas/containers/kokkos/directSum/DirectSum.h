
#pragma once
#include "Kokkos_Core.hpp"
#include "autopas/containers/ParticleContainerInterface.h"

namespace autopas::kokkos {

template <class Particle_T>
class DirectSum : public ParticleContainerInterface<Particle_T> {
  /** Minimal corner of the container box */
  std::array<double, 3> _boxMin;
  /** Maximal corner of the container box */
  std::array<double, 3> _boxMax;
  /** Cutoff-Radius of the container */
  double _cutoff;
  /** Skin Length */
  double _skin;

  Kokkos::View<Particle_T *> _particles;
  Kokkos::View<Particle_T *[6], Kokkos::LayoutLeft> _haloParticles;

  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle_T;

  /**
   * Constructor of DirectSum.
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   */
  DirectSum(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double cutoff, double skin)
      : _boxMin{boxMin}, _boxMax{boxMax}, _cutoff{cutoff}, _skin{skin} {}

  /*
   * Getters and Setters
   */

  /** @copydoc autopas::ParticleContainerInterface::getBoxMax() */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override { return _boxMax; }

  /** @copydoc autopas::ParticleContainerInterface::getBoxMin() */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override { return _boxMin; }

  /** @copydoc autopas::ParticleContainerInterface::getCutoff() */
  [[nodiscard]] double getCutoff() const override { return _cutoff; }

  /** @copydoc autopas::ParticleContainerInterface::setCutoff() */
  void setCutoff(double cutoff) override { _cutoff = cutoff; }

  /** @copydoc autopas::ParticleContainerInterface::getVerletSkin() */
  [[nodiscard]] double getVerletSkin() const override { return _skin; }

  /** @copydoc autopas::ParticleContainerInterface::getInteractionLength() */
  [[nodiscard]] double getInteractionLength() const override { return _skin + _cutoff; }

 public:
  /**
   * The Kokkos DirectSum holds the particles directly
   * @return CellType::FullParticleCell
   */
  CellType getParticleCellTypeEnum() const override { return CellType::FullParticleCell; }

  [[nodiscard]] ContainerOption getContainerType() const override;
  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override;
  bool updateHaloParticle(const Particle_T &haloParticle) override;
  void rebuildNeighborLists(TraversalInterface *traversal) override;
  void deleteHaloParticles() override;
  void deleteAllParticles() override;
  [[nodiscard]] size_t getNumberOfParticles(IteratorBehavior behavior) const override;
  [[nodiscard]] size_t size() const override;
  void computeInteractions(TraversalInterface *traversal) override;

  [[nodiscard]] std::vector<ParticleContainerInterface<Particle_T>::ParticleType> updateContainer(
      bool keepNeighborListsValid) override;
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override;
  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                             IteratorBehavior iteratorBehavior) const override;
  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                             IteratorBehavior iteratorBehavior,
                                                             const std::array<double, 3> &boxMin,
                                                             const std::array<double, 3> &boxMax) const override;
  bool deleteParticle(Particle_T &particle) override;
  bool deleteParticle(size_t cellIndex, size_t particleIndex) override;

 protected:
  void addParticleImpl(const Particle_T &p) override;
  void addHaloParticleImpl(const Particle_T &haloParticle) override;
};
}  // namespace autopas::kokkos
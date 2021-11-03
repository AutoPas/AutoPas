/**
 * @file KokkosCellBasedParticleContainer.h
 *
 * @date 2 Nov 2021
 * @author lgaertner
 */

#pragma once

#include "autopas/containers/KokkosCellBasedParticleContainer.h"

namespace autopas {

/**
 * DirectSumKokkos class.
 */
template <class Particle>
class DirectSumKokkos : public KokkosCellBasedParticleContainer<Particle> {
 public:
  CellType getParticleCellTypeEnum () override {return CellType::KokkosCell;}

  /**
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   */
  DirectSumKokkos(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin)
      : KokkosCellBasedParticleContainer<Particle>(boxMin, boxMax, cutoff, skin) {}


  ~DirectSumKokkos() = default;

  ContainerOption getContainerType() const override {
      return ContainerOption::kokkosDirectSum;
  };

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const Particle &p) override {
    this->_particles.addParticle(p);
  }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addHaloParticleImpl(const Particle &p) override {
    this->_particles.addParticle(p);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const Particle &haloParticle) override {
    // TODO lgaertner
    return false;
  }

  /**
   * @copydoc ParticleContainerInterface::rebuildNeighborLists()
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // TODO lgaertner
  }

  /**
   * @copydoc ParticleContainerInterface::deleteHaloParticles()
   */
  void deleteHaloParticles() override {
    // TODO lgaertner
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    //TODO lgaertner
  }

  std::vector<Particle> updateContainer(bool keepNeighborListsValid) override {
    //TODO lgaertner
    autopas::utils::ExceptionHandler::exception("TODO");
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // TODO lgaertner
    return TraversalSelectorInfo(this->getCellBlock().getCellsPerDimensionWithHalo(), this->getInteractionLength(),
                                 this->getCellBlock().getCellLength(), 0);
  }

 private:
  Kokkos::View<size_t> start;
  Kokkos::View<size_t> cellSize;

  size_t _OWNED{0};
  size_t _HALO{1};
};

}
/**
 * @file KokkosCellBasedParticleContainer.h
 *
 * @date 2 Nov 2021
 * @author lgaertner
 */

#pragma once

#include <autopas/cells/KokkosParticleCell.h>
#include <autopas/containers/ParticleContainerInterface.h>
#include <autopas/kokkosContainers/ParticleView.h>

namespace autopas {

// consider multiple inheritance or delegation to avoid virtual call to Functor
/**
 * The KokkosCellBasedParticleContainer class stores particles in a ParticleView and provides
 * methods to iterate over its particles.
 * @tparam Particle Class for the particles
 */
template <class Particle>
class KokkosCellBasedParticleContainer : public ParticleContainerInterface<Particle> {
  using ParticleCell = KokkosParticleCell<Particle>;

 public:
  /**
   * Constructor of KokkosCellBasedParticleContainer
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   */
  KokkosCellBasedParticleContainer(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax,
                                   const double cutoff, const double skin, const size_t cellcount)
      : _particles(), _cells("cellsView", cellcount), _boxMin(boxMin), _boxMax(boxMax), _cutoff(cutoff), _skin(skin) {
    // change to create_view to avoid persisting changes when host= devicememoryspace, at cost of storing additional
    // view
    _cellsHostMirror = Kokkos::create_mirror_view(_cells);
  }

  /**
   * Destructor of KokkosCellBasedParticleContainer.
   */
  ~KokkosCellBasedParticleContainer() override = default;

  /**
   * Delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  KokkosCellBasedParticleContainer(const KokkosCellBasedParticleContainer &obj) = delete;

  /**
   * Delete the copy assignment operator to prevent unwanted copies
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  KokkosCellBasedParticleContainer &operator=(const KokkosCellBasedParticleContainer &other) = delete;

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMax()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override final { return _boxMax; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setBoxMax()
   */
  void setBoxMax(const std::array<double, 3> &boxMax) override final { _boxMax = boxMax; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMin()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override final { return _boxMin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setBoxMin()
   */
  void setBoxMin(const std::array<double, 3> &boxMin) override final { _boxMin = boxMin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getCutoff()
   */
  [[nodiscard]] double getCutoff() const override final { return _cutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(double cutoff) override final { _cutoff = cutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getSkin()
   */
  [[nodiscard]] double getSkin() const override final { return _skin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setSkin()
   */
  void setSkin(double skin) override final { _skin = skin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getInteractionLength()
   */
  [[nodiscard]] double getInteractionLength() const override final { return _cutoff + _skin; }

  /**
   * Deletes all particles from the container.
   */
  void deleteAllParticles() override { _particles.deleteAll(); }

  /**
   * Fetches all cells from the DeviceMemory. !!!Use carefully, as changes are not persisent when using differing host-
   * and device memory, but DO persist when host- = devicememoryspace. Fix by changing create_mirror_view to
   * create_view!!!
   */
  Kokkos::View<KokkosParticleCell<Particle> *> getCellsHost() {
    Kokkos::resize(_cellsHostMirror, _cells.size());
    Kokkos::deep_copy(_cellsHostMirror, _cells);
    return _cellsHostMirror;
  }

  /**
   * Get the number of particles saved in the container.
   * @return Number of particles in the container.
   */
  [[nodiscard]] unsigned long getNumParticles() const override { return _particles.getSize(); }

  ParticleIteratorWrapper<Particle, true> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) override {
    autopas::utils::ExceptionHandler::exception("Deprecated function 'begin' used. Use forEach instead!");
    return ParticleIteratorWrapper<Particle, true>();
  }

  ParticleIteratorWrapper<Particle, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const override {
    autopas::utils::ExceptionHandler::exception("Deprecated function 'begin' used. Use forEach instead!");
    return ParticleIteratorWrapper<Particle, false>();
  }

  ParticleIteratorWrapper<Particle, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                            const std::array<double, 3> &higherCorner,
                                                            IteratorBehavior behavior) override {
    autopas::utils::ExceptionHandler::exception(
        "Deprecated function 'getRegionIterator' used. Use forEachInRegion instead!");
    return ParticleIteratorWrapper<Particle, true>();
  }

  ParticleIteratorWrapper<Particle, false> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                             const std::array<double, 3> &higherCorner,
                                                             IteratorBehavior behavior) const override {
    autopas::utils::ExceptionHandler::exception(
        "Deprecated function 'getRegionIterator' used. Use forEachInRegion instead!");
    return ParticleIteratorWrapper<Particle, false>();
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticles()
   */
  void updateHaloParticles(const std::vector<Particle> &haloParticles, std::vector<Particle> &notUpdatedParticles) override {
    size_t numHaloParticles = haloParticles.size();

    //Kokkos: There are only two separate views defined if host and device memory space differ.
    Kokkos::View<Particle *> haloParticlesDevice("", numHaloParticles);
    typename Kokkos::View<Particle *>::HostMirror haloParticlesHost = Kokkos::create_mirror_view(haloParticlesDevice);

    Kokkos::View<bool *> wasUpdatedDevice("", numHaloParticles);
    typename Kokkos::View<bool *>::HostMirror wasUpdatedHost = Kokkos::create_mirror_view(wasUpdatedDevice);

    Kokkos::parallel_for(Kokkos::RangePolicy<>(0, numHaloParticles), [=] (const size_t i) {
      Particle p = haloParticles[i];
      p.setOwnershipState(OwnershipState::halo);
      haloParticlesHost[i] = p;
    });
    Kokkos::fence();

    //Kokkos: This is a noop if host and memory space are the same.
    Kokkos::deep_copy(haloParticlesDevice, haloParticlesHost);

    Kokkos::RangePolicy<> rangePolicy(0ul, numHaloParticles);
    auto particles = _particles.getParticles();
    Kokkos::parallel_for("update halo particles", rangePolicy, KOKKOS_LAMBDA (const size_t i) {
      Particle haloParticle = haloParticlesDevice[i];
      size_t begin, end;
      if (_isDirty) {
        begin = 0ul;
        end = _particles.getSize();
      } else {
        size_t possibleCell = assignCellToParticle(haloParticle);
        begin = _cells[possibleCell].begin;
        end = _cells[possibleCell].cellSize + begin;
      }

      bool updated = false;
      for(size_t j = begin; j < end; j++) {
        Particle p = particles(j);
        if (p.getID() == haloParticle.getID()) {
          auto distanceVec = autopas::utils::ArrayMath::sub(p.getR(), haloParticle.getR());
          auto distanceSqr = autopas::utils::ArrayMath::dot(distanceVec, distanceVec);
          if (distanceSqr < this->getSkin() * this->getSkin()) {
            // found the particle -> update, set updated to true and break for loop
            _particles.getParticles()[j] = haloParticle;
            updated = true;
            break;
          }
        }
      }

      wasUpdatedDevice[i] = updated;
    });
    Kokkos::fence();

    Kokkos::deep_copy(wasUpdatedHost, wasUpdatedDevice);

    for (size_t i = 0; i < numHaloParticles; i++) {
      if (not wasUpdatedDevice[i]) {
        Particle p = haloParticles[i];
        notUpdatedParticles.push_back(p);
      }
    }
  }

  virtual size_t assignCellToParticle(Particle &p) = 0;

  bool getIsDirty() override {
    return _isDirty;
  }

 protected:
  /**
   * Managed view of all particles. This view is then sorted according to desired layout and particle properties.
   */
  ParticleView<Particle> _particles;

  Kokkos::View<KokkosParticleCell<Particle> *> _cells;
  typename Kokkos::View<KokkosParticleCell<Particle> *>::HostMirror _cellsHostMirror;

  bool _isDirty;

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
  double _skin;
};
}  // namespace autopas
#ifndef AUTOPAS_AUTOPAS_H
#define AUTOPAS_AUTOPAS_H

#include <iostream>
#include <memory>
#include "autopasIncludes.h"

/**
 * The AutoPas class is intended to be the main point of Interaction for the
 * user. It puts a layer of abstraction over the container and handles the
 * autotuning. (->TODO)
 * @tparam Particle Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle, class ParticleCell>
class AutoPas {
 public:
  enum ContainerOption { directSum, linkedCells };
  enum DataLayoutOption { aos, soa };

  /**
   * Initialize container
   */
  void init(ContainerOption containerOption, std::array<double, 3> boxSize,
            double cutoff) {
    switch (containerOption) {
      case directSum: {
        container = std::unique_ptr<ContainerType>(new autopas::DirectSum<Particle, ParticleCell>(
            {0., 0., 0.}, boxSize, cutoff));
        break;
      }
      case linkedCells: {
        container = std::unique_ptr<ContainerType>(new autopas::LinkedCells<Particle, ParticleCell>(
            {0., 0., 0.}, boxSize, cutoff));
        break;
      }
      default: {
        std::cerr << "AutoPas.init(): Unknown container Option! "
                  << containerOption << std::endl;
        exit(1);
      }
    }
  }

  /**
   * Returns a pointer to the actual container.
   * @return container
   */
  // TODO: do we need the whole container functionality available to the outside
  // or whould wrapper over some container functions be better?
  autopas::ParticleContainer<Particle, ParticleCell> *getContainer() const {
    return container.get();
  }

  /**
   * Adds a particle to the container.
   * @param p Reference to the particle to be added
   */
  void addParticle(Particle &p) { container->addParticle(p); }

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @tparam useSoA Bool to decide if SoA or AoS should be used.
   * @param f Functor that describes the pair-potential
   */
  template <bool useSoA>
  void iteratePairwise(autopas::Functor<Particle, ParticleCell> *f) {
    if (useSoA) {
      container->iteratePairwiseSoA(f);
    } else {
      container->iteratePairwiseAoS(f);
    }
  }

  /**
   * iterate over all particles by using
   * for(auto iter = container.begin(); iter.isValid(); ++iter)
   * @return iterator to the first particle
   */
  autopas::ParticleIterator<Particle, ParticleCell> begin() {
      return container->begin();
  }


  /**
   * iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner,
   * highCorner);iter.isValid();++iter)
   * @param lowerCorner lower corner of the region
   * @param higherCorner higher corner of the region
   * @return iterator to iterate over all particles in a specific region
   */
  autopas::RegionParticleIterator<Particle, ParticleCell> getRegionIterator(
          std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner) {
      return container->getRegionIterator(lowerCorner, higherCorner);
  }
 private:
  typedef autopas::ParticleContainer<Particle, ParticleCell> ContainerType;
  std::unique_ptr<ContainerType> container;
};

#endif  // AUTOPAS_AUTOPAS_H

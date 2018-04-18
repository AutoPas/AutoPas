#ifndef AUTOPAS_AUTOPAS_H
#define AUTOPAS_AUTOPAS_H

#include <iostream>
#include <memory>
#include "autopasIncludes.h"

namespace autopas {
/**
 * Possible choices for the particle container type.
 */
enum ContainerOption { directSum, linkedCells };

/**
 * Provides a way to iterate over the possible choices of ContainerOption.
 */
std::array<ContainerOption, 2> possibleContainerOptions = {ContainerOption::directSum, ContainerOption::linkedCells};

/**
 * Possible Choices for the particle data layout.
 */
enum DataLayoutOption { aos, soa };
}  // namespace autopas

/**
 * The AutoPas class is intended to be the main point of Interaction for the
 * user. It puts a layer of abstraction over the container and handles the
 * autotuning.
 * @todo autotuning
 * @tparam Particle Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle, class ParticleCell>
class AutoPas {
 public:
  /**
   * Initialize the particle container.
   *
   * For possible container choices see AutoPas::ContainerOption.
   *
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param containerOption Type of the container.
   */
  void init(std::array<double, 3> boxMin, std::array<double, 3> boxMax,
            double cutoff, autopas::ContainerOption containerOption) {
    switch (containerOption) {
      case autopas::directSum: {
        container = std::unique_ptr<ContainerType>(
            new autopas::DirectSum<Particle, ParticleCell>(boxMin, boxMax,
                                                           cutoff));
        break;
      }
      case autopas::linkedCells: {
        container = std::unique_ptr<ContainerType>(
            new autopas::LinkedCells<Particle, ParticleCell>(boxMin, boxMax,
                                                             cutoff));
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
   * @overload
   *
   * @param boxSize Size of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param containerOption Type of the container.
   */
  void init(std::array<double, 3> boxSize, double cutoff,
            autopas::ContainerOption containerOption) {
    init({0, 0, 0}, boxSize, cutoff, containerOption);
  }

  /**
   * Returns a pointer to the actual container.
   * @todo do we need the whole container functionality available to the outside
   * @return container
   */
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
   * @param f Functor that describes the pair-potential
   * @param dataLayoutOption useSoA Bool to decide if SoA or AoS should be used.
   */
  void iteratePairwise(autopas::Functor<Particle, ParticleCell> *f,
                       autopas::DataLayoutOption dataLayoutOption) {
    switch (dataLayoutOption) {
      case autopas::aos: {
        container->iteratePairwiseAoS(f);
        break;
      }
      case autopas::soa: {
        container->iteratePairwiseSoA(f);
      }
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

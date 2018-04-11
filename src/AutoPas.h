#ifndef AUTOPAS_AUTOPAS_H
#define AUTOPAS_AUTOPAS_H

#include <iostream>
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
        container = new autopas::DirectSum<Particle, ParticleCell>(
            {0., 0., 0.}, boxSize, cutoff);
        break;
      }
      case linkedCells: {
        container = new autopas::LinkedCells<Particle, ParticleCell>(
            {0., 0., 0.}, boxSize, cutoff);
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
    return container;
  }

 private:
  autopas::ParticleContainer<Particle, ParticleCell> *container;
};

#endif  // AUTOPAS_AUTOPAS_H

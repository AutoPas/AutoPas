/**
 * @file KokkosParticleCell.h
 * @date 20.10.2021
 * @author lgaertner
 */

#include "autopas/cells/ParticleCell.h"

namespace autopas {

/**
 * Dummy particle cell to use CellBlock3D with Kokkos Containers
 * @tparam Particle
 */
template<class Particle>
class KokkosParticleCell : public ParticleCell<Particle> {

};
}
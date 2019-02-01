/**
 * @file CellFunctorCuda.h
 *
 * @date 13.12.2018
 * @author jspahl
 */

#pragma once

#include "LJFunctorCuda.h"
#include "cuda_runtime.h"

namespace autopas {

/**
 * A cell functor. This functor is build from the normal Functor of the template
 * type ParticleFunctor. It is an internal object to handle interactions between
 * two cells of particles.
 * @todo: currently always used newton3!
 * @tparam Particle
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class Particle, class ParticleCell, class ParticleFunctor, bool useSoA, bool useNewton3 = true>
class CellFunctorCuda {
 public:
  /**
   * The constructor of CellFunctor
   * @param f the particlefunctor which should be used for the interaction.
   */
  explicit CellFunctorCuda(ParticleFunctor *f) : _functor(f) {}

  /**
   * process the interactions inside one cell
   * @param cell all pairwise interactions of particles inside this cell are
   * calculated
   */
  void processCell(ParticleCell &cell);

  /**
   * process the interactions between the particles of cell1 with particles of
   * cell2.
   * @param cell1
   * @param cell2
   */
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2);

 private:
  ParticleFunctor *_functor;
};

template <class Particle, class ParticleCell, class ParticleFunctor, bool useSoA, bool useNewton3>
void CellFunctorCuda<Particle, ParticleCell, ParticleFunctor, useSoA, useNewton3>::processCell(ParticleCell &cell) {
  static_assert(not useNewton3, "Newton3 not available for Cuda Cell Functor");

  if (useSoA) {
    _functor->SoAFunctorNoN3(cell._particleSoABuffer.getNumParticles(), cell._particleSoABufferDevice);

    return;
  } else {
    _functor->AoSFunctorNoN3(cell.numParticles(), cell._particlesDevice);
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, bool useSoA, bool useNewton3>
void CellFunctorCuda<Particle, ParticleCell, ParticleFunctor, useSoA, useNewton3>::processCellPair(
    ParticleCell &cell1, ParticleCell &cell2) {
  if (useSoA) {
    return;
  } else {
    _functor->AoSFunctorNoN3(cell1.numParticles(), cell2.numParticles(), cell1._particlesDevice,
                             cell2._particlesDevice);
  }
}

}  // namespace autopas

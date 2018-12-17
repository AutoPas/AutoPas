/**
 * @file CellFunctorCuda.h
 *
 * @date 13.12.2018
 * @author jspahl
 */

#pragma once

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
template<class Particle, class ParticleCell, class ParticleFunctor, bool useSoA,
		bool useNewton3 = true>
class CellFunctorCuda {
public:
	/**
	 * The constructor of CellFunctor
	 * @param f the particlefunctor which should be used for the interaction.
	 */
	explicit CellFunctorCuda(ParticleFunctor *f) :
			_functor(f) {
	}

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

template<class Particle, class ParticleCell, class ParticleFunctor, bool useSoA,
		bool useNewton3>
void CellFunctorCuda<Particle, ParticleCell, ParticleFunctor, useSoA, useNewton3>::processCell(
		ParticleCell &cell) {
	static_assert(not useNewton3, "Newton3 not available for Cuda Cell Functor");

	if (not useSoA) {
		const size_t c0_size = cell._particles.size();

		double3* cell0_r = new double3[c0_size];
		double3* cell0_f = new double3[c0_size];
		double3* d_cell0_r, d_cell0_f;

		for (size_t i = 0; i < c0_size; ++i) {
			cell0_r[i] = _functor.serializePosition(cell._particles[i]);
			cell0_f[i] = {0,0,0};
		}
		cudaError_t e;
		e = cudaMalloc((void **) &d_cell0_r, sizeof(double3) * c0_size);
		e = cudaMalloc((void **) &d_cell0_f, sizeof(double3) * c0_size);

		e = cudaMemcpy(d_cell0_r, cell0_r, c0_size, cudaMemcpyHostToDevice);

		_functor->CudaFunctor(c0_size, d_cell0_r, d_cell0_f, useNewton3);

		e = cudaMemcpy(cell0_f, d_cell0_f, c0_size, cudaMemcpyDeviceToHost);
		cudaFree(d_cell0_r);
		cudaFree(d_cell0_f);

		for (size_t i = 0; i < c0_size; ++i) {
			updateForce(cell._particles[i], cell0_f[i]);;
		}

		return;
	}
}

template<class Particle, class ParticleCell, class ParticleFunctor, bool useSoA,
		bool useNewton3>
void CellFunctorCuda<Particle, ParticleCell, ParticleFunctor, useSoA, useNewton3>::processCellPair(
		ParticleCell &cell1, ParticleCell &cell2) {
}

} // namespace autopas


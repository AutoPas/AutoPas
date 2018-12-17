/**
 * @file LJFunctorCuda.h
 *
 * @date 14.12.2018
 * @author jspahl
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"

namespace autopas {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * @tparam Particle The type of particle.
 * @tparam ParticleCell The type of particlecell.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, class ParticleCell, bool calculateGlobals = false, bool relevantForTuning = true>
class LJFunctorCuda : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorCuda() = delete;

  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift.
   * @param cutoff
   * @param epsilon
   * @param sigma
   * @param shift
   * @param lowCorner Lower corner of the local simulation domain.
   * @param highCorner Upper corner of the local simulation domain.
   * @param duplicatedCalculation Defines whether duplicated calculations are happening across processes / over the
   * simulation boundary. e.g. eightShell: false, fullShell: true.
   */
  explicit LJFunctorCuda(double cutoff, double epsilon, double sigma, double shift,
                     std::array<double, 3> lowCorner = {0., 0., 0.}, std::array<double, 3> highCorner = {0., 0., 0.},
                     bool duplicatedCalculation = true)
      : _cutoffsquare{cutoff * cutoff},
        _epsilon24{epsilon * 24.0},
        _sigmasquare{sigma * sigma},
        _shift6{shift * 6.0},
        _upotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _duplicatedCalculations{duplicatedCalculation},
        _lowCorner{lowCorner},
        _highCorner{highCorner},
        _postProcessed{false} {
    if (calculateGlobals and duplicatedCalculation) {
      if (lowCorner == highCorner) {
        throw utils::ExceptionHandler::AutoPasException(
            "Please specify the lowCorner and highCorner properly if calculateGlobals and duplicatedCalculation are "
            "set to true.");
      }
    }
    if (calculateGlobals) {
      _aosThreadData.resize(autopas_get_max_threads());
    }
  }

  struct constants{
	  double cutoffsquare;
	  double epsilon24;
	  double sigmasquare;
  };

  bool isRelevantForTuning() override { return relevantForTuning; }

  __host__
  double3 serializePosition(Particle& p){
	  double3 d_r;
	  d_r.x = p.getR()[0];
	  d_r.y = p.getR()[1];
	  d_r.z = p.getR()[2];

	  return d_r;
  }

  __host__
  void updateForce(Particle& p, double3& d_f){
	  p.getF()[0] = d_f.x;
	  p.getF()[1] = d_f.y;
	  p.getF()[2] = d_f.z;
  }

  void CudaFunctor(int N, double3* rs, double3* fs, bool newton3);


 private:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, upotSum{0.} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      upotSum = 0.;
    }

    // variables
    std::array<double, 3> virialSum;
    double upotSum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[4];
  };
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  double _cutoffsquare, _epsilon24, _sigmasquare, _shift6;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _upotSum;
  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // bool that defines whether duplicate calculations are happening
  bool _duplicatedCalculations;
  // lower and upper corner of the domain of the current process
  std::array<double, 3> _lowCorner, _highCorner;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

};  // class LJFunctor


__device__
double3 bodyBodyF(double3 i,
		double3 j, double3 fi, bool newton3) {
  double drx = i.x - j.x;
  double dry = i.y - j.y;
  double drz = i.z - j.z;

  double dr2 = drx* drx + dry * dry + drz* drz;
  double cutoffsquare = 1.0;
  double sigmasquare = 1.0;
  double epsilon24 = 1.0;

  if (dr2 > cutoffsquare) {
  	return fi;
  }

  double invdr2 = 1. / dr2;
  double lj6 = sigmasquare * invdr2;
  lj6 = lj6 * lj6 * lj6;
  double lj12 = lj6 * lj6;
  double lj12m6 = lj12 - lj6;
  double fac = epsilon24 * (lj12 + lj12m6) * invdr2;

  fi.x += drx * fac;
  fi.y += dry * fac;
  fi.z += drz * fac;

  if (newton3) {
    // only if we use newton 3 here, we want to
    //j.subF(f);
  }
  return fi;
}

__device__
double3 tile_calculation(double3 myposition, double3 f, bool newton3){
	  int i;
	  extern __shared__ double3 all_pos[];
	  for(i = 0; i < blockDim.x; ++i){
		  f = bodyBodyF(myposition, all_pos[i], f, newton3);
	  }
	  return f;
}

__global__
void AoSFunctorCuda(int N, int p, double3* rs, double3* fs, bool newton3){
	 extern __shared__ double3 all_pos[];
	 int i, tile;
	 int tid = blockIdx.x * blockDim.x + threadIdx.x;
	 double3 myposition = rs[tid];
	 double3 myf = fs[tid];
	 for(i = 0, tile = 0; i < N; i+=p, ++tile){
		 int idx = tile * blockDim.x + threadIdx.x;
		 all_pos[threadIdx.x] = rs[idx];
		 __syncthreads();
		 myf = tile_calculation(myposition, myf, newton3);
		 __syncthreads();
	 }

	 fs[tid] = myf;
}

template <class Particle, class ParticleCell, bool calculateGlobals, bool relevantForTuning>
void LJFunctorCuda<Particle, ParticleCell, calculateGlobals, relevantForTuning>::CudaFunctor(int N, double3* rs, double3* fs, bool newton3){
	constants c;
	c.cutoffsquare = _cutoffsquare;
	c.epsilon24 = _epsilon24;
	c.sigmasquare = _sigmasquare;

	__constant__ constants global_constants;

	cudaMemcpyToSymbol(global_constants, c, sizeof(constants));

	AoSFunctorCuda<<<1,1, N * sizeof(double3)>>>(N, 1, rs, fs, newton3);
}

}  // namespace autopas

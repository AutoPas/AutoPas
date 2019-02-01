/**
 * @file LJFunctorCuda.h
 *
 * @date 14.12.2018
 * @author jspahl
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"

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
                         std::array<double, 3> lowCorner = {0., 0., 0.},
                         std::array<double, 3> highCorner = {0., 0., 0.}, bool duplicatedCalculation = true)
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

    loadConstants(_cutoffsquare, _epsilon24, _sigmasquare);
  }

  bool isRelevantForTuning() override { return relevantForTuning; }

  void AoSFunctorNoN3(int N, double *particles) { AoSFunctorNoN3Wrapper(N, particles); }

  void AoSFunctorNoN3(int N, int M, double *particles1, double *particles2) {
    AoSFunctorNoN3Wrapper(N, M, particles1, particles2);
  }

  void SoAFunctorNoN3(int N, typename Particle::SoADevice &device_handle) {
    SoAFunctorNoN3Wrapper(N, device_handle.posX, device_handle.posY, device_handle.posZ, device_handle.forceX,
                          device_handle.forceY, device_handle.forceZ);
  }

  void deviceAoSLoader(ParticleCell &cell, double *device_buffer) {
    size_t num_particles = cell.numParticles();
    double *particles = new double[num_particles * 6];

    cudaMalloc((void **)device_buffer, sizeof(double) * 6 * num_particles);

    auto cellIter = cell.begin();
    for (size_t i = 0; cellIter.isValid(); ++cellIter) {
      particles[i++] = cellIter->getR()[0];
      particles[i++] = cellIter->getR()[1];
      particles[i++] = cellIter->getR()[2];
      particles[i++] = cellIter->getF()[0];
      particles[i++] = cellIter->getF()[1];
      particles[i++] = cellIter->getF()[2];
    }

    cudaMemcpy(device_buffer, particles, num_particles * 6, cudaMemcpyHostToDevice);
    delete[] particles;
  }

  void deviceAoSExtractor(ParticleCell &cell, double *device_buffer) {
    size_t num_particles = cell.numParticles();
    double *particles = new double[num_particles * 6];

    cudaMemcpy(particles, device_buffer, num_particles * 6, cudaMemcpyDeviceToHost);
    cudaFree(device_buffer);

    auto cellIter = cell.begin();

    for (size_t i = 3; cellIter.isValid(); i += 4, ++cellIter) {
      cellIter->addF({particles[i], particles[++i], particles[++i]});
    }

    delete[] particles;
  }

  void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa, typename Particle::SoADevice &device_handle) {
    size_t size = soa.getNumParticles();
    if (size == 0) return;
    // cudaMalloc((void **)device_handle.ids, sizeof(typename Particle::SoADevice::ids) * size);
    cudaMalloc((void **)&device_handle.posX, sizeof(typename Particle::SoADevice::posX) * size);
    cudaMalloc((void **)&device_handle.posY, sizeof(typename Particle::SoADevice::posY) * size);
    cudaMalloc((void **)&device_handle.posZ, sizeof(typename Particle::SoADevice::posZ) * size);
    cudaMalloc((void **)&device_handle.forceX, sizeof(typename Particle::SoADevice::forceX) * size);
    cudaMalloc((void **)&device_handle.forceY, sizeof(typename Particle::SoADevice::forceY) * size);
    cudaMalloc((void **)&device_handle.forceZ, sizeof(typename Particle::SoADevice::forceZ) * size);

    cudaError_t e;
    e = cudaMemcpy(device_handle.posX, soa.template begin<Particle::AttributeNames::posX>(), size,
                   cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.posY, soa.template begin<Particle::AttributeNames::posY>(), size,
                   cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.posZ, soa.template begin<Particle::AttributeNames::posZ>(), size,
                   cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.forceX, soa.template begin<Particle::AttributeNames::forceX>(), size,
                   cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.forceY, soa.template begin<Particle::AttributeNames::forceY>(), size,
                   cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.forceZ, soa.template begin<Particle::AttributeNames::forceZ>(), size,
                   cudaMemcpyHostToDevice);
  }

  void deviceSoAExtractor(::autopas::SoA<SoAArraysType> &soa, typename Particle::SoADevice &device_handle) {
    size_t size = soa.getNumParticles();
    if (size == 0) return;

    cudaMemcpy(soa.template begin<Particle::AttributeNames::forceX>(), device_handle.forceX, size,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(soa.template begin<Particle::AttributeNames::forceY>(), device_handle.forceY, size,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(soa.template begin<Particle::AttributeNames::forceZ>(), device_handle.forceZ, size,
               cudaMemcpyDeviceToHost);
    cudaFree(device_handle.forceX);
    cudaFree(device_handle.forceY);
    cudaFree(device_handle.forceZ);
    cudaFree(device_handle.posX);
    cudaFree(device_handle.posY);
    cudaFree(device_handle.posZ);
  }

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

}  // namespace autopas

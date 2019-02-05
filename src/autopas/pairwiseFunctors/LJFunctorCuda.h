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
    AoSFunctorNoN3PairWrapper(N, M, particles1, particles2);
  }

  void SoAFunctorNoN3(int N, typename Particle::SoADevice &device_handle) {
    SoAFunctorNoN3Wrapper(N, device_handle.posX, device_handle.posY, device_handle.posZ, device_handle.forceX,
                          device_handle.forceY, device_handle.forceZ);
  }

  void SoAFunctorNoN3(int N, typename Particle::SoADevice &device_handle1, int M,
                      typename Particle::SoADevice &device_handle2) {
    SoAFunctorNoN3PairWrapper(N, device_handle1.posX, device_handle1.posY, device_handle1.posZ, device_handle1.forceX,
                              device_handle1.forceY, device_handle1.forceZ, M, device_handle2.posX, device_handle2.posY,
                              device_handle2.posZ);
  }

  void deviceAoSLoader(ParticleCell &cell, double **device_buffer) {
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
    cudaMemcpy(*device_buffer, particles, num_particles * 6 * sizeof(double), cudaMemcpyHostToDevice);

    delete[] particles;
  }

  void deviceAoSExtractor(ParticleCell &cell, double **device_buffer) {
    size_t num_particles = cell.numParticles();
    double *particles = new double[num_particles * 6];

    cudaMemcpy(particles, *device_buffer, num_particles * 6 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(*device_buffer);

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
    cudaMalloc((void **)&device_handle.posX, sizeof(Particle::SoADevice::posX) * size);
    cudaMalloc((void **)&device_handle.posY, sizeof(Particle::SoADevice::posY) * size);
    cudaMalloc((void **)&device_handle.posZ, sizeof(Particle::SoADevice::posZ) * size);
    cudaMalloc((void **)&device_handle.forceX, sizeof(Particle::SoADevice::forceX) * size);
    cudaMalloc((void **)&device_handle.forceY, sizeof(Particle::SoADevice::forceY) * size);
    cudaMalloc((void **)&device_handle.forceZ, sizeof(Particle::SoADevice::forceZ) * size);

    cudaError_t e;
    e = cudaMemcpy(device_handle.posX, soa.template begin<Particle::AttributeNames::posX>(),
                   sizeof(Particle::SoADevice::posX) * size, cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.posY, soa.template begin<Particle::AttributeNames::posY>(),
                   sizeof(Particle::SoADevice::posY) * size, cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.posZ, soa.template begin<Particle::AttributeNames::posZ>(),
                   sizeof(Particle::SoADevice::posZ) * size, cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.forceX, soa.template begin<Particle::AttributeNames::forceX>(),
                   sizeof(Particle::SoADevice::forceX) * size, cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.forceY, soa.template begin<Particle::AttributeNames::forceY>(),
                   sizeof(Particle::SoADevice::forceY) * size, cudaMemcpyHostToDevice);
    e = cudaMemcpy(device_handle.forceZ, soa.template begin<Particle::AttributeNames::forceZ>(),
                   sizeof(Particle::SoADevice::forceZ) * size, cudaMemcpyHostToDevice);
  }

  void deviceSoAExtractor(::autopas::SoA<SoAArraysType> &soa, typename Particle::SoADevice &device_handle) {
    size_t size = soa.getNumParticles();
    if (size == 0) return;

    cudaMemcpy(soa.template begin<Particle::AttributeNames::forceX>(), device_handle.forceX,
               sizeof(Particle::SoADevice::forceX) * size, cudaMemcpyDeviceToHost);
    cudaMemcpy(soa.template begin<Particle::AttributeNames::forceY>(), device_handle.forceY,
               sizeof(Particle::SoADevice::forceY) * size, cudaMemcpyDeviceToHost);
    cudaMemcpy(soa.template begin<Particle::AttributeNames::forceZ>(), device_handle.forceZ,
               sizeof(Particle::SoADevice::forceZ) * size, cudaMemcpyDeviceToHost);
    cudaFree(device_handle.forceX);
    cudaFree(device_handle.forceY);
    cudaFree(device_handle.forceZ);
    cudaFree(device_handle.posX);
    cudaFree(device_handle.posY);
    cudaFree(device_handle.posZ);
  }

  /**
   * SoALoader
   * @param cell
   * @param soa
   * @param offset
   */
  AUTOPAS_FUNCTOR_SOALOADER(
      cell, soa, offset,
      // @todo it is probably better to resize the soa only once, before calling
      // SoALoader (verlet-list only)
      soa.resizeArrays(offset + cell.numParticles());

      if (cell.numParticles() == 0) return;

      unsigned long *const __restrict__ idptr = soa.template begin<Particle::AttributeNames::id>();
      double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();
      double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
      double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
      double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

      auto cellIter = cell.begin();
      // load particles in SoAs
      for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) {
        idptr[i] = cellIter->getID();
        xptr[i] = cellIter->getR()[0];
        yptr[i] = cellIter->getR()[1];
        zptr[i] = cellIter->getR()[2];
        fxptr[i] = cellIter->getF()[0];
        fyptr[i] = cellIter->getF()[1];
        fzptr[i] = cellIter->getF()[2];
      })
  /**
   * soaextractor
   * @param cell
   * @param soa
   * @param offset
   */
  AUTOPAS_FUNCTOR_SOAEXTRACTOR(
      cell, soa, offset,
      // body start
      if (soa.getNumParticles() == 0) return;

      auto cellIter = cell.begin();

#ifndef NDEBUG
      unsigned long *const __restrict__ idptr = soa.template begin<Particle::AttributeNames::id>();
#endif

      double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
      double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
      double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

      for (size_t i = offset; cellIter.isValid(); ++i, ++cellIter) {
        assert(idptr[i] == cellIter->getID());
        cellIter->setF({fxptr[i], fyptr[i], fzptr[i]});
      })

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

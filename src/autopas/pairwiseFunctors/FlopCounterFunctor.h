/**
 * @file FlopCounterFunctor.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
/**
 * This class helps in getting the number of performed floating point
 * operations. It is a functor that only calculated the amount of floating point
 * operations.
 * @todo this class currently is limited to the following case:
 *  - constant cutoff radius
 *  - constant amount of floating point operations for one kernel call (distance
 * < cutoff)
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle, class ParticleCell>
class FlopCounterFunctor : public Functor<Particle, ParticleCell> {
  typedef typename Particle::SoAArraysType SoAArraysType;

 public:
  bool isRelevantForTuning() override { return false; }

  /**
   * constructor of FlopCounterFunctor
   * @param cutoffRadius the cutoff radius
   */
  explicit FlopCounterFunctor<Particle, ParticleCell>(double cutoffRadius)
      : autopas::Functor<Particle, ParticleCell>(),
        _cutoffSquare(cutoffRadius * cutoffRadius),
        _distanceCalculations(0ul),
        _kernelCalls(0ul) {}

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    auto dr = ArrayMath::sub(i.getR(), j.getR());
    double dr2 = ArrayMath::dot(dr, dr);
#ifdef AUTOPAS_OPENMP
#pragma omp critical
#endif
    {
      ++_distanceCalculations;

      if (dr2 <= _cutoffSquare) ++_kernelCalls;
    };
  }


  void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ x1ptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y1ptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z1ptr = soa.template begin<Particle::AttributeNames::posZ>();

    for (unsigned int i = 0; i < soa.getNumParticles(); ++i) {
      unsigned long distanceCalculationsAcc = 0;
      unsigned long kernelCallsAcc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : kernelCallsAcc, distanceCalculationsAcc)
      for (unsigned int j = i + 1; j < soa.getNumParticles(); ++j) {
        ++distanceCalculationsAcc;

        const double drx = x1ptr[i] - x1ptr[j];
        const double dry = y1ptr[i] - y1ptr[j];
        const double drz = z1ptr[i] - z1ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 <= _cutoffSquare) ++kernelCallsAcc;
      }
#ifdef AUTOPAS_OPENMP
#pragma omp critical
#endif
      {
        _distanceCalculations += distanceCalculationsAcc;
        _kernelCalls += kernelCallsAcc;
      }
    }
  }

  void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3) override {
    double *const __restrict__ x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    double *const __restrict__ x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      unsigned long distanceCalculationsAcc = 0;
      unsigned long kernelCallsAcc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : kernelCallsAcc, distanceCalculationsAcc)
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        ++distanceCalculationsAcc;

        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 <= _cutoffSquare) {
          ++kernelCallsAcc;
        }
      }
#ifdef AUTOPAS_OPENMP
#pragma omp critical
#endif
      {
        _distanceCalculations += distanceCalculationsAcc;
        _kernelCalls += kernelCallsAcc;
      }
    }
  }

  void SoAFunctor(SoA<SoAArraysType> &soa,
                  const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom,
                  size_t iTo, bool newton3) override {
    auto numParts = soa.getNumParticles();

    if (numParts == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    for (size_t i = iFrom; i < iTo; ++i) {
      const size_t listSizeI = neighborList[i].size();
      const size_t *const __restrict__ currentList = neighborList[i].data();

      // this is a magic number, that should correspond to at least
      // vectorization width*N have testet multiple sizes:
      // 4: small speedup compared to AoS
      // 8: small speedup compared to AoS
      // 12: small but best speedup compared to Aos
      // 16: smaller speedup
      // in theory this is a variable, we could auto-tune over...
#ifdef __AVX512F__
      // use a multiple of 8 for avx
      const size_t vecsize = 16;
#else
      // for everything else 12 is faster
      const size_t vecsize = 12;
#endif
      size_t joff = 0;

      // if the size of the verlet list is larger than the given size vecsize,
      // we will use a vectorized version.
      if (listSizeI >= vecsize) {
        alignas(64) std::array<double, vecsize> xtmp{}, ytmp{}, ztmp{}, xArr{}, yArr{}, zArr{};
        // broadcast of the position of particle i
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          xtmp[tmpj] = xptr[i];
          ytmp[tmpj] = yptr[i];
          ztmp[tmpj] = zptr[i];
        }
        // loop over the verlet list from 0 to x*vecsize
        for (; joff < listSizeI - vecsize + 1; joff += vecsize) {
          unsigned long distanceCalculationsAcc = 0;
          unsigned long kernelCallsAcc = 0;
          // in each iteration we calculate the interactions of particle i with
          // vecsize particles in the neighborlist of particle i starting at
          // particle joff

          // gather position of particle j
#pragma omp simd safelen(vecsize)
          for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
            xArr[tmpj] = xptr[currentList[joff + tmpj]];
            yArr[tmpj] = yptr[currentList[joff + tmpj]];
            zArr[tmpj] = zptr[currentList[joff + tmpj]];
          }

          // do omp simd with reduction of the interaction
#pragma omp simd reduction(+ : kernelCallsAcc, distanceCalculationsAcc) safelen(vecsize)
          for (size_t j = 0; j < vecsize; j++) {
            ++distanceCalculationsAcc;
            const double drx = xtmp[j] - xArr[j];
            const double dry = ytmp[j] - yArr[j];
            const double drz = ztmp[j] - zArr[j];

            const double drx2 = drx * drx;
            const double dry2 = dry * dry;
            const double drz2 = drz * drz;

            const double dr2 = drx2 + dry2 + drz2;

            const unsigned long mask = (dr2 <= _cutoffSquare) ? 1 : 0;

            kernelCallsAcc += mask;
          }
#ifdef AUTOPAS_OPENMP
#pragma omp critical
#endif
          {
            _distanceCalculations += distanceCalculationsAcc;
            _kernelCalls += kernelCallsAcc;
          }
        }
      }
      unsigned long distanceCalculationsAcc = 0;
      unsigned long kernelCallsAcc = 0;
      // this loop goes over the remainder and uses no optimizations
      for (size_t jNeighIndex = joff; jNeighIndex < listSizeI; ++jNeighIndex) {
        size_t j = neighborList[i][jNeighIndex];
        if (i == j) continue;

        ++distanceCalculationsAcc;
        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 <= _cutoffSquare) {
          ++kernelCallsAcc;
        }
      }
#ifdef AUTOPAS_OPENMP
#pragma omp critical
#endif
      {
        _distanceCalculations += distanceCalculationsAcc;
        _kernelCalls += kernelCallsAcc;
      }
    }
  }

  void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3) override {
#if defined(AUTOPAS_CUDA)
    // estimate flops on GPU
    size_t size = device_handle.template get<Particle::AttributeNames::posX>().size();
    _distanceCalculations += size * size;
    _kernelCalls += size * size;
#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }

  void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                   CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3) override {
#if defined(AUTOPAS_CUDA)
    // estimate flops on GPU
    size_t size1 = device_handle1.template get<Particle::AttributeNames::posX>().size();
    size_t size2 = device_handle2.template get<Particle::AttributeNames::posX>().size();

    _distanceCalculations += size1 * size2;
    _kernelCalls += size1 * size2;
#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }

  /**
   * @copydoc Functor::deviceSoALoader
   */
  void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
                       CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {
#if defined(AUTOPAS_CUDA)
    size_t size = soa.getNumParticles();
    if (size == 0) return;

    device_handle.template get<Particle::AttributeNames::posX>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::posX>());
    device_handle.template get<Particle::AttributeNames::posY>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::posY>());
    device_handle.template get<Particle::AttributeNames::posZ>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::posZ>());

    device_handle.template get<Particle::AttributeNames::forceX>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::forceX>());
    device_handle.template get<Particle::AttributeNames::forceY>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::forceY>());
    device_handle.template get<Particle::AttributeNames::forceZ>().copyHostToDevice(
        size, soa.template begin<Particle::AttributeNames::forceZ>());
#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }

  /**
   * @copydoc Functor::deviceSoAExtractor
   */
  void deviceSoAExtractor(::autopas::SoA<SoAArraysType> &soa,
                          CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {
#if defined(AUTOPAS_CUDA)
    size_t size = soa.getNumParticles();
    if (size == 0) return;

    device_handle.template get<Particle::AttributeNames::forceX>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceX>());
    device_handle.template get<Particle::AttributeNames::forceY>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceY>());
    device_handle.template get<Particle::AttributeNames::forceZ>().copyDeviceToHost(
        size, soa.template begin<Particle::AttributeNames::forceZ>());
#else
    utils::ExceptionHandler::exception("AutoPas was compiled without CUDA support!");
#endif
  }

  AUTOPAS_FUNCTOR_SOALOADER(
      cell, soa, offset,
      // body start
      soa.resizeArrays(offset + cell.numParticles());

      if (cell.numParticles() == 0) return;

      double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

      auto cellIter = cell.begin();
      // load particles in SoAs
      for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) {
        xptr[i] = cellIter->getR()[0];
        yptr[i] = cellIter->getR()[1];
        zptr[i] = cellIter->getR()[2];
      })

  /**
   * Empty SoAExtractor.
   * Nothing to be done yet.
   */
  AUTOPAS_FUNCTOR_SOAEXTRACTOR(, , , )

  /**
   * get the hit rate of the pair-wise interaction, i.e. the ratio of the number
   * of kernel calls compared to the number of distance calculations
   * @return the hit rate
   */
  double getHitRate() { return static_cast<double>(_kernelCalls) / static_cast<double>(_distanceCalculations); }

  /**
   * get the total number of flops
   * @param numFlopsPerKernelCall
   * @return
   */
  double getFlops(unsigned long numFlopsPerKernelCall) const {
    const double distFlops = numFlopsPerDistanceCalculation * static_cast<double>(_distanceCalculations);
    const double kernFlops = numFlopsPerKernelCall * static_cast<double>(_kernelCalls);
    return distFlops + kernFlops;
  }

  /**
   * get the number of calculated distance operations
   * @return
   */
  unsigned long getDistanceCalculations() const { return _distanceCalculations; }

  /**
   * get the number of kernel calls, i.e. the number of pairs of particles with
   * a distance not larger than the cutoff
   * @return
   */
  unsigned long getKernelCalls() const { return _kernelCalls; }

  /**
   * number of flops for one distance calculation.
   * 3 sub + 3 square + 2 add
   */
  static constexpr double numFlopsPerDistanceCalculation = 8.0;

 private:
  double _cutoffSquare;
  unsigned long _distanceCalculations, _kernelCalls;
};

}  // namespace autopas

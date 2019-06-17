/**
 * @file LJFunctorAVX.h
 *
 * @date 17 Jan 2018
 * @author F. Gratl
 */
#pragma once

#include <immintrin.h>
#include <array>
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This Version is implemented using AVX intrinsics.
 * @tparam Particle The type of particle.
 * @tparam ParticleCell The type of particlecell.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, class ParticleCell, bool calculateGlobals = false, bool relevantForTuning = true>
class LJFunctorAVX : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorAVX() = delete;

  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift.
   * @param cutoff
   * @param epsilon
   * @param sigma
   * @param shift
   * @param duplicatedCalculation Defines whether duplicated calculations are happening across processes / over the
   * simulation boundary. e.g. eightShell: false, fullShell: true.
   */
  explicit LJFunctorAVX(double cutoff, double epsilon, double sigma, double shift, bool duplicatedCalculation = false)
#ifdef __AVX__
      : _one{_mm256_set1_pd(1.)},
        _masks{
            _mm256_set_epi64x(0, 0, 0, -1),
            _mm256_set_epi64x(0, 0, -1, -1),
            _mm256_set_epi64x(0, -1, -1, -1),
        },
        _cutoffsquare{_mm256_set1_pd(cutoff * cutoff)},
        _epsilon24{_mm256_set1_pd(epsilon * 24.0)},
        _sigmasquare{_mm256_set1_pd(sigma * sigma)},
        _shift6{shift * 6.0},
        _upotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _duplicatedCalculations{duplicatedCalculation},
        _postProcessed{false} {
    if (calculateGlobals or duplicatedCalculation) {
      utils::ExceptionHandler::exception("Calculation of global values and duplicated calculations not supported.");
    }
    //    if (calculateGlobals) {
    //      _aosthreaddata.resize(autopas_get_max_threads());
    //    }
  }
#else
      : _one{0}, _masks{0, 0, 0}, _cutoffsquare{0}, _epsilon24{0}, _sigmasquare{0} {
    utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
  }
#endif

  bool isRelevantForTuning() override { return relevantForTuning; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    utils::ExceptionHandler::exception("LJFunctorAVX.AoSFunctor() not implemented!");
  }

  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa, bool newton3)
   * This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {
#ifdef __AVX__
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    // @TODO: Globlas
    //    bool isHaloCell1 = false;
    //    bool isHaloCell2 = false;

    // reverse outer loop s.th. inner loop always beginns at aligned array start
    // typecast to detect underflow
    for (size_t i = soa.getNumParticles() - 1; (long)i >= 0; --i) {
      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();

      const __m256d x1 = _mm256_broadcast_sd(&xptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&yptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&zptr[i]);

      // floor soa numParticles to multiple of vecLength
      unsigned int j = 0;
      for (; j < (i & ~(vecLength - 1)); j += 4) {
        SoAKernel<true, false>(j, x1, y1, z1, xptr, yptr, zptr, fxptr, fyptr, fzptr, fxacc, fyacc, fzacc);
      }
      const int rest = (int)(i & (vecLength - 1));
      if (rest > 0)
        SoAKernel<true, true>(j, x1, y1, z1, xptr, yptr, zptr, fxptr, fyptr, fzptr, fxacc, fyacc, fzacc, rest);

      // horizontally reduce fDacc to sumfD
      const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
      const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

      const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
      const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

      const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
      const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

      const union {
        __m128d reg;
        double arr[2];
      } sumfxfyVEC = {.reg = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh)};
      const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC.arr[0];
      const double sumfy = sumfxfyVEC.arr[1];
      const double sumfz = _mm_cvtsd_f64(sumfzVEC);

      fxptr[i] += sumfx;
      fyptr[i] += sumfy;
      fzptr[i] += sumfz;
    }
#endif
  }

 private:
  template <bool newton3, bool masked>
  inline void SoAKernel(size_t j, const __m256d &x1, const __m256d &y1, const __m256d &z1,
                        double *const __restrict__ &x2ptr, double *const __restrict__ &y2ptr,
                        double *const __restrict__ &z2ptr, double *const __restrict__ &fx2ptr,
                        double *const __restrict__ &fy2ptr, double *const __restrict__ &fz2ptr, __m256d &fxacc,
                        __m256d &fyacc, __m256d &fzacc, const unsigned int rest = 0) {
#ifdef __AVX__
    const __m256d x2 = masked ? _mm256_maskload_pd(&x2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&x2ptr[j]);
    const __m256d y2 = masked ? _mm256_maskload_pd(&y2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&y2ptr[j]);
    const __m256d z2 = masked ? _mm256_maskload_pd(&z2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&z2ptr[j]);

    const __m256d drx = _mm256_sub_pd(x1, x2);
    const __m256d dry = _mm256_sub_pd(y1, y2);
    const __m256d drz = _mm256_sub_pd(z1, z2);

    const __m256d drx2 = _mm256_mul_pd(drx, drx);
    const __m256d dry2 = _mm256_mul_pd(dry, dry);
    const __m256d drz2 = _mm256_mul_pd(drz, drz);

    const __m256d dr2PART = _mm256_add_pd(drx2, dry2);
    const __m256d dr2 = _mm256_add_pd(dr2PART, drz2);

    // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
    // signaling = throw error if NaN is encountered
    // dr2 <= _cutoffsquare ? 0xFFFFFFFFFFFFFFFF : 0
    const __m256d cutoffMask = _mm256_cmp_pd(dr2, _cutoffsquare, _CMP_LE_OS);

    const __m256d invdr2 = _mm256_div_pd(_one, dr2);
    const __m256d lj2 = _mm256_mul_pd(_sigmasquare, invdr2);
    const __m256d lj4 = _mm256_mul_pd(lj2, lj2);
    const __m256d lj6 = _mm256_mul_pd(lj2, lj4);
    const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
    const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
    const __m256d lj12m6alj12 = _mm256_add_pd(lj12m6, lj12);
    const __m256d lj12m6alj12e = _mm256_mul_pd(lj12m6alj12, _epsilon24);
    const __m256d fac = _mm256_mul_pd(lj12m6alj12e, invdr2);

    const __m256d facMasked = masked
                                  ? _mm256_and_pd(fac, _mm256_and_pd(cutoffMask, _mm256_castsi256_pd(_masks[rest - 1])))
                                  : _mm256_and_pd(fac, cutoffMask);

    const __m256d fx = _mm256_mul_pd(drx, facMasked);
    const __m256d fy = _mm256_mul_pd(dry, facMasked);
    const __m256d fz = _mm256_mul_pd(drz, facMasked);

    fxacc = _mm256_add_pd(fxacc, fx);
    fyacc = _mm256_add_pd(fyacc, fy);
    fzacc = _mm256_add_pd(fzacc, fz);

    // if newton 3 is used subtract fD from particle j
    if (newton3) {
      const __m256d fx2 = masked ? _mm256_maskload_pd(&fx2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fx2ptr[j]);
      const __m256d fy2 = masked ? _mm256_maskload_pd(&fy2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fy2ptr[j]);
      const __m256d fz2 = masked ? _mm256_maskload_pd(&fz2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fz2ptr[j]);

      const __m256d fx2new = _mm256_sub_pd(fx2, fx);
      const __m256d fy2new = _mm256_sub_pd(fy2, fy);
      const __m256d fz2new = _mm256_sub_pd(fz2, fz);

      masked ? _mm256_maskstore_pd(&fx2ptr[j], _masks[rest - 1], fx2new) : _mm256_store_pd(&fx2ptr[j], fx2new);
      masked ? _mm256_maskstore_pd(&fy2ptr[j], _masks[rest - 1], fy2new) : _mm256_store_pd(&fy2ptr[j], fy2new);
      masked ? _mm256_maskstore_pd(&fz2ptr[j], _masks[rest - 1], fz2new) : _mm256_store_pd(&fz2ptr[j], fz2new);
    }
#endif
  }

 public:
  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3)
   */
  void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, const bool newton3) override {
#ifdef __AVX__
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

    double *const __restrict__ x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    double *const __restrict__ x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    double *const __restrict__ fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    // @@TODO: Globlas
    //    bool isHaloCell1 = false;
    //    bool isHaloCell2 = false;

    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();

      const __m256d x1 = _mm256_broadcast_sd(&x1ptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&y1ptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&z1ptr[i]);

      // floor soa2 numParticles to multiple of vecLength
      if (newton3) {
        unsigned int j = 0;
        for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {
          SoAKernel<true, false>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc);
        }
        const int rest = (int)(soa2.getNumParticles() & (vecLength - 1));
        if (rest > 0)
          SoAKernel<true, true>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc, rest);
      } else {
        unsigned int j = 0;
        for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {
          SoAKernel<false, false>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc);
        }
        const int rest = (int)(soa2.getNumParticles() & (vecLength - 1));
        if (rest > 0)
          SoAKernel<false, true>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc, rest);
      }

      // horizontally reduce fDacc to sumfD
      const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
      const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

      const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
      const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

      const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
      const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

      const union {
        __m128d reg;
        double arr[2];
      } sumfxfyVEC = {.reg = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh)};
      const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC.arr[0];
      const double sumfy = sumfxfyVEC.arr[1];
      const double sumfz = _mm_cvtsd_f64(sumfzVEC);

      fx1ptr[i] += sumfx;
      fy1ptr[i] += sumfy;
      fz1ptr[i] += sumfz;
    }
#endif
  }

  // clang-format off
  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa, const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
  // clang-format on
  void SoAFunctor(SoA<SoAArraysType> &soa,
                  const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom,
                  size_t iTo, bool newton3) override {
    utils::ExceptionHandler::exception("Verlet SoA functor not implemented!");
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

  /**
   * Get the number of flops used per kernel call. This should count the
   * floating point operations needed for two particles that lie within a cutoff
   * radius.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall() {
    // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply
    // scale) sum Forces: 6 (forces) kernel total = 12 + 6 = 18
    return 18ul;
  }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() override {
    _upotSum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
    }
  }

  /**
   * Accumulates global values, e.g. upot and virial.
   * @param newton3
   */
  void endTraversal(bool newton3) override {
    if (_postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _upotSum += _aosThreadData[i].upotSum;
      _virialSum = ArrayMath::add(_virialSum, _aosThreadData[i].virialSum);
    }
    if (not newton3) {
      // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2 here.
      _upotSum *= 0.5;
      _virialSum = ArrayMath::mulScalar(_virialSum, 0.5);
    }
    // we have always calculated 6*upot, so we divide by 6 here!
    _upotSum /= 6.;
    _postProcessed = true;
  }

  /**
   * Get the potential Energy
   * @return the potential Energy
   */
  double getUpot() {
    if (not calculateGlobals) {
      throw utils::ExceptionHandler::AutoPasException(
          "Trying to get upot even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw utils::ExceptionHandler::AutoPasException("Cannot get upot, because endTraversal was not called.");
    }
    return _upotSum;
  }

  /**
   * Get the virial
   * @return the virial
   */
  double getVirial() {
    if (not calculateGlobals) {
      throw utils::ExceptionHandler::AutoPasException(
          "Trying to get virial even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw utils::ExceptionHandler::AutoPasException("Cannot get virial, because endTraversal was not called.");
    }
    return _virialSum[0] + _virialSum[1] + _virialSum[2];
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

  //  double _cutoffsquare, _epsilon24, _sigmasquare, _shift6;
  const __m256d _one;
  //  const std::vector<__m256i> _masks;
  const __m256i _masks[3];
  const __m256d _cutoffsquare;
  const __m256d _epsilon24;
  const __m256d _sigmasquare;
  double _shift6;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _upotSum;
  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // bool that defines whether duplicate calculations are happening
  bool _duplicatedCalculations;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // number of double values that fit into a vector register.
  constexpr static size_t vecLength = 4;

};  // namespace autopas
}  // namespace autopas

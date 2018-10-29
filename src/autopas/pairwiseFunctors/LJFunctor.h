/**
 * @file LJFunctor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * A functor to handle lennard-jones interactions between two particles
 * (molecules).
 * @tparam Particle the type of particle
 * @tparam ParticleCell the type of particlecell
 * @tparam calculateGlobals defines whether the global values are to be calculated (energy, virial)
 */
template <class Particle, class ParticleCell, bool calculateGlobals = false>
class LJFunctor : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctor() = delete;

  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift
   * @param cutoff
   * @param epsilon
   * @param sigma
   * @param shift
   * @param lowCorner lower corner of the local simulation domain
   * @param highCorner upper corner of the local simulation domain
   * @param duplicatedCalculation defines whether duplicated calculations are happening across processes / over the
   * simulation boundary. e.g. eightShell: false, fullShell: true
   */
  explicit LJFunctor(double cutoff, double epsilon, double sigma, double shift,
                     std::array<double, 3> lowCorner = {0., 0., 0.}, std::array<double, 3> highCorner = {0., 0., 0.},
                     bool duplicatedCalculation = true)
      : _virialSum{0., 0., 0.} {
    _cutoffsquare = cutoff * cutoff;
    _epsilon24 = epsilon * 24.0;
    _sigmasquare = sigma * sigma;
    _shift6 = shift * 6.0;
    _upotSum = 0.;
    _duplicatedCalculations = duplicatedCalculation;
    _lowCorner = lowCorner;
    _highCorner = highCorner;
    _postProcessed = false;
  }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    auto dr = ArrayMath::sub(i.getR(), j.getR());
    double dr2 = ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffsquare) return;

    double invdr2 = 1. / dr2;
    double lj6 = _sigmasquare * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = _epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = ArrayMath::mulScalar(dr, fac);
    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }
    if (calculateGlobals) {
      auto virial = ArrayMath::mul(dr, f);
      double upot = _epsilon24 * lj12m6 + _shift6;

      if (newton3 and _duplicatedCalculations) {
        double upot_half = upot * 0.5;
        auto virial_half = ArrayMath::mulScalar(virial, 0.5);

        if (autopas::inBox(i.getR(), _lowCorner, _highCorner)) {
          _upotSum += upot_half;
          _virialSum = ArrayMath::add(_virialSum, virial_half);
        }
        if (autopas::inBox(j.getR(), _lowCorner, _highCorner)) {
          _upotSum += upot_half;
          _virialSum = ArrayMath::add(_virialSum, virial_half);
        }
      } else {
        _upotSum += upot;
        _virialSum = ArrayMath::add(_virialSum, virial);
      }
    }
  }

  void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    for (unsigned int i = 0; i < soa.getNumParticles(); ++i) {
      double fxacc = 0;
      double fyacc = 0;
      double fzacc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = i + 1; j < soa.getNumParticles(); ++j) {
        // if (i == j) continue;

        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        const double mask = (dr2 > _cutoffsquare) ? 0. : 1.;

        const double invdr2 = 1. / dr2;
        const double lj2 = _sigmasquare * invdr2;
        const double lj6 = lj2 * lj2 * lj2;
        const double lj12 = lj6 * lj6;
        const double lj12m6 = lj12 - lj6;
        const double fac = _epsilon24 * (lj12 + lj12m6) * invdr2 * mask;

        const double fx = drx * fac;
        const double fy = dry * fac;
        const double fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
        if (newton3) {
          fxptr[j] -= fx;
          fyptr[j] -= fy;
          fzptr[j] -= fz;
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
  }

  void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3) override {
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

    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      double fxacc = 0;
      double fyacc = 0;
      double fzacc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        const double mask = (dr2 > _cutoffsquare) ? 0. : 1.;

        const double invdr2 = 1. / dr2;
        const double lj2 = _sigmasquare * invdr2;
        const double lj6 = lj2 * lj2 * lj2;
        const double lj12 = lj6 * lj6;
        const double lj12m6 = lj12 - lj6;
        const double fac = _epsilon24 * (lj12 + lj12m6) * invdr2 * mask;

        const double fx = drx * fac;
        const double fy = dry * fac;
        const double fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
        if (newton3) {
          fx2ptr[j] -= fx;
          fy2ptr[j] -= fy;
          fz2ptr[j] -= fz;
        }
      }

      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
  }

  // clang-format off
  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa, const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
   * @note if you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly
   */
  // clang-format on
  void SoAFunctor(SoA<SoAArraysType> &soa,
                  const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom,
                  size_t iTo, bool newton3) override {
    auto numParts = soa.getNumParticles();
    AutoPasLog(debug, "SoAFunctorVerlet: {}", soa.getNumParticles());

    if (numParts == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    for (size_t i = iFrom; i < iTo; ++i) {
      double fxacc = 0;
      double fyacc = 0;
      double fzacc = 0;
      const size_t listSizeI = neighborList[i].size();
      const size_t *const __restrict__ currentList = neighborList[i].data();

      // this is a magic number, that should correspond to at least
      // vectorization width*N have testet multiple sizes:
      // 4: does not give a speedup, slower than original AoSFunctor
      // 8: small speedup compared to AoS
      // 12: highest speedup compared to Aos
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
        alignas(64) std::array<double, vecsize> xtmp, ytmp, ztmp, xArr, yArr, zArr, fxArr, fyArr, fzArr;
        // broadcast of the position of particle i
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          xtmp[tmpj] = xptr[i];
          ytmp[tmpj] = yptr[i];
          ztmp[tmpj] = zptr[i];
        }
        // loop over the verlet list from 0 to x*vecsize
        for (; joff < listSizeI - vecsize + 1; joff += vecsize) {
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
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc) safelen(vecsize)
          for (size_t j = 0; j < vecsize; j++) {
            // const size_t j = currentList[jNeighIndex];

            const double drx = xtmp[j] - xArr[j];
            const double dry = ytmp[j] - yArr[j];
            const double drz = ztmp[j] - zArr[j];

            const double drx2 = drx * drx;
            const double dry2 = dry * dry;
            const double drz2 = drz * drz;

            const double dr2 = drx2 + dry2 + drz2;

            const double mask = (dr2 <= _cutoffsquare) ? 1. : 0.;

            const double invdr2 = 1. / dr2 * mask;
            const double lj2 = _sigmasquare * invdr2;
            const double lj6 = lj2 * lj2 * lj2;
            const double lj12 = lj6 * lj6;
            const double lj12m6 = lj12 - lj6;
            const double fac = _epsilon24 * (lj12 + lj12m6) * invdr2;

            const double fx = drx * fac;
            const double fy = dry * fac;
            const double fz = drz * fac;

            fxacc += fx;
            fyacc += fy;
            fzacc += fz;
            if (newton3) {
              fxArr[j] = fx;
              fyArr[j] = fy;
              fzArr[j] = fz;
            }
          }
          // scatter the forces to where they belong, this is only needed for newton3
          if (newton3) {
#pragma omp simd safelen(vecsize)
            for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
              const size_t j = currentList[joff + tmpj];
              fxptr[j] -= fxArr[tmpj];
              fyptr[j] -= fyArr[tmpj];
              fzptr[j] -= fzArr[tmpj];
            }
          }
        }
      }
      // this loop goes over the remainder and uses no optimizations
      for (size_t jNeighIndex = joff; jNeighIndex < listSizeI; ++jNeighIndex) {
        size_t j = neighborList[i][jNeighIndex];
        if (i == j) continue;

        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 > _cutoffsquare) continue;

        const double invdr2 = 1. / dr2;
        const double lj2 = _sigmasquare * invdr2;
        const double lj6 = lj2 * lj2 * lj2;
        const double lj12 = lj6 * lj6;
        const double lj12m6 = lj12 - lj6;
        const double fac = _epsilon24 * (lj12 + lj12m6) * invdr2;

        const double fx = drx * fac;
        const double fy = dry * fac;
        const double fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
        if (newton3) {
          fxptr[j] -= fx;
          fyptr[j] -= fy;
          fzptr[j] -= fz;
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
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
   * get the number of flops used per kernel call. This should count the
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
  void resetGlobalValues() {
    _upotSum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
  }

  /**
   * postprocesses global values, e.g. upot and virial
   * @param newton3
   */
  void postProcessGlobalValues(bool newton3) {
    if (not newton3) {
      _upotSum *= 0.5;
      _virialSum = ArrayMath::mulScalar(0.5, _virialSum);
    }
  }

  /**
   * Get the potential Energy
   * @return the potential Energy
   */
  double getUpot() {
    if (not _postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
          "not yet postprocessed, please call postProcessGlobalValues first.");
    }
    return _upotSum;
  }

  /**
   * Get the virial
   * @return the virial
   */
  double getVirial() {
    if (not _postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
          "not yet postprocessed, please call postProcessGlobalValues first.");
    }
    return _virialSum[0] + _virialSum[1] + _virialSum[2];
  }

 private:
  double _cutoffsquare, _epsilon24, _sigmasquare, _shift6;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _upotSum;
  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // bool that defines whether duplicate calculations are happening
  bool _duplicatedCalculations;
  // lower and upper corner of the domain of the current process
  std::array<double, 3> _lowCorner, _highCorner;

  //
  bool _postProcessed;

};  // class LJFunctor
}  // namespace autopas

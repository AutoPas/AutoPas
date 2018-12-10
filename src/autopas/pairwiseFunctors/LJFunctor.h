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
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * @tparam Particle The type of particle.
 * @tparam ParticleCell The type of particlecell.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, class ParticleCell, bool calculateGlobals = false, bool relevantForTuning = true>
class LJFunctor : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctor() = delete;

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
  explicit LJFunctor(double cutoff, double epsilon, double sigma, double shift,
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

  bool isRelevantForTuning() override { return relevantForTuning; }

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

      const int threadnum = autopas_get_thread_num();
      if (_duplicatedCalculations) {
        // for non-newton3 the division is in the post-processing step.
        if (newton3) {
          upot *= 0.5;
          virial = ArrayMath::mulScalar(virial, 0.5);
        }
        if (autopas::utils::inBox(i.getR(), _lowCorner, _highCorner)) {
          _aosThreadData[threadnum].upotSum += upot;
          _aosThreadData[threadnum].virialSum = ArrayMath::add(_aosThreadData[threadnum].virialSum, virial);
        }
        // for non-newton3 the second particle will be considered in a separate calculation
        if (newton3 and autopas::utils::inBox(j.getR(), _lowCorner, _highCorner)) {
          _aosThreadData[threadnum].upotSum += upot;
          _aosThreadData[threadnum].virialSum = ArrayMath::add(_aosThreadData[threadnum].virialSum, virial);
        }
      } else {
        // for non-newton3 we divide by 2 only in the postprocess step!
        _aosThreadData[threadnum].upotSum += upot;
        _aosThreadData[threadnum].virialSum = ArrayMath::add(_aosThreadData[threadnum].virialSum, virial);
      }
    }
  }

  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa, bool newton3)
   * This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();
    // the local redeclaration of the following values helps the auto-generation of various compilers.
    const double cutoffsquare = _cutoffsquare, epsilon24 = _epsilon24, sigmasquare = _sigmasquare, shift6 = _shift6;
    if (calculateGlobals) {
      // Checks if the cell is a halo cell, if it is, we skip it.
      // We cannot do this in normal cases (where we do not calculate globals), as _lowCorner and _highCorner are not
      // set. (as of 23.11.2018)
      bool isHaloCell = false;
      isHaloCell |= xptr[0] < _lowCorner[0] || xptr[0] >= _highCorner[0];
      isHaloCell |= yptr[0] < _lowCorner[1] || yptr[0] >= _highCorner[1];
      isHaloCell |= zptr[0] < _lowCorner[2] || zptr[0] >= _highCorner[2];
      if (isHaloCell) {
        return;
      }
    }

    double upotSum = 0.;
    double virialSumX = 0.;
    double virialSumY = 0.;
    double virialSumZ = 0.;

    for (unsigned int i = 0; i < soa.getNumParticles(); ++i) {
      double fxacc = 0.;
      double fyacc = 0.;
      double fzacc = 0.;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, upotSum, virialSumX, virialSumY, virialSumZ)
      for (unsigned int j = i + 1; j < soa.getNumParticles(); ++j) {
        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        const double mask = (dr2 > cutoffsquare) ? 0. : 1.;

        const double invdr2 = 1. / dr2;
        const double lj2 = sigmasquare * invdr2;
        const double lj6 = lj2 * lj2 * lj2;
        const double lj12 = lj6 * lj6;
        const double lj12m6 = lj12 - lj6;
        const double fac = epsilon24 * (lj12 + lj12m6) * invdr2 * mask;

        const double fx = drx * fac;
        const double fy = dry * fac;
        const double fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        // newton 3
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;

        if (calculateGlobals) {
          const double virialx = drx * fx;
          const double virialy = dry * fy;
          const double virialz = drz * fz;
          const double upot = (epsilon24 * lj12m6 + shift6) * mask;

          // these calculations assume that this functor is not called for halo cells!
          upotSum += upot;
          virialSumX += virialx;
          virialSumY += virialy;
          virialSumZ += virialz;
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
    if (calculateGlobals) {
      const int threadnum = autopas_get_thread_num();
      // if newton3 is false, then we divide by 2 later on, so we multiply by two here (very hacky, but needed for
      // AoS)
      _aosThreadData[threadnum].upotSum += upotSum * (newton3 ? 1 : 2);
      _aosThreadData[threadnum].virialSum[0] += virialSumX * (newton3 ? 1 : 2);
      _aosThreadData[threadnum].virialSum[1] += virialSumY * (newton3 ? 1 : 2);
      _aosThreadData[threadnum].virialSum[2] += virialSumZ * (newton3 ? 1 : 2);
    }
  }

  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3)
   */
  void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, const bool newton3) override {
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

    bool isHaloCell1 = false;
    bool isHaloCell2 = false;
    // Checks whether the cells are halo cells.
    // This check cannot be done if _lowCorner and _highCorner are not set. So we do this only if calculateGlobals is
    // defined. (as of 23.11.2018)
    if (calculateGlobals) {
      isHaloCell1 |= x1ptr[0] < _lowCorner[0] || x1ptr[0] >= _highCorner[0];
      isHaloCell1 |= y1ptr[0] < _lowCorner[1] || y1ptr[0] >= _highCorner[1];
      isHaloCell1 |= z1ptr[0] < _lowCorner[2] || z1ptr[0] >= _highCorner[2];
      isHaloCell2 |= x2ptr[0] < _lowCorner[0] || x2ptr[0] >= _highCorner[0];
      isHaloCell2 |= y2ptr[0] < _lowCorner[1] || y2ptr[0] >= _highCorner[1];
      isHaloCell2 |= z2ptr[0] < _lowCorner[2] || z2ptr[0] >= _highCorner[2];

      // This if is commented out because the AoS vs SoA test would fail otherwise. Even though it is physically
      // correct!
      /*if(_duplicatedCalculations and isHaloCell1 and isHaloCell2){
        return;
      }*/
    }
    double upotSum = 0.;
    double virialSumX = 0.;
    double virialSumY = 0.;
    double virialSumZ = 0.;

    const double cutoffsquare = _cutoffsquare, epsilon24 = _epsilon24, sigmasquare = _sigmasquare, shift6 = _shift6;
    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      double fxacc = 0;
      double fyacc = 0;
      double fzacc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, upotSum, virialSumX, virialSumY, virialSumZ)
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        const double mask = (dr2 > cutoffsquare) ? 0. : 1.;

        const double invdr2 = 1. / dr2;
        const double lj2 = sigmasquare * invdr2;
        const double lj6 = lj2 * lj2 * lj2;
        const double lj12 = lj6 * lj6;
        const double lj12m6 = lj12 - lj6;
        const double fac = epsilon24 * (lj12 + lj12m6) * invdr2 * mask;

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

        if (calculateGlobals) {
          double virialx = drx * fx;
          double virialy = dry * fy;
          double virialz = drz * fz;
          double upot = (epsilon24 * lj12m6 + shift6) * mask;

          upotSum += upot;
          virialSumX += virialx;
          virialSumY += virialy;
          virialSumZ += virialz;
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
    if (calculateGlobals) {
      double energyfactor = 1.;
      if (_duplicatedCalculations) {
        // if we have duplicated calculations, i.e., we calculate interactions multiple times, we have to take care
        // that we do not add the energy multiple times!
        energyfactor = isHaloCell1 ? 0. : 1.;
        if (newton3) {
          energyfactor += isHaloCell2 ? 0. : 1.;
          energyfactor *= 0.5;  // we count the energies partly to one of the two cells!
        }
      }
      const int threadnum = autopas_get_thread_num();

      _aosThreadData[threadnum].upotSum += upotSum * energyfactor;
      _aosThreadData[threadnum].virialSum[0] += virialSumX * energyfactor;
      _aosThreadData[threadnum].virialSum[1] += virialSumY * energyfactor;
      _aosThreadData[threadnum].virialSum[2] += virialSumZ * energyfactor;
    }
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
    auto numParts = soa.getNumParticles();
    AutoPasLog(debug, "SoAFunctorVerlet: {}", soa.getNumParticles());

    if (numParts == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const double cutoffsquare = _cutoffsquare, epsilon24 = _epsilon24, sigmasquare = _sigmasquare, shift6 = _shift6;

    double upotSum = 0.;
    double virialSumX = 0.;
    double virialSumY = 0.;
    double virialSumZ = 0.;

    for (size_t i = iFrom; i < iTo; ++i) {
      double fxacc = 0;
      double fyacc = 0;
      double fzacc = 0;
      const size_t listSizeI = neighborList[i].size();
      const size_t *const __restrict__ currentList = neighborList[i].data();

      // checks whether particle 1 is in the domain box, unused if _duplicatedCalculations is false!
      bool inbox1 = false;
      if (_duplicatedCalculations) {  // only for duplicated calculations we need this value
        inbox1 = autopas::utils::inBox({xptr[i], yptr[i], zptr[i]}, _lowCorner, _highCorner);
      }

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
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, upotSum, virialSumX, virialSumY, virialSumZ) safelen(vecsize)
          for (size_t j = 0; j < vecsize; j++) {
            // const size_t j = currentList[jNeighIndex];

            const double drx = xtmp[j] - xArr[j];
            const double dry = ytmp[j] - yArr[j];
            const double drz = ztmp[j] - zArr[j];

            const double drx2 = drx * drx;
            const double dry2 = dry * dry;
            const double drz2 = drz * drz;

            const double dr2 = drx2 + dry2 + drz2;

            const double mask = (dr2 <= cutoffsquare) ? 1. : 0.;

            const double invdr2 = 1. / dr2 * mask;
            const double lj2 = sigmasquare * invdr2;
            const double lj6 = lj2 * lj2 * lj2;
            const double lj12 = lj6 * lj6;
            const double lj12m6 = lj12 - lj6;
            const double fac = epsilon24 * (lj12 + lj12m6) * invdr2;

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
            if (calculateGlobals) {
              double virialx = drx * fx;
              double virialy = dry * fy;
              double virialz = drz * fz;
              double upot = (epsilon24 * lj12m6 + shift6) * mask;

              upotSum += upot;
              virialSumX += virialx;
              virialSumY += virialy;
              virialSumZ += virialz;

              if (_duplicatedCalculations) {
                // for non-newton3 the division is in the post-processing step.
                if (newton3) {
                  upot *= 0.5;
                  virialx *= 0.5;
                  virialy *= 0.5;
                  virialz *= 0.5;
                }
                if (inbox1) {
                  upotSum += upot;
                  virialSumX += virialx;
                  virialSumY += virialy;
                  virialSumZ += virialz;
                }
                // for non-newton3 the second particle will be considered in a separate calculation
                if (newton3 and autopas::utils::inBox({xArr[j], yArr[j], zArr[j]}, _lowCorner, _highCorner)) {
                  upotSum += upot;
                  virialSumX += virialx;
                  virialSumY += virialy;
                  virialSumZ += virialz;
                }
              } else {
                // for non-newton3 we divide by 2 only in the postprocess step!
                upotSum += upot;
                virialSumX += virialx;
                virialSumY += virialy;
                virialSumZ += virialz;
              }
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
        if (calculateGlobals) {
          double virialx = drx * fx;
          double virialy = dry * fy;
          double virialz = drz * fz;
          double upot = (epsilon24 * lj12m6 + shift6);

          upotSum += upot;
          virialSumX += virialx;
          virialSumY += virialy;
          virialSumZ += virialz;

          if (_duplicatedCalculations) {
            // for non-newton3 the division is in the post-processing step.
            if (newton3) {
              upot *= 0.5;
              virialx *= 0.5;
              virialy *= 0.5;
              virialz *= 0.5;
            }
            if (inbox1) {
              upotSum += upot;
              virialSumX += virialx;
              virialSumY += virialy;
              virialSumZ += virialz;
            }
            // for non-newton3 the second particle will be considered in a separate calculation
            if (newton3 and autopas::utils::inBox({xptr[j], xptr[j], xptr[j]}, _lowCorner, _highCorner)) {
              upotSum += upot;
              virialSumX += virialx;
              virialSumY += virialy;
              virialSumZ += virialz;
            }
          } else {
            // for non-newton3 we divide by 2 only in the postprocess step!
            upotSum += upot;
            virialSumX += virialx;
            virialSumY += virialy;
            virialSumZ += virialz;
          }
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }

    if (calculateGlobals) {
      const int threadnum = autopas_get_thread_num();

      _aosThreadData[threadnum].upotSum += upotSum;
      _aosThreadData[threadnum].virialSum[0] += virialSumX;
      _aosThreadData[threadnum].virialSum[1] += virialSumY;
      _aosThreadData[threadnum].virialSum[2] += virialSumZ;
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
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
    }
  }

  /**
   * postprocesses global values, e.g. upot and virial
   * @param newton3
   */
  void postProcessGlobalValues(bool newton3) {
    if (_postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, please don't call postProcessGlobalValues(bool newton3) twice without calling "
          "resetGlobalValues().");
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
      throw utils::ExceptionHandler::AutoPasException(
          "Not yet postprocessed, please call postProcessGlobalValues first.");
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
      throw utils::ExceptionHandler::AutoPasException(
          "Not yet postprocessed, please call postProcessGlobalValues first.");
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

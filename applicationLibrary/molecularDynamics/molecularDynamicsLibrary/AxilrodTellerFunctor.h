/**
 * @file AxilrodTellerFunctor.h
 * @author M. Muehlhaeusser
 * @date 25/07/23
 */

#pragma once

#include <simde/x86/avx512.h>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool relevantForTuning = true>
class AxilrodTellerFunctor
    : public autopas::TriwiseFunctor<
          Particle, AxilrodTellerFunctor<Particle, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  AxilrodTellerFunctor() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit AxilrodTellerFunctor(double cutoff, void * /*dummy*/)
      : autopas::TriwiseFunctor<
            Particle, AxilrodTellerFunctor<Particle, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _cutoffSquaredPd{simde_mm512_set1_pd(cutoff * cutoff)},
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
  }

 public:
  /**
   * Constructor for Functor with mixing disabled. When using this functor it is necessary to call
   * setParticleProperties() to set internal constants because it does not use a particle properties library.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit AxilrodTellerFunctor(double cutoff) : AxilrodTellerFunctor(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like nu.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit AxilrodTellerFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : AxilrodTellerFunctor(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  /**
   * Returns name of functor. Intended for use with the iteration logger, to differentiate between calls to
   * iterateTriwise using different functors in the logs.
   * @return name of functor.
   */
  virtual std::string getName() { return "AxilrodTellerFunctorAutoVec"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle &i, Particle &j, Particle &k, bool newton3) {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy() or k.isDummy()) {
      return;
    }
    auto nu = _nu;
    if constexpr (useMixing) {
      nu = _PPLibrary->getMixingNu(i.getTypeId(), j.getTypeId(), k.getTypeId());
    }
    auto drij = j.getR() - i.getR();
    auto drjk = k.getR() - j.getR();
    auto drki = i.getR() - k.getR();

    double dr2ij = autopas::utils::ArrayMath::dot(drij, drij);
    double dr2jk = autopas::utils::ArrayMath::dot(drjk, drjk);
    double dr2ki = autopas::utils::ArrayMath::dot(drki, drki);

    // Check cutoff
    if (dr2ij > _cutoffSquared or dr2jk > _cutoffSquared or dr2ki > _cutoffSquared) {
      return;
    }

    // Dot products of distances belonging to one particle
    double dr2i = autopas::utils::ArrayMath::dot(drij, drki);
    double dr2j = autopas::utils::ArrayMath::dot(drij, drjk);
    double dr2k = autopas::utils::ArrayMath::dot(drjk, drki);

    double dr2ijk = dr2i * dr2j * dr2k;

    double dr2 = dr2ij * dr2jk * dr2ki;
    double dr5 = dr2 * dr2 * std::sqrt(dr2);
    double invdr5 = nu / dr5;

    auto fi = drjk * dr2i * (dr2j - dr2k) + drij * (dr2j * dr2k - dr2jk * dr2ki + 5.0 * dr2ijk / dr2ij) +
              drki * (-dr2j * dr2k + dr2ij * dr2jk - 5.0 * dr2ijk / dr2ki);
    fi *= 3.0 * invdr5;
    i.addF(fi);

    auto fj = fi;
    auto fk = fi;
    if (newton3) {
      fj = drki * dr2j * (dr2k - dr2i) + drij * (-dr2i * dr2k + dr2jk * dr2ki - 5.0 * dr2ijk / dr2ij) +
           drjk * (dr2i * dr2k - dr2ij * dr2ki + 5.0 * dr2ijk / dr2jk);
      fj *= 3.0 * invdr5;
      j.addF(fj);

      /* auto fk = drij * dr2k * (dr2i - dr2j)
                + drjk * (- dr2i * dr2j + dr2ij * dr2ki - 5.0 * dr2ijk / dr2jk)
                + drki * (dr2i * dr2j - dr2ij * dr2jk + 5.0 * dr2ijk / dr2ki);
      fk *= 3.0 * invdr5; */
      fk = (fi + fj) * (-1.0);
      k.addF(fk);
    }

    if (calculateGlobals) {
      // Virial is calculated as f_i * r_i
      auto virialI = fi * i.getR();
      // Calculate third of total potential energy from 3-body interaction
      double potentialEnergy = invdr5 * (dr2 - 3.0 * dr2ijk) / 3.0;

      const int threadnum = autopas::autopas_get_thread_num();
      if (i.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialI;
      }
      // for non-newton3 particles j and/or k will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        auto virialJ = fj * j.getR();
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialJ;
      }
      if (newton3 and k.isOwned()) {
        auto virialK = fk * k.getR();
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialK;
      }
    }
  }

  // TODO simde

  inline double simde_mm512_reduce_add_pd(__m512d a) {
#ifdef __AVX512F__
    _mm512_reduce_add_pd(a);
#else
    std::array<double, vecLength> a_vec = {0.0};
    simde_mm512_storeu_pd(a_vec.data(), a);
    double acc = 0;
    for (int i = 0; i < vecLength; ++i) {
      acc += a_vec[i];
    }
    return acc;
#endif
  }

  inline void simde_mm512_mask_i64scatter_pd(void *base_addr, __mmask8 k, __m512i vindex, __m512d a, const int scale) {
#ifdef __AVX512F__
    _mm512_mask_i64scatter_pd(base_addr, vindex, a, scale);
#else
    for (int i = 0; i < vecLength; ++i) {
      if (k & (1 << i)) {
        static_cast<double *>(base_addr)[vindex[i]] = a[i];
      }
    }
#endif
  }

  inline void simde_mm512_i64scatter_pd(void *base_addr, __m512i vindex, __m512d a, const int scale) {
#ifdef __AVX512F__
    _mm512_i64scatter_pd(base_addr, vindex, a, scale);
#else
    for (int i = 0; i < vecLength; ++i) {
      static_cast<double *>(base_addr)[vindex[i]] = a[i];
    }
#endif
  }

  inline simde__m512i simde_mm512_alignr_epi64(simde__m512i a, simde__m512i b, const int imm8) {
#ifdef __AVX512F__
    _mm512_alignr_epi64(a, b, imm8);
#else
    std::array<int64_t, 2 *vecLength> buffer = {0};
    simde_mm512_storeu_epi64(buffer.data() + vecLength, a);
    simde_mm512_storeu_epi64(buffer.data(), b);
    return simde_mm512_loadu_epi64(buffer.data() + imm8);
#endif
  }

  inline bool SoAParticlesInCutoff(const double *const x1ptr, const double *const y1ptr, const double *const z1ptr,
                                   const double *const x2ptr, const double *const y2ptr, const double *const z2ptr,
                                   size_t index1, size_t index2) {
    // TODO vectorize this maybe
    const SoAFloatPrecision drxij = x2ptr[index2] - x1ptr[index1];
    const SoAFloatPrecision dryij = y2ptr[index2] - y1ptr[index1];
    const SoAFloatPrecision drzij = z2ptr[index2] - z1ptr[index1];

    const SoAFloatPrecision drxij2 = drxij * drxij;
    const SoAFloatPrecision dryij2 = dryij * dryij;
    const SoAFloatPrecision drzij2 = drzij * drzij;

    const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
    return drij2 <= _cutoffSquared;
  }

  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImplOld<true>(soa);
    } else {
      SoAFunctorSingleImplOld<false>(soa);
    }
  }

  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImplOld<true>(soa1, soa2);
    } else {
      SoAFunctorPairImplOld<false>(soa1, soa2);
    }
  }

  /**
   * Functor for structure of arrays (SoA)
   *
   * This functor calculates the forces
   * between all particles of soa1 and soa2 and soa3.
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
   * @param soa3 Third structure of arrays.
   * @param newton3 defines whether or whether not to use newton 3
   */
  void SoAFunctorTriple(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                        autopas::SoAView<SoAArraysType> soa3, const bool newton3) {
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorTriple() is not implemented.");
    if (newton3) {
      SoAFunctorTripleImplOld<true>(soa1, soa2, soa3);
    } else {
      SoAFunctorTripleImplOld<false>(soa1, soa2, soa3);
    }
  }

 private:
  template <bool newton3j, bool newton3k, bool remainderIsMasked>
  void SoAKernel(const std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> &indicesK, const size_t kStart,
                 const simde__m512d &x1, const simde__m512d &y1, const simde__m512d &z1, const simde__m512d &x2,
                 const simde__m512d &y2, const simde__m512d &z2, const double *const __restrict x3ptr,
                 const double *const __restrict y3ptr, const double *const __restrict z3ptr, const size_t type1,
                 const size_t type2, const size_t *const __restrict type3ptr, const simde__m512d &drxij,
                 const simde__m512d &dryij, const simde__m512d &drzij, const simde__m512d &drij2, simde__m512d &fxiacc,
                 simde__m512d &fyiacc, simde__m512d &fziacc, simde__m512d &fxjacc, simde__m512d &fyjacc,
                 simde__m512d &fzjacc, double *const __restrict fx3ptr, double *const __restrict fy3ptr,
                 double *const __restrict fz3ptr, unsigned int rest = 0) {
    const simde__m512i vindex = remainderIsMasked ? simde_mm512_maskz_loadu_epi64(_masks[rest], &indicesK[kStart])
                                                  : simde_mm512_load_si512(&indicesK[kStart]);

    const simde__m512d x3 = remainderIsMasked ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], vindex, x3ptr, 8)
                                              : simde_mm512_i64gather_pd(vindex, x3ptr, 8);
    const simde__m512d y3 = remainderIsMasked ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], vindex, y3ptr, 8)
                                              : simde_mm512_i64gather_pd(vindex, y3ptr, 8);
    const simde__m512d z3 = remainderIsMasked ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], vindex, z3ptr, 8)
                                              : simde_mm512_i64gather_pd(vindex, z3ptr, 8);

    const simde__m512d drxjk = simde_mm512_sub_pd(x3, x2);
    const simde__m512d dryjk = simde_mm512_sub_pd(y3, y2);
    const simde__m512d drzjk = simde_mm512_sub_pd(z3, z2);

    const simde__m512d drxjk2 = simde_mm512_mul_pd(drxjk, drxjk);
    /*const simde__m512d dryjk2 = simde_mm512_mul_pd(dryjk, dryjk);
    const simde__m512d drzjk2 = simde_mm512_mul_pd(drzjk, drzjk);
    const simde__m512d drjk2PART = simde_mm512_add_pd(drxjk2, dryjk2);
    const simde__m512d drjk2 = simde_mm512_add_pd(drjk2PART, drzjk2);*/
    const simde__m512d drjk2PART = simde_mm512_fmadd_pd(dryjk, dryjk, drxjk2);
    const simde__m512d drjk2 = simde_mm512_fmadd_pd(drzjk, drzjk, drjk2PART);

    const simde__m512d drxki = simde_mm512_sub_pd(x1, x3);
    const simde__m512d dryki = simde_mm512_sub_pd(y1, y3);
    const simde__m512d drzki = simde_mm512_sub_pd(z1, z3);

    const simde__m512d drxki2 = simde_mm512_mul_pd(drxki, drxki);
    /*const simde__m512d dryki2 = simde_mm512_mul_pd(dryki, dryki);
    const simde__m512d drzki2 = simde_mm512_mul_pd(drzki, drzki);
    const simde__m512d drki2PART = simde_mm512_add_pd(drxki2, dryki2);
    const simde__m512d drki2 = simde_mm512_add_pd(drki2PART, drzki2);*/
    const simde__m512d drki2PART = simde_mm512_fmadd_pd(dryki, dryki, drxki2);
    const simde__m512d drki2 = simde_mm512_fmadd_pd(drzki, drzki, drki2PART);

    const simde__m512d drxi2 = simde_mm512_mul_pd(drxij, drxki);
    /*const simde__m512d dryi2 = simde_mm512_mul_pd(dryij, dryki);
    const simde__m512d drzi2 = simde_mm512_mul_pd(drzij, drzki);
    const simde__m512d dri2PART = simde_mm512_add_pd(drxi2, dryi2);
    const simde__m512d dri2 = simde_mm512_add_pd(dri2PART, drzi2);*/
    const simde__m512d dri2PART = simde_mm512_fmadd_pd(dryij, dryki, drxi2);
    const simde__m512d dri2 = simde_mm512_fmadd_pd(drzij, drzki, dri2PART);

    const simde__m512d drxj2 = simde_mm512_mul_pd(drxij, drxjk);
    /*const simde__m512d dryj2 = simde_mm512_mul_pd(dryij, dryjk);
    const simde__m512d drzj2 = simde_mm512_mul_pd(drzij, drzjk);
    const simde__m512d drj2PART = simde_mm512_add_pd(drxj2, dryj2);
    const simde__m512d drj2 = simde_mm512_add_pd(drj2PART, drzj2);*/
    const simde__m512d drj2PART = simde_mm512_fmadd_pd(dryij, dryjk, drxj2);
    const simde__m512d drj2 = simde_mm512_fmadd_pd(drzij, drzjk, drj2PART);

    const simde__m512d drxk2 = simde_mm512_mul_pd(drxjk, drxki);
    /*const simde__m512d dryk2 = simde_mm512_mul_pd(dryjk, dryki);
    const simde__m512d drzk2 = simde_mm512_mul_pd(drzjk, drzki);
    const simde__m512d drk2PART = simde_mm512_add_pd(drxk2, dryk2);
    const simde__m512d drk2 = simde_mm512_add_pd(drk2PART, drzk2);*/
    const simde__m512d drk2PART = simde_mm512_fmadd_pd(dryjk, dryki, drxk2);
    const simde__m512d drk2 = simde_mm512_fmadd_pd(drzjk, drzki, drk2PART);

    const simde__m512d drijk2PART = simde_mm512_mul_pd(dri2, drj2);
    const simde__m512d drijk2 = simde_mm512_mul_pd(drijk2PART, drk2);

    const simde__m512d dr2PART = simde_mm512_mul_pd(drij2, drjk2);
    const simde__m512d dr2 = simde_mm512_mul_pd(dr2PART, drki2);
    const simde__m512d dr2sqrt = simde_mm512_sqrt_pd(dr2);
    const simde__m512d dr5PART = simde_mm512_mul_pd(dr2, dr2);
    const simde__m512d dr5 = simde_mm512_mul_pd(dr5PART, dr2sqrt);

    if constexpr (useMixing) {
      _nuPd = simde_mm512_set_pd(
          not remainderIsMasked or rest > 7 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[kStart + 7]]) : 0,
          not remainderIsMasked or rest > 6 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[kStart + 6]]) : 0,
          not remainderIsMasked or rest > 5 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[kStart + 5]]) : 0,
          not remainderIsMasked or rest > 4 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[kStart + 4]]) : 0,
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[kStart + 3]]) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[kStart + 2]]) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[kStart + 1]]) : 0,
          _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[kStart + 0]]));
    }
    const simde__m512d invdr5 =
        remainderIsMasked ? simde_mm512_maskz_div_pd(_masks[rest], _nuPd, dr5) : simde_mm512_div_pd(_nuPd, dr5);
    const simde__m512d invdr53 = simde_mm512_mul_pd(_three, invdr5);

    const simde__m512d drj2_drk2 = simde_mm512_sub_pd(drj2, drk2);
    const simde__m512d i_jkFactor = simde_mm512_mul_pd(dri2, drj2_drk2);
    const simde__m512d drj2drk2 = simde_mm512_mul_pd(drj2, drk2);
    const simde__m512d drjk2drki2 = simde_mm512_mul_pd(drjk2, drki2);
    const simde__m512d drijk2BYdrij2 =
        remainderIsMasked ? simde_mm512_maskz_div_pd(_masks[rest], drijk2, drij2) : simde_mm512_div_pd(drijk2, drij2);
    const simde__m512d i_ijFactorPART = simde_mm512_sub_pd(drj2drk2, drjk2drki2);
    const simde__m512d i_ijFactor = simde_mm512_fmadd_pd(drijk2BYdrij2, _five, i_ijFactorPART);
    const simde__m512d drij2drjk2 = simde_mm512_mul_pd(drij2, drjk2);
    const simde__m512d drijk2BYdrki2 =
        remainderIsMasked ? simde_mm512_maskz_div_pd(_masks[rest], drijk2, drki2) : simde_mm512_div_pd(drijk2, drki2);
    const simde__m512d i_kiFactorPART = simde_mm512_sub_pd(drj2drk2, drij2drjk2);
    const simde__m512d i_kiFactor = simde_mm512_fmadd_pd(drijk2BYdrki2, _five, i_kiFactorPART);

    const simde__m512d fxiPART = simde_mm512_mul_pd(drxjk, i_jkFactor);
    const simde__m512d fxiPART2 = simde_mm512_fmadd_pd(drxij, i_ijFactor, fxiPART);
    const simde__m512d fxiPART3 = simde_mm512_fnmadd_pd(drxki, i_kiFactor, fxiPART2);
    const simde__m512d fxi = simde_mm512_mul_pd(invdr53, fxiPART3);
    fxiacc = simde_mm512_add_pd(fxi, fxiacc);
    const simde__m512d fyiPART = simde_mm512_mul_pd(dryjk, i_jkFactor);
    const simde__m512d fyiPART2 = simde_mm512_fmadd_pd(dryij, i_ijFactor, fyiPART);
    const simde__m512d fyiPART3 = simde_mm512_fnmadd_pd(dryki, i_kiFactor, fyiPART2);
    const simde__m512d fyi = simde_mm512_mul_pd(invdr53, fyiPART3);
    fyiacc = simde_mm512_add_pd(fyi, fyiacc);
    const simde__m512d fziPART = simde_mm512_mul_pd(drzjk, i_jkFactor);
    const simde__m512d fziPART2 = simde_mm512_fmadd_pd(drzij, i_ijFactor, fziPART);
    const simde__m512d fziPART3 = simde_mm512_fnmadd_pd(drzki, i_kiFactor, fziPART2);
    const simde__m512d fzi = simde_mm512_mul_pd(invdr53, fziPART3);
    fziacc = simde_mm512_add_pd(fzi, fziacc);

    /*const auto fxiPART3 = drxjk * dri2 * (drj2 - drk2) +
                     drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) -
                     drxki * (drj2 * drk2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
    fxiacc += fxiPART3 * 3.0 * invdr5;
    const auto fyiPART3 = dryjk * dri2 * (drj2 - drk2) +
                     dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) -
                     dryki * (drj2 * drk2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
    fyiacc += fyiPART3 * 3.0 * invdr5;
    const auto fziPART3 = drzjk * dri2 * (drj2 - drk2) +
                     drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) -
                     drzki * (drj2 * drk2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
    fziacc += fziPART3 * 3.0 * invdr5;*/

    if (newton3j) {
      const simde__m512d drk2_dri2 = simde_mm512_sub_pd(drk2, dri2);
      const simde__m512d j_kiFactor = simde_mm512_mul_pd(drj2, drk2_dri2);
      const simde__m512d dri2drk2 = simde_mm512_mul_pd(dri2, drk2);
      const simde__m512d j_ijFactorPART = simde_mm512_sub_pd(dri2drk2, drjk2drki2);
      const simde__m512d j_ijFactor = simde_mm512_fmadd_pd(drijk2BYdrij2, _five, j_ijFactorPART);
      const simde__m512d drij2drki2 = simde_mm512_mul_pd(drij2, drki2);
      const simde__m512d drijk2BYdrjk2 =
          remainderIsMasked ? simde_mm512_maskz_div_pd(_masks[rest], drijk2, drjk2) : simde_mm512_div_pd(drijk2, drjk2);
      const simde__m512d j_jkFactorPART = simde_mm512_sub_pd(dri2drk2, drij2drki2);
      const simde__m512d j_jkFactor = simde_mm512_fmadd_pd(drijk2BYdrjk2, _five, j_jkFactorPART);

      const simde__m512d fxjPART = simde_mm512_mul_pd(drxki, j_kiFactor);
      const simde__m512d fxjPART2 = simde_mm512_fnmadd_pd(drxij, j_ijFactor, fxjPART);
      const simde__m512d fxjPART3 = simde_mm512_fmadd_pd(drxjk, j_jkFactor, fxjPART2);
      const simde__m512d fxj = simde_mm512_mul_pd(invdr53, fxjPART3);
      fxjacc = simde_mm512_add_pd(fxj, fxjacc);
      const simde__m512d fyjPART = simde_mm512_mul_pd(dryki, j_kiFactor);
      const simde__m512d fyjPART2 = simde_mm512_fnmadd_pd(dryij, j_ijFactor, fyjPART);
      const simde__m512d fyjPART3 = simde_mm512_fmadd_pd(dryjk, j_jkFactor, fyjPART2);
      const simde__m512d fyj = simde_mm512_mul_pd(invdr53, fyjPART3);
      fyjacc = simde_mm512_add_pd(fyj, fyjacc);
      const simde__m512d fzjPART = simde_mm512_mul_pd(drzki, j_kiFactor);
      const simde__m512d fzjPART2 = simde_mm512_fnmadd_pd(drzij, j_ijFactor, fzjPART);
      const simde__m512d fzjPART3 = simde_mm512_fmadd_pd(drzjk, j_jkFactor, fzjPART2);
      const simde__m512d fzj = simde_mm512_mul_pd(invdr53, fzjPART3);
      fzjacc = simde_mm512_add_pd(fzj, fzjacc);

      /*const auto fxjPART3 = drxki * drj2 * (drk2 - dri2) -
                       drxij * (dri2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                       drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
      fxjacc += fxjPART3 * 3.0 * invdr5;
      const auto fyjPART3 = dryki * drj2 * (drk2 - dri2) -
                       dryij * (dri2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                       dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
      fyjacc += fyjPART3 * 3.0 * invdr5;
      const auto fzjPART3 = drzki * drj2 * (drk2 - dri2) -
                       drzij * (dri2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                       drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
      fzjacc += fzjPART3 * 3.0 * invdr5;*/

      if constexpr (newton3k) {
        const simde__m512d fxk = simde_mm512_add_pd(fxi, fxj);
        const simde__m512d fyk = simde_mm512_add_pd(fyi, fyj);
        const simde__m512d fzk = simde_mm512_add_pd(fzi, fzj);

        const simde__m512d fxk_old = remainderIsMasked
                                         ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], vindex, fx3ptr, 8)
                                         : simde_mm512_i64gather_pd(vindex, fx3ptr, 8);
        const simde__m512d fyk_old = remainderIsMasked
                                         ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], vindex, fy3ptr, 8)
                                         : simde_mm512_i64gather_pd(vindex, fy3ptr, 8);
        const simde__m512d fzk_old = remainderIsMasked
                                         ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], vindex, fz3ptr, 8)
                                         : simde_mm512_i64gather_pd(vindex, fz3ptr, 8);

        const simde__m512d fxk_new = simde_mm512_sub_pd(fxk_old, fxk);
        const simde__m512d fyk_new = simde_mm512_sub_pd(fyk_old, fyk);
        const simde__m512d fzk_new = simde_mm512_sub_pd(fzk_old, fzk);

#ifndef __AVX512F__
        std::array<double, vecLength> buffer_fxk = {0.0};
        std::array<double, vecLength> buffer_fyk = {0.0};
        std::array<double, vecLength> buffer_fzk = {0.0};
#endif
        if constexpr (remainderIsMasked) {
#ifdef __AVX512F__
          simde_mm512_mask_i64scatter_pd(fx3ptr, _masks[rest], vindex, fxk_new, 8);
          simde_mm512_mask_i64scatter_pd(fy3ptr, _masks[rest], vindex, fyk_new, 8);
          simde_mm512_mask_i64scatter_pd(fy3ptr, _masks[rest], vindex, fzk_new, 8);
#endif
          simde_mm512_mask_storeu_pd(buffer_fxk.data(), _masks[rest], fxk_new);
          simde_mm512_mask_storeu_pd(buffer_fyk.data(), _masks[rest], fyk_new);
          simde_mm512_mask_storeu_pd(buffer_fzk.data(), _masks[rest], fzk_new);
        } else {
#ifdef __AVX512F__
          simde_mm512_i64scatter_pd(fx3ptr, vindex, fxk_new, 8);
          simde_mm512_i64scatter_pd(fy3ptr, vindex, fyk_new, 8);
          simde_mm512_i64scatter_pd(fz3ptr, vindex, fzk_new, 8);
#endif
          simde_mm512_storeu_pd(buffer_fxk.data(), fxk_new);
          simde_mm512_storeu_pd(buffer_fyk.data(), fyk_new);
          simde_mm512_storeu_pd(buffer_fzk.data(), fzk_new);
        }
#ifndef __AVX512F__
        for (int index = 0; index < (remainderIsMasked ? rest : vecLength); ++index) {
          fx3ptr[indicesK[kStart + index]] = buffer_fxk[index];
          fy3ptr[indicesK[kStart + index]] = buffer_fyk[index];
          fz3ptr[indicesK[kStart + index]] = buffer_fzk[index];
        }
#endif
      }
    }
    if constexpr (calculateGlobals) {
      autopas::utils::ExceptionHandler::exception("Globals not yet implemented.");  // TODO
    }
  }

  template <bool newton3j, bool newton3k, bool remainderIsMasked>
  void SoAKernel(const simde__m512i indicesK, const simde__m512d &x1, const simde__m512d &y1, const simde__m512d &z1,
                 const simde__m512d &x2, const simde__m512d &y2, const simde__m512d &z2,
                 const double *const __restrict x3ptr, const double *const __restrict y3ptr,
                 const double *const __restrict z3ptr, const size_t type1, const size_t type2,
                 const size_t *const __restrict type3ptr, const simde__m512d &drxij, const simde__m512d &dryij,
                 const simde__m512d &drzij, const simde__m512d &drij2, simde__m512d &fxiacc, simde__m512d &fyiacc,
                 simde__m512d &fziacc, simde__m512d &fxjacc, simde__m512d &fyjacc, simde__m512d &fzjacc,
                 double *const __restrict fx3ptr, double *const __restrict fy3ptr, double *const __restrict fz3ptr,
                 unsigned int rest = 0) {
    const simde__m512d x3 = remainderIsMasked ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, x3ptr, 8)
                                              : simde_mm512_i64gather_pd(indicesK, x3ptr, 8);
    const simde__m512d y3 = remainderIsMasked ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, y3ptr, 8)
                                              : simde_mm512_i64gather_pd(indicesK, y3ptr, 8);
    const simde__m512d z3 = remainderIsMasked ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, z3ptr, 8)
                                              : simde_mm512_i64gather_pd(indicesK, z3ptr, 8);

    const simde__m512d drxjk = simde_mm512_sub_pd(x3, x2);
    const simde__m512d dryjk = simde_mm512_sub_pd(y3, y2);
    const simde__m512d drzjk = simde_mm512_sub_pd(z3, z2);

    const simde__m512d drxjk2 = simde_mm512_mul_pd(drxjk, drxjk);
    /*const simde__m512d dryjk2 = simde_mm512_mul_pd(dryjk, dryjk);
    const simde__m512d drzjk2 = simde_mm512_mul_pd(drzjk, drzjk);
    const simde__m512d drjk2PART = simde_mm512_add_pd(drxjk2, dryjk2);
    const simde__m512d drjk2 = simde_mm512_add_pd(drjk2PART, drzjk2);*/
    const simde__m512d drjk2PART = simde_mm512_fmadd_pd(dryjk, dryjk, drxjk2);
    const simde__m512d drjk2 = simde_mm512_fmadd_pd(drzjk, drzjk, drjk2PART);

    const simde__m512d drxki = simde_mm512_sub_pd(x1, x3);
    const simde__m512d dryki = simde_mm512_sub_pd(y1, y3);
    const simde__m512d drzki = simde_mm512_sub_pd(z1, z3);

    const simde__m512d drxki2 = simde_mm512_mul_pd(drxki, drxki);
    /*const simde__m512d dryki2 = simde_mm512_mul_pd(dryki, dryki);
    const simde__m512d drzki2 = simde_mm512_mul_pd(drzki, drzki);
    const simde__m512d drki2PART = simde_mm512_add_pd(drxki2, dryki2);
    const simde__m512d drki2 = simde_mm512_add_pd(drki2PART, drzki2);*/
    const simde__m512d drki2PART = simde_mm512_fmadd_pd(dryki, dryki, drxki2);
    const simde__m512d drki2 = simde_mm512_fmadd_pd(drzki, drzki, drki2PART);

    const simde__m512d drxi2 = simde_mm512_mul_pd(drxij, drxki);
    /*const simde__m512d dryi2 = simde_mm512_mul_pd(dryij, dryki);
    const simde__m512d drzi2 = simde_mm512_mul_pd(drzij, drzki);
    const simde__m512d dri2PART = simde_mm512_add_pd(drxi2, dryi2);
    const simde__m512d dri2 = simde_mm512_add_pd(dri2PART, drzi2);*/
    const simde__m512d dri2PART = simde_mm512_fmadd_pd(dryij, dryki, drxi2);
    const simde__m512d dri2 = simde_mm512_fmadd_pd(drzij, drzki, dri2PART);

    const simde__m512d drxj2 = simde_mm512_mul_pd(drxij, drxjk);
    /*const simde__m512d dryj2 = simde_mm512_mul_pd(dryij, dryjk);
    const simde__m512d drzj2 = simde_mm512_mul_pd(drzij, drzjk);
    const simde__m512d drj2PART = simde_mm512_add_pd(drxj2, dryj2);
    const simde__m512d drj2 = simde_mm512_add_pd(drj2PART, drzj2);*/
    const simde__m512d drj2PART = simde_mm512_fmadd_pd(dryij, dryjk, drxj2);
    const simde__m512d drj2 = simde_mm512_fmadd_pd(drzij, drzjk, drj2PART);

    const simde__m512d drxk2 = simde_mm512_mul_pd(drxjk, drxki);
    /*const simde__m512d dryk2 = simde_mm512_mul_pd(dryjk, dryki);
    const simde__m512d drzk2 = simde_mm512_mul_pd(drzjk, drzki);
    const simde__m512d drk2PART = simde_mm512_add_pd(drxk2, dryk2);
    const simde__m512d drk2 = simde_mm512_add_pd(drk2PART, drzk2);*/
    const simde__m512d drk2PART = simde_mm512_fmadd_pd(dryjk, dryki, drxk2);
    const simde__m512d drk2 = simde_mm512_fmadd_pd(drzjk, drzki, drk2PART);

    const simde__m512d drijk2PART = simde_mm512_mul_pd(dri2, drj2);
    const simde__m512d drijk2 = simde_mm512_mul_pd(drijk2PART, drk2);

    const simde__m512d dr2PART = simde_mm512_mul_pd(drij2, drjk2);
    const simde__m512d dr2 = simde_mm512_mul_pd(dr2PART, drki2);
    const simde__m512d dr2sqrt = simde_mm512_sqrt_pd(dr2);
    const simde__m512d dr5PART = simde_mm512_mul_pd(dr2, dr2);
    const simde__m512d dr5 = simde_mm512_mul_pd(dr5PART, dr2sqrt);

    if constexpr (useMixing) {
      _nuPd = simde_mm512_set_pd(
          not remainderIsMasked or rest > 7 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[7]]) : 0,
          not remainderIsMasked or rest > 6 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[6]]) : 0,
          not remainderIsMasked or rest > 5 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[5]]) : 0,
          not remainderIsMasked or rest > 4 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[4]]) : 0,
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[3]]) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[2]]) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[1]]) : 0,
          _PPLibrary->getMixingNu(type1, type2, type3ptr[indicesK[0]]));
    }
    const simde__m512d invdr5 =
        remainderIsMasked ? simde_mm512_maskz_div_pd(_masks[rest], _nuPd, dr5) : simde_mm512_div_pd(_nuPd, dr5);
    const simde__m512d invdr53 = simde_mm512_mul_pd(_three, invdr5);

    const simde__m512d drj2_drk2 = simde_mm512_sub_pd(drj2, drk2);
    const simde__m512d i_jkFactor = simde_mm512_mul_pd(dri2, drj2_drk2);
    const simde__m512d drj2drk2 = simde_mm512_mul_pd(drj2, drk2);
    const simde__m512d drjk2drki2 = simde_mm512_mul_pd(drjk2, drki2);
    const simde__m512d drijk2BYdrij2 =
        remainderIsMasked ? simde_mm512_maskz_div_pd(_masks[rest], drijk2, drij2) : simde_mm512_div_pd(drijk2, drij2);
    const simde__m512d i_ijFactorPART = simde_mm512_sub_pd(drj2drk2, drjk2drki2);
    const simde__m512d i_ijFactor = simde_mm512_fmadd_pd(drijk2BYdrij2, _five, i_ijFactorPART);
    const simde__m512d drij2drjk2 = simde_mm512_mul_pd(drij2, drjk2);
    const simde__m512d drijk2BYdrki2 =
        remainderIsMasked ? simde_mm512_maskz_div_pd(_masks[rest], drijk2, drki2) : simde_mm512_div_pd(drijk2, drki2);
    const simde__m512d i_kiFactorPART = simde_mm512_sub_pd(drj2drk2, drij2drjk2);
    const simde__m512d i_kiFactor = simde_mm512_fmadd_pd(drijk2BYdrki2, _five, i_kiFactorPART);

    const simde__m512d fxiPART = simde_mm512_mul_pd(drxjk, i_jkFactor);
    const simde__m512d fxiPART2 = simde_mm512_fmadd_pd(drxij, i_ijFactor, fxiPART);
    const simde__m512d fxiPART3 = simde_mm512_fnmadd_pd(drxki, i_kiFactor, fxiPART2);
    const simde__m512d fxi = simde_mm512_mul_pd(invdr53, fxiPART3);
    fxiacc = simde_mm512_add_pd(fxi, fxiacc);
    const simde__m512d fyiPART = simde_mm512_mul_pd(dryjk, i_jkFactor);
    const simde__m512d fyiPART2 = simde_mm512_fmadd_pd(dryij, i_ijFactor, fyiPART);
    const simde__m512d fyiPART3 = simde_mm512_fnmadd_pd(dryki, i_kiFactor, fyiPART2);
    const simde__m512d fyi = simde_mm512_mul_pd(invdr53, fyiPART3);
    fyiacc = simde_mm512_add_pd(fyi, fyiacc);
    const simde__m512d fziPART = simde_mm512_mul_pd(drzjk, i_jkFactor);
    const simde__m512d fziPART2 = simde_mm512_fmadd_pd(drzij, i_ijFactor, fziPART);
    const simde__m512d fziPART3 = simde_mm512_fnmadd_pd(drzki, i_kiFactor, fziPART2);
    const simde__m512d fzi = simde_mm512_mul_pd(invdr53, fziPART3);
    fziacc = simde_mm512_add_pd(fzi, fziacc);

    /*const auto fxiPART3 = drxjk * dri2 * (drj2 - drk2) +
                     drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) -
                     drxki * (drj2 * drk2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
    fxiacc += fxiPART3 * 3.0 * invdr5;
    const auto fyiPART3 = dryjk * dri2 * (drj2 - drk2) +
                     dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) -
                     dryki * (drj2 * drk2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
    fyiacc += fyiPART3 * 3.0 * invdr5;
    const auto fziPART3 = drzjk * dri2 * (drj2 - drk2) +
                     drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) -
                     drzki * (drj2 * drk2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
    fziacc += fziPART3 * 3.0 * invdr5;*/

    if constexpr (newton3j) {
      const simde__m512d drk2_dri2 = simde_mm512_sub_pd(drk2, dri2);
      const simde__m512d j_kiFactor = simde_mm512_mul_pd(drj2, drk2_dri2);
      const simde__m512d dri2drk2 = simde_mm512_mul_pd(dri2, drk2);
      const simde__m512d j_ijFactorPART = simde_mm512_sub_pd(dri2drk2, drjk2drki2);
      const simde__m512d j_ijFactor = simde_mm512_fmadd_pd(drijk2BYdrij2, _five, j_ijFactorPART);
      const simde__m512d drij2drki2 = simde_mm512_mul_pd(drij2, drki2);
      const simde__m512d drijk2BYdrjk2 =
          remainderIsMasked ? simde_mm512_maskz_div_pd(_masks[rest], drijk2, drjk2) : simde_mm512_div_pd(drijk2, drjk2);
      const simde__m512d j_jkFactorPART = simde_mm512_sub_pd(dri2drk2, drij2drki2);
      const simde__m512d j_jkFactor = simde_mm512_fmadd_pd(drijk2BYdrjk2, _five, j_jkFactorPART);

      const simde__m512d fxjPART = simde_mm512_mul_pd(drxki, j_kiFactor);
      const simde__m512d fxjPART2 = simde_mm512_fnmadd_pd(drxij, j_ijFactor, fxjPART);
      const simde__m512d fxjPART3 = simde_mm512_fmadd_pd(drxjk, j_jkFactor, fxjPART2);
      const simde__m512d fxj = simde_mm512_mul_pd(invdr53, fxjPART3);
      fxjacc = simde_mm512_add_pd(fxj, fxjacc);
      const simde__m512d fyjPART = simde_mm512_mul_pd(dryki, j_kiFactor);
      const simde__m512d fyjPART2 = simde_mm512_fnmadd_pd(dryij, j_ijFactor, fyjPART);
      const simde__m512d fyjPART3 = simde_mm512_fmadd_pd(dryjk, j_jkFactor, fyjPART2);
      const simde__m512d fyj = simde_mm512_mul_pd(invdr53, fyjPART3);
      fyjacc = simde_mm512_add_pd(fyj, fyjacc);
      const simde__m512d fzjPART = simde_mm512_mul_pd(drzki, j_kiFactor);
      const simde__m512d fzjPART2 = simde_mm512_fnmadd_pd(drzij, j_ijFactor, fzjPART);
      const simde__m512d fzjPART3 = simde_mm512_fmadd_pd(drzjk, j_jkFactor, fzjPART2);
      const simde__m512d fzj = simde_mm512_mul_pd(invdr53, fzjPART3);
      fzjacc = simde_mm512_add_pd(fzj, fzjacc);

      /*const auto fxjPART3 = drxki * drj2 * (drk2 - dri2) -
                       drxij * (dri2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                       drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
      fxjacc += fxjPART3 * 3.0 * invdr5;
      const auto fyjPART3 = dryki * drj2 * (drk2 - dri2) -
                       dryij * (dri2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                       dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
      fyjacc += fyjPART3 * 3.0 * invdr5;
      const auto fzjPART3 = drzki * drj2 * (drk2 - dri2) -
                       drzij * (dri2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                       drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
      fzjacc += fzjPART3 * 3.0 * invdr5;*/

      if constexpr (newton3k) {
        const simde__m512d fxk = simde_mm512_add_pd(fxi, fxj);
        const simde__m512d fyk = simde_mm512_add_pd(fyi, fyj);
        const simde__m512d fzk = simde_mm512_add_pd(fzi, fzj);

        const simde__m512d fxk_old = remainderIsMasked
                                         ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, fx3ptr, 8)
                                         : simde_mm512_i64gather_pd(indicesK, fx3ptr, 8);
        const simde__m512d fyk_old = remainderIsMasked
                                         ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, fy3ptr, 8)
                                         : simde_mm512_i64gather_pd(indicesK, fy3ptr, 8);
        const simde__m512d fzk_old = remainderIsMasked
                                         ? simde_mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, fz3ptr, 8)
                                         : simde_mm512_i64gather_pd(indicesK, fz3ptr, 8);

        const simde__m512d fxk_new = simde_mm512_sub_pd(fxk_old, fxk);
        const simde__m512d fyk_new = simde_mm512_sub_pd(fyk_old, fyk);
        const simde__m512d fzk_new = simde_mm512_sub_pd(fzk_old, fzk);

#ifndef __AVX512F__
        std::array<double, vecLength> buffer_fxk = {0.0};
        std::array<double, vecLength> buffer_fyk = {0.0};
        std::array<double, vecLength> buffer_fzk = {0.0};
#endif
        if constexpr (remainderIsMasked) {
          // TODO scatter
#ifdef __AVX512F__
          simde_mm512_mask_i64scatter_pd(fx3ptr, _masks[rest], indicesK, fxk_new, 8);
          simde_mm512_mask_i64scatter_pd(fy3ptr, _masks[rest], indicesK, fyk_new, 8);
          simde_mm512_mask_i64scatter_pd(fy3ptr, _masks[rest], indicesK, fzk_new, 8);
#endif
          simde_mm512_mask_storeu_pd(buffer_fxk.data(), _masks[rest], fxk_new);
          simde_mm512_mask_storeu_pd(buffer_fyk.data(), _masks[rest], fyk_new);
          simde_mm512_mask_storeu_pd(buffer_fzk.data(), _masks[rest], fzk_new);
        } else {
          // TODO scatter
#ifdef __AVX512F__
          simde_mm512_i64scatter_pd(fx3ptr, indicesK, fxk_new, 8);
          simde_mm512_i64scatter_pd(fy3ptr, indicesK, fyk_new, 8);
          simde_mm512_i64scatter_pd(fz3ptr, indicesK, fzk_new, 8);
#endif
          simde_mm512_storeu_pd(buffer_fxk.data(), fxk_new);
          simde_mm512_storeu_pd(buffer_fyk.data(), fyk_new);
          simde_mm512_storeu_pd(buffer_fzk.data(), fzk_new);
        }
#ifndef __AVX512F__
        for (int index = 0; index < (remainderIsMasked ? rest : vecLength); ++index) {  // TODO mask_store
          fx3ptr[indicesK[index]] = buffer_fxk[index];
          fy3ptr[indicesK[index]] = buffer_fyk[index];
          fz3ptr[indicesK[index]] = buffer_fzk[index];
        }
#endif
      }
    }
  }

  template <bool newton3j, bool newton3k, bool remainderIsMasked>
  void SoACompressAlignr(simde__m512i &interactionIndices, int &numAssignedRegisters, const size_t k,
                         const simde__m512d x1, const simde__m512d y1, const simde__m512d z1, const simde__m512d x2,
                         const simde__m512d y2, const simde__m512d z2, const double *const __restrict x3ptr,
                         const double *const __restrict y3ptr, const double *const __restrict z3ptr,
                         const size_t *const ownedState3ptr, const size_t type1, const size_t type2,
                         const size_t *const __restrict type3ptr, const simde__m512d &drxij, const simde__m512d &dryij,
                         const simde__m512d &drzij, const simde__m512d &drij2, simde__m512d &fxiacc,
                         simde__m512d &fyiacc, simde__m512d &fziacc, simde__m512d &fxjacc, simde__m512d &fyjacc,
                         simde__m512d &fzjacc, double *const __restrict fx3ptr, double *const __restrict fy3ptr,
                         double *const __restrict fz3ptr, unsigned int rest = 0) {
    const simde__m512i loopIndices = simde_mm512_add_epi64(simde_mm512_set1_epi64(k), _ascendingIndices);

    const simde__m512d x3 =
        remainderIsMasked ? simde_mm512_maskz_loadu_pd(_masks[rest], &x3ptr[k]) : simde_mm512_loadu_pd(&x3ptr[k]);
    const simde__m512d y3 =
        remainderIsMasked ? simde_mm512_maskz_loadu_pd(_masks[rest], &y3ptr[k]) : simde_mm512_loadu_pd(&y3ptr[k]);
    const simde__m512d z3 =
        remainderIsMasked ? simde_mm512_maskz_loadu_pd(_masks[rest], &z3ptr[k]) : simde_mm512_loadu_pd(&z3ptr[k]);

    const simde__m512d drxki = simde_mm512_sub_pd(x1, x3);
    const simde__m512d dryki = simde_mm512_sub_pd(y1, y3);
    const simde__m512d drzki = simde_mm512_sub_pd(z1, z3);

    const simde__m512d drxki2 = simde_mm512_mul_pd(drxki, drxki);
    const simde__m512d dryki2 = simde_mm512_mul_pd(dryki, dryki);
    const simde__m512d drzki2 = simde_mm512_mul_pd(drzki, drzki);

    const simde__m512d drki2PART = simde_mm512_add_pd(drxki2, dryki2);
    const simde__m512d drki2 = simde_mm512_add_pd(drki2PART, drzki2);

    const simde__mmask8 cutoffMask_ki = simde_mm512_cmp_pd_mask(drki2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

    const simde__m512d drxjk = simde_mm512_sub_pd(x3, x2);
    const simde__m512d dryjk = simde_mm512_sub_pd(y3, y2);
    const simde__m512d drzjk = simde_mm512_sub_pd(z3, z2);

    const simde__m512d drxjk2 = simde_mm512_mul_pd(drxjk, drxjk);
    const simde__m512d dryjk2 = simde_mm512_mul_pd(dryjk, dryjk);
    const simde__m512d drzjk2 = simde_mm512_mul_pd(drzjk, drzjk);

    const simde__m512d drjk2PART = simde_mm512_add_pd(drxjk2, dryjk2);
    const simde__m512d drjk2 = simde_mm512_add_pd(drjk2PART, drzjk2);

    const simde__mmask8 cutoffMask_jk = simde_mm512_cmp_pd_mask(drjk2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

    const simde__mmask8 cutoffMask = simde_mm512_kand(cutoffMask_jk, cutoffMask_ki);

    const simde__m512i ownershipState3 = simde_mm512_loadu_epi64(&ownedState3ptr[k]);
    const simde__mmask8 dummyMask =
        simde_mm512_cmp_epi64_mask(ownershipState3, _ownedStateDummyEpi64, SIMDE_MM_CMPINT_NE);

    const simde__mmask8 mask = remainderIsMasked
                                   ? simde_mm512_kand(_masks[rest], simde_mm512_kand(cutoffMask, dummyMask))
                                   : simde_mm512_kand(cutoffMask, dummyMask);

    const int popCountMask = std::bitset<8>(mask).count();

    const simde__m512i newInteractionIndices = simde_mm512_maskz_compress_epi64(mask, loopIndices);
    if (numAssignedRegisters + popCountMask < 8) {
      interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, popCountMask);
      numAssignedRegisters += popCountMask;
    } else {
      interactionIndices =
          simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, 8 - numAssignedRegisters);

      SoAKernel<newton3j, newton3k, false>(interactionIndices, x1, y1, z1, x2, y2, z2, x3ptr, y3ptr, z3ptr, type1,
                                           type2, type3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc,
                                           fyjacc, fzjacc, fx3ptr, fy3ptr, fz3ptr);

      interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, _zeroI, popCountMask);
      numAssignedRegisters += popCountMask - 8;
    }
  }

  /**
   * Implementation function of SoAFunctorSingle(soa, newton3)
   *
   * @tparam newton3
   * @param soa
   */
  template <bool newton3>
  void SoAFunctorSingleImplOld(autopas::SoAView<SoAArraysType> soa) {
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorSingle() is not implemented.");
    if (soa.size() == 0) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesK;
      indicesK.reserve(soa.size());

    if constexpr (not useMixing) {
      _nuPd = simde_mm512_set1_pd(_nu);
    }

    for (size_t i = 0; i < soa.size(); ++i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      const simde__m512d x1 = simde_mm512_set1_pd(xptr[i]);
      const simde__m512d y1 = simde_mm512_set1_pd(yptr[i]);
      const simde__m512d z1 = simde_mm512_set1_pd(zptr[i]);

      simde__m512d fxiacc = simde_mm512_setzero_pd();
      simde__m512d fyiacc = simde_mm512_setzero_pd();
      simde__m512d fziacc = simde_mm512_setzero_pd();

      // std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesJ;
      for (size_t j = i + 1; j < soa.size(); ++j) {
        if (ownedStatePtr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(xptr, yptr, zptr, xptr, yptr, zptr, i, j)) {
          continue;
        }

        const simde__m512d x2 = simde_mm512_set1_pd(xptr[j]);
        const simde__m512d y2 = simde_mm512_set1_pd(yptr[j]);
        const simde__m512d z2 = simde_mm512_set1_pd(zptr[j]);

        simde__m512d fxjacc = simde_mm512_setzero_pd();
        simde__m512d fyjacc = simde_mm512_setzero_pd();
        simde__m512d fzjacc = simde_mm512_setzero_pd();

        const simde__m512d drxij = simde_mm512_sub_pd(x2, x1);
        const simde__m512d dryij = simde_mm512_sub_pd(y2, y1);
        const simde__m512d drzij = simde_mm512_sub_pd(z2, z1);

        const simde__m512d drxij2 = simde_mm512_mul_pd(drxij, drxij);
        const simde__m512d dryij2 = simde_mm512_mul_pd(dryij, dryij);
        const simde__m512d drzij2 = simde_mm512_mul_pd(drzij, drzij);

        const simde__m512d drij2PART = simde_mm512_add_pd(drxij2, dryij2);
        const simde__m512d drij2 = simde_mm512_add_pd(drij2PART, drzij2);

        indicesK.clear();
        for (size_t k = j + 1; k < soa.size(); ++k) {
          if (ownedStatePtr[k] == autopas::OwnershipState::dummy or
              not SoAParticlesInCutoff(xptr, yptr, zptr, xptr, yptr, zptr, j, k) or
              not SoAParticlesInCutoff(xptr, yptr, zptr, xptr, yptr, zptr, k, i)) {
            continue;
          }

          indicesK.push_back(k);
        }

        size_t kIndex = 0;
        for (; kIndex < (indicesK.size() & ~(vecLength - 1)); kIndex += vecLength) {
          SoAKernel<true, true, false>(indicesK, kIndex, x1, y1, z1, x2, y2, z2, xptr, yptr, zptr, typeptr[i],
                                       typeptr[j], typeptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc,
                                       fyjacc, fzjacc, fxptr, fyptr, fzptr);
        }
        if (kIndex < indicesK.size()) {
          SoAKernel<true, true, true>(indicesK, kIndex, x1, y1, z1, x2, y2, z2, xptr, yptr, zptr, typeptr[i],
                                      typeptr[j], typeptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc,
                                      fyjacc, fzjacc, fxptr, fyptr, fzptr, indicesK.size() - kIndex);
        }
        fxptr[j] += simde_mm512_reduce_add_pd(fxjacc);
        fyptr[j] += simde_mm512_reduce_add_pd(fyjacc);
        fzptr[j] += simde_mm512_reduce_add_pd(fzjacc);
      }
      fxptr[i] += simde_mm512_reduce_add_pd(fxiacc);
      fyptr[i] += simde_mm512_reduce_add_pd(fyiacc);
      fzptr[i] += simde_mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      autopas::utils::ExceptionHandler::exception("Globals not yet implemented.");  // TODO
    }
  }

  template <bool newton3, bool compress = false>
  void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorSingle() is not implemented.");
    if (soa.size() == 0) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesK;
    if (not compress) {
      indicesK.reserve(soa.size());
    }

    if constexpr (not useMixing) {
      _nuPd = simde_mm512_set1_pd(_nu);
    }

    for (size_t i = soa.size(); i > 0; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      const simde__m512d x1 = simde_mm512_set1_pd(xptr[i]);
      const simde__m512d y1 = simde_mm512_set1_pd(yptr[i]);
      const simde__m512d z1 = simde_mm512_set1_pd(zptr[i]);

      simde__m512d fxiacc = simde_mm512_setzero_pd();
      simde__m512d fyiacc = simde_mm512_setzero_pd();
      simde__m512d fziacc = simde_mm512_setzero_pd();

      // std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesJ;
      for (size_t j = i - 1; j > 0; --j) {
        if (ownedStatePtr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(xptr, yptr, zptr, xptr, yptr, zptr, i, j)) {
          continue;
        }

        const simde__m512d x2 = simde_mm512_set1_pd(xptr[j]);
        const simde__m512d y2 = simde_mm512_set1_pd(yptr[j]);
        const simde__m512d z2 = simde_mm512_set1_pd(zptr[j]);

        simde__m512d fxjacc = simde_mm512_setzero_pd();
        simde__m512d fyjacc = simde_mm512_setzero_pd();
        simde__m512d fzjacc = simde_mm512_setzero_pd();

        const simde__m512d drxij = simde_mm512_sub_pd(x2, x1);
        const simde__m512d dryij = simde_mm512_sub_pd(y2, y1);
        const simde__m512d drzij = simde_mm512_sub_pd(z2, z1);

        const simde__m512d drxij2 = simde_mm512_mul_pd(drxij, drxij);
        const simde__m512d dryij2 = simde_mm512_mul_pd(dryij, dryij);
        const simde__m512d drzij2 = simde_mm512_mul_pd(drzij, drzij);

        const simde__m512d drij2PART = simde_mm512_add_pd(drxij2, dryij2);
        const simde__m512d drij2 = simde_mm512_add_pd(drij2PART, drzij2);
//
        fxptr[j] += simde_mm512_reduce_add_pd(fxjacc);
        fyptr[j] += simde_mm512_reduce_add_pd(fyjacc);
        fzptr[j] += simde_mm512_reduce_add_pd(fzjacc);
      }
      fxptr[i] += simde_mm512_reduce_add_pd(fxiacc);
      fyptr[i] += simde_mm512_reduce_add_pd(fyiacc);
      fzptr[i] += simde_mm512_reduce_add_pd(fziacc);
    }
  }

  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImplOld(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesK;
    indicesK.reserve(soa2.size());

    if constexpr (not useMixing) {
      _nuPd = simde_mm512_set1_pd(_nu);
    }

    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      const simde__m512d x1 = simde_mm512_set1_pd(x1ptr[i]);
      const simde__m512d y1 = simde_mm512_set1_pd(y1ptr[i]);
      const simde__m512d z1 = simde_mm512_set1_pd(z1ptr[i]);

      simde__m512d fxiacc = simde_mm512_setzero_pd();
      simde__m512d fyiacc = simde_mm512_setzero_pd();
      simde__m512d fziacc = simde_mm512_setzero_pd();

      // TODO maybe change loop order to always pick particle 2 from soa2 and soa3 from either soa1 or soa2
      // particle 2 from soa1 and 3 from soa2

      // std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesJ;
      for (size_t j = i + 1; j < soa1.size(); ++j) {
        if (ownedState1ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x1ptr, y1ptr, z1ptr, i, j)) {
          continue;
        }

        const simde__m512d x2 = simde_mm512_set1_pd(x1ptr[j]);
        const simde__m512d y2 = simde_mm512_set1_pd(y1ptr[j]);
        const simde__m512d z2 = simde_mm512_set1_pd(z1ptr[j]);

        simde__m512d fxjacc = simde_mm512_setzero_pd();
        simde__m512d fyjacc = simde_mm512_setzero_pd();
        simde__m512d fzjacc = simde_mm512_setzero_pd();

        const simde__m512d drxij = simde_mm512_sub_pd(x2, x1);
        const simde__m512d dryij = simde_mm512_sub_pd(y2, y1);
        const simde__m512d drzij = simde_mm512_sub_pd(z2, z1);

        const simde__m512d drxij2 = simde_mm512_mul_pd(drxij, drxij);
        const simde__m512d dryij2 = simde_mm512_mul_pd(dryij, dryij);
        const simde__m512d drzij2 = simde_mm512_mul_pd(drzij, drzij);

        const simde__m512d drij2PART = simde_mm512_add_pd(drxij2, dryij2);
        const simde__m512d drij2 = simde_mm512_add_pd(drij2PART, drzij2);

        indicesK.clear();
        for (size_t k = 0; k < soa2.size(); ++k) {
          if (ownedState2ptr[k] == autopas::OwnershipState::dummy or
              not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, j, k) or
              not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, k)) {
            continue;
          }

          indicesK.push_back(k);
        }

        size_t kIndex = 0;
        for (; kIndex < (indicesK.size() & ~(vecLength - 1)); kIndex += vecLength) {
          SoAKernel<true, newton3, false>(indicesK, kIndex, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr, type1ptr[i],
                                          type1ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                          fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr);
        }
        if (kIndex < indicesK.size()) {
          SoAKernel<true, newton3, true>(indicesK, kIndex, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr, type1ptr[i],
                                         type1ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                         fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr, indicesK.size() - kIndex);
        }
        fx1ptr[j] += simde_mm512_reduce_add_pd(fxjacc);
        fy1ptr[j] += simde_mm512_reduce_add_pd(fyjacc);
        fz1ptr[j] += simde_mm512_reduce_add_pd(fzjacc);
      }

      // both particles 2 and 3 from soa2

      for (size_t j = 0; j < soa2.size(); ++j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, j)) {
          continue;
        }

        const simde__m512d x2 = simde_mm512_set1_pd(x2ptr[j]);
        const simde__m512d y2 = simde_mm512_set1_pd(y2ptr[j]);
        const simde__m512d z2 = simde_mm512_set1_pd(z2ptr[j]);

        simde__m512d fxjacc = simde_mm512_setzero_pd();
        simde__m512d fyjacc = simde_mm512_setzero_pd();
        simde__m512d fzjacc = simde_mm512_setzero_pd();

        const simde__m512d drxij = simde_mm512_sub_pd(x2, x1);
        const simde__m512d dryij = simde_mm512_sub_pd(y2, y1);
        const simde__m512d drzij = simde_mm512_sub_pd(z2, z1);

        const simde__m512d drxij2 = simde_mm512_mul_pd(drxij, drxij);
        const simde__m512d dryij2 = simde_mm512_mul_pd(dryij, dryij);
        const simde__m512d drzij2 = simde_mm512_mul_pd(drzij, drzij);

        const simde__m512d drij2PART = simde_mm512_add_pd(drxij2, dryij2);
        const simde__m512d drij2 = simde_mm512_add_pd(drij2PART, drzij2);

        // particle 3 from soa 2
        indicesK.clear();
        for (size_t k = j + 1; k < soa2.size(); ++k) {
          if (ownedState2ptr[k] == autopas::OwnershipState::dummy or
              not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, k) or
              not SoAParticlesInCutoff(x2ptr, y2ptr, z2ptr, x2ptr, y2ptr, z2ptr, j, k)) {
            continue;
          }

          indicesK.push_back(k);
        }

        size_t kIndex = 0;
        for (; kIndex < (indicesK.size() & ~(vecLength - 1)); kIndex += vecLength) {
          SoAKernel<newton3, newton3, false>(indicesK, kIndex, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr, type1ptr[i],
                                             type2ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                             fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr);
        }
        if (kIndex < indicesK.size()) {
          SoAKernel<newton3, newton3, true>(indicesK, kIndex, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr, type1ptr[i],
                                            type2ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                            fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr, indicesK.size() - kIndex);
        }

        if constexpr (newton3) {
          fx2ptr[j] += simde_mm512_reduce_add_pd(fxjacc);
          fy2ptr[j] += simde_mm512_reduce_add_pd(fyjacc);
          fz2ptr[j] += simde_mm512_reduce_add_pd(fzjacc);
        }
      }

      fx1ptr[i] += simde_mm512_reduce_add_pd(fxiacc);
      fy1ptr[i] += simde_mm512_reduce_add_pd(fyiacc);
      fz1ptr[i] += simde_mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      autopas::utils::ExceptionHandler::exception("Globals not yet implemented.");  // TODO
    }
  }

  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesK;
    indicesK.reserve(soa2.size());

    if constexpr (not useMixing) {
      _nuPd = simde_mm512_set1_pd(_nu);
    }

    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      const simde__m512d x1 = simde_mm512_set1_pd(x1ptr[i]);
      const simde__m512d y1 = simde_mm512_set1_pd(y1ptr[i]);
      const simde__m512d z1 = simde_mm512_set1_pd(z1ptr[i]);

      simde__m512d fxiacc = simde_mm512_setzero_pd();
      simde__m512d fyiacc = simde_mm512_setzero_pd();
      simde__m512d fziacc = simde_mm512_setzero_pd();

      // TODO maybe change loop order to always pick particle 2 from soa2 and soa3 from either soa1 or soa2
      // particle 2 from soa1 and 3 from soa2

      // std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesJ;
      for (size_t j = i + 1; j < soa1.size(); ++j) {
        if (ownedState1ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x1ptr, y1ptr, z1ptr, i, j)) {
          continue;
        }

        const simde__m512d x2 = simde_mm512_set1_pd(x1ptr[j]);
        const simde__m512d y2 = simde_mm512_set1_pd(y1ptr[j]);
        const simde__m512d z2 = simde_mm512_set1_pd(z1ptr[j]);

        simde__m512d fxjacc = simde_mm512_setzero_pd();
        simde__m512d fyjacc = simde_mm512_setzero_pd();
        simde__m512d fzjacc = simde_mm512_setzero_pd();

        const simde__m512d drxij = simde_mm512_sub_pd(x2, x1);
        const simde__m512d dryij = simde_mm512_sub_pd(y2, y1);
        const simde__m512d drzij = simde_mm512_sub_pd(z2, z1);

        const simde__m512d drxij2 = simde_mm512_mul_pd(drxij, drxij);
        const simde__m512d dryij2 = simde_mm512_mul_pd(dryij, dryij);
        const simde__m512d drzij2 = simde_mm512_mul_pd(drzij, drzij);

        const simde__m512d drij2PART = simde_mm512_add_pd(drxij2, dryij2);
        const simde__m512d drij2 = simde_mm512_add_pd(drij2PART, drzij2);

        simde__m512i interactionIndices = _zeroI;
        int numAssignedRegisters = 0;

        size_t k = 0;
        for (; k < (soa2.size() & ~(vecLength - 1)); k += vecLength) {
          const simde__m512i loopIndices = simde_mm512_add_epi64(simde_mm512_set1_epi64(k), _ascendingIndices);

          const simde__m512d x3 = simde_mm512_loadu_pd(&x2ptr[k]);
          const simde__m512d y3 = simde_mm512_loadu_pd(&y2ptr[k]);
          const simde__m512d z3 = simde_mm512_loadu_pd(&z2ptr[k]);

          const simde__m512d drxki = simde_mm512_sub_pd(x1, x3);
          const simde__m512d dryki = simde_mm512_sub_pd(y1, y3);
          const simde__m512d drzki = simde_mm512_sub_pd(z1, z3);

          const simde__m512d drxki2 = simde_mm512_mul_pd(drxki, drxki);
          const simde__m512d dryki2 = simde_mm512_mul_pd(dryki, dryki);
          const simde__m512d drzki2 = simde_mm512_mul_pd(drzki, drzki);

          const simde__m512d drki2PART = simde_mm512_add_pd(drxki2, dryki2);
          const simde__m512d drki2 = simde_mm512_add_pd(drki2PART, drzki2);

          const simde__mmask8 cutoffMask_ki = simde_mm512_cmp_pd_mask(drki2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

          const simde__m512d drxjk = simde_mm512_sub_pd(x3, x2);
          const simde__m512d dryjk = simde_mm512_sub_pd(y3, y2);
          const simde__m512d drzjk = simde_mm512_sub_pd(z3, z2);

          const simde__m512d drxjk2 = simde_mm512_mul_pd(drxjk, drxjk);
          const simde__m512d dryjk2 = simde_mm512_mul_pd(dryjk, dryjk);
          const simde__m512d drzjk2 = simde_mm512_mul_pd(drzjk, drzjk);

          const simde__m512d drjk2PART = simde_mm512_add_pd(drxjk2, dryjk2);
          const simde__m512d drjk2 = simde_mm512_add_pd(drjk2PART, drzjk2);

          const simde__mmask8 cutoffMask_jk = simde_mm512_cmp_pd_mask(drjk2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

          const simde__mmask8 cutoffMask = simde_mm512_kand(cutoffMask_jk, cutoffMask_ki);

          const simde__m512i ownershipState3 = simde_mm512_loadu_epi64(&ownedState2ptr[k]);
          const simde__mmask8 dummyMask =
              simde_mm512_cmp_epi64_mask(ownershipState3, _ownedStateDummyEpi64, SIMDE_MM_CMPINT_NE);

          const simde__mmask8 mask = simde_mm512_kand(cutoffMask, dummyMask);

          const int popCountMask = std::bitset<8>(mask).count();

          const simde__m512i newInteractionIndices = simde_mm512_maskz_compress_epi64(mask, loopIndices);
          if (numAssignedRegisters + popCountMask < 8) {
            interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, popCountMask);
            numAssignedRegisters += popCountMask;
          } else {
            interactionIndices =
                simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, 8 - numAssignedRegisters);

            SoAKernel<true, newton3, false>(interactionIndices, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr,
                                            type1ptr[i], type1ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc,
                                            fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr);

            interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, _zeroI, popCountMask);
            numAssignedRegisters += popCountMask - 8;
          }
        }
        if (k < soa2.size()) {
          const simde__mmask8 remainderMask = _masks[soa2.size() - k];
          const simde__m512i loopIndices = simde_mm512_add_epi64(simde_mm512_set1_epi64(k), _ascendingIndices);

          const simde__m512d x3 = simde_mm512_maskz_loadu_pd(remainderMask, &x2ptr[k]);
          const simde__m512d y3 = simde_mm512_maskz_loadu_pd(remainderMask, &y2ptr[k]);
          const simde__m512d z3 = simde_mm512_maskz_loadu_pd(remainderMask, &z2ptr[k]);

          const simde__m512d drxki = simde_mm512_sub_pd(x1, x3);
          const simde__m512d dryki = simde_mm512_sub_pd(y1, y3);
          const simde__m512d drzki = simde_mm512_sub_pd(z1, z3);

          const simde__m512d drxki2 = simde_mm512_mul_pd(drxki, drxki);
          const simde__m512d dryki2 = simde_mm512_mul_pd(dryki, dryki);
          const simde__m512d drzki2 = simde_mm512_mul_pd(drzki, drzki);

          const simde__m512d drki2PART = simde_mm512_add_pd(drxki2, dryki2);
          const simde__m512d drki2 = simde_mm512_add_pd(drki2PART, drzki2);

          const simde__mmask8 cutoffMask_ki = simde_mm512_cmp_pd_mask(drki2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

          const simde__m512d drxjk = simde_mm512_sub_pd(x3, x2);
          const simde__m512d dryjk = simde_mm512_sub_pd(y3, y2);
          const simde__m512d drzjk = simde_mm512_sub_pd(z3, z2);

          const simde__m512d drxjk2 = simde_mm512_mul_pd(drxjk, drxjk);
          const simde__m512d dryjk2 = simde_mm512_mul_pd(dryjk, dryjk);
          const simde__m512d drzjk2 = simde_mm512_mul_pd(drzjk, drzjk);

          const simde__m512d drjk2PART = simde_mm512_add_pd(drxjk2, dryjk2);
          const simde__m512d drjk2 = simde_mm512_add_pd(drjk2PART, drzjk2);

          const simde__mmask8 cutoffMask_jk = simde_mm512_cmp_pd_mask(drjk2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

          const simde__mmask8 cutoffMask = simde_mm512_kand(cutoffMask_jk, cutoffMask_ki);

          const simde__m512i ownershipState3 = simde_mm512_loadu_epi64(&ownedState2ptr[k]);
          const simde__mmask8 dummyMask =
              simde_mm512_cmp_epi64_mask(ownershipState3, _ownedStateDummyEpi64, SIMDE_MM_CMPINT_NE);

          const simde__mmask8 mask = simde_mm512_kand(remainderMask, simde_mm512_kand(cutoffMask, dummyMask));

          const int popCountMask = std::bitset<8>(mask).count();

          const simde__m512i newInteractionIndices = simde_mm512_maskz_compress_epi64(mask, loopIndices);
          if (numAssignedRegisters + popCountMask < 8) {
            interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, popCountMask);
            numAssignedRegisters += popCountMask;
          } else {
            interactionIndices =
                simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, 8 - numAssignedRegisters);

            SoAKernel<true, newton3, false>(interactionIndices, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr,
                                            type1ptr[i], type1ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc,
                                            fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr);

            interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, _zeroI, popCountMask);
            numAssignedRegisters += popCountMask - 8;
          }
        }
        if (numAssignedRegisters > 0) {
          // shift remainder to right
          interactionIndices = simde_mm512_alignr_epi64(_zeroI, interactionIndices, 8 - numAssignedRegisters);
          SoAKernel<true, newton3, true>(interactionIndices, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr, type1ptr[i],
                                         type1ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                         fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr, numAssignedRegisters);
        }
        fx1ptr[j] += simde_mm512_reduce_add_pd(fxjacc);
        fy1ptr[j] += simde_mm512_reduce_add_pd(fyjacc);
        fz1ptr[j] += simde_mm512_reduce_add_pd(fzjacc);
      }

      // both particles 2 and 3 from soa2

      for (size_t j = soa2.size() - 1; j > 0; --j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, j)) {
          continue;
        }

        const simde__m512d x2 = simde_mm512_set1_pd(x2ptr[j]);
        const simde__m512d y2 = simde_mm512_set1_pd(y2ptr[j]);
        const simde__m512d z2 = simde_mm512_set1_pd(z2ptr[j]);

        simde__m512d fxjacc = simde_mm512_setzero_pd();
        simde__m512d fyjacc = simde_mm512_setzero_pd();
        simde__m512d fzjacc = simde_mm512_setzero_pd();

        const simde__m512d drxij = simde_mm512_sub_pd(x2, x1);
        const simde__m512d dryij = simde_mm512_sub_pd(y2, y1);
        const simde__m512d drzij = simde_mm512_sub_pd(z2, z1);

        const simde__m512d drxij2 = simde_mm512_mul_pd(drxij, drxij);
        const simde__m512d dryij2 = simde_mm512_mul_pd(dryij, dryij);
        const simde__m512d drzij2 = simde_mm512_mul_pd(drzij, drzij);

        const simde__m512d drij2PART = simde_mm512_add_pd(drxij2, dryij2);
        const simde__m512d drij2 = simde_mm512_add_pd(drij2PART, drzij2);

        // particle 3 from soa 2

        simde__m512i interactionIndices = _zeroI;
        int numAssignedRegisters = 0;

        size_t k = 0;
        for (; k < (j & ~(vecLength - 1)); k += vecLength) {
          const simde__m512i loopIndices = simde_mm512_add_epi64(simde_mm512_set1_epi64(k), _ascendingIndices);

          const simde__m512d x3 = simde_mm512_loadu_pd(&x2ptr[k]);
          const simde__m512d y3 = simde_mm512_loadu_pd(&y2ptr[k]);
          const simde__m512d z3 = simde_mm512_loadu_pd(&z2ptr[k]);

          const simde__m512d drxki = simde_mm512_sub_pd(x1, x3);
          const simde__m512d dryki = simde_mm512_sub_pd(y1, y3);
          const simde__m512d drzki = simde_mm512_sub_pd(z1, z3);

          const simde__m512d drxki2 = simde_mm512_mul_pd(drxki, drxki);
          const simde__m512d dryki2 = simde_mm512_mul_pd(dryki, dryki);
          const simde__m512d drzki2 = simde_mm512_mul_pd(drzki, drzki);

          const simde__m512d drki2PART = simde_mm512_add_pd(drxki2, dryki2);
          const simde__m512d drki2 = simde_mm512_add_pd(drki2PART, drzki2);

          const simde__mmask8 cutoffMask_ki = simde_mm512_cmp_pd_mask(drki2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

          const simde__m512d drxjk = simde_mm512_sub_pd(x3, x2);
          const simde__m512d dryjk = simde_mm512_sub_pd(y3, y2);
          const simde__m512d drzjk = simde_mm512_sub_pd(z3, z2);

          const simde__m512d drxjk2 = simde_mm512_mul_pd(drxjk, drxjk);
          const simde__m512d dryjk2 = simde_mm512_mul_pd(dryjk, dryjk);
          const simde__m512d drzjk2 = simde_mm512_mul_pd(drzjk, drzjk);

          const simde__m512d drjk2PART = simde_mm512_add_pd(drxjk2, dryjk2);
          const simde__m512d drjk2 = simde_mm512_add_pd(drjk2PART, drzjk2);

          const simde__mmask8 cutoffMask_jk = simde_mm512_cmp_pd_mask(drjk2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

          const simde__mmask8 cutoffMask = simde_mm512_kand(cutoffMask_jk, cutoffMask_ki);

          const simde__m512i ownershipState3 = simde_mm512_loadu_epi64(&ownedState2ptr[k]);
          const simde__mmask8 dummyMask =
              simde_mm512_cmp_epi64_mask(ownershipState3, _ownedStateDummyEpi64, SIMDE_MM_CMPINT_NE);

          const simde__mmask8 mask = simde_mm512_kand(cutoffMask, dummyMask);

          const int popCountMask = std::bitset<8>(mask).count();

          const simde__m512i newInteractionIndices = simde_mm512_maskz_compress_epi64(mask, loopIndices);
          if (numAssignedRegisters + popCountMask < 8) {
            interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, popCountMask);
            numAssignedRegisters += popCountMask;
          } else {
            interactionIndices =
                simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, 8 - numAssignedRegisters);

            SoAKernel<newton3, newton3, false>(interactionIndices, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr,
                                               type1ptr[i], type2ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc,
                                               fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr);

            interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, _zeroI, popCountMask);
            numAssignedRegisters += popCountMask - 8;
          }
        }
        if (k < j) {
          const simde__mmask8 remainderMask = _masks[j - k];
          const simde__m512i loopIndices = simde_mm512_add_epi64(simde_mm512_set1_epi64(k), _ascendingIndices);

          const simde__m512d x3 = simde_mm512_maskz_loadu_pd(remainderMask, &x2ptr[k]);
          const simde__m512d y3 = simde_mm512_maskz_loadu_pd(remainderMask, &y2ptr[k]);
          const simde__m512d z3 = simde_mm512_maskz_loadu_pd(remainderMask, &z2ptr[k]);

          const simde__m512d drxki = simde_mm512_sub_pd(x1, x3);
          const simde__m512d dryki = simde_mm512_sub_pd(y1, y3);
          const simde__m512d drzki = simde_mm512_sub_pd(z1, z3);

          const simde__m512d drxki2 = simde_mm512_mul_pd(drxki, drxki);
          const simde__m512d dryki2 = simde_mm512_mul_pd(dryki, dryki);
          const simde__m512d drzki2 = simde_mm512_mul_pd(drzki, drzki);

          const simde__m512d drki2PART = simde_mm512_add_pd(drxki2, dryki2);
          const simde__m512d drki2 = simde_mm512_add_pd(drki2PART, drzki2);

          const simde__mmask8 cutoffMask_ki = simde_mm512_cmp_pd_mask(drki2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

          const simde__m512d drxjk = simde_mm512_sub_pd(x3, x2);
          const simde__m512d dryjk = simde_mm512_sub_pd(y3, y2);
          const simde__m512d drzjk = simde_mm512_sub_pd(z3, z2);

          const simde__m512d drxjk2 = simde_mm512_mul_pd(drxjk, drxjk);
          const simde__m512d dryjk2 = simde_mm512_mul_pd(dryjk, dryjk);
          const simde__m512d drzjk2 = simde_mm512_mul_pd(drzjk, drzjk);

          const simde__m512d drjk2PART = simde_mm512_add_pd(drxjk2, dryjk2);
          const simde__m512d drjk2 = simde_mm512_add_pd(drjk2PART, drzjk2);

          const simde__mmask8 cutoffMask_jk = simde_mm512_cmp_pd_mask(drjk2, _cutoffSquaredPd, SIMDE_CMP_LE_OS);

          const simde__mmask8 cutoffMask = simde_mm512_kand(cutoffMask_jk, cutoffMask_ki);

          const simde__m512i ownershipState3 = simde_mm512_loadu_epi64(&ownedState2ptr[k]);
          const simde__mmask8 dummyMask =
              simde_mm512_cmp_epi64_mask(ownershipState3, _ownedStateDummyEpi64, SIMDE_MM_CMPINT_NE);

          const simde__mmask8 mask = simde_mm512_kand(remainderMask, simde_mm512_kand(cutoffMask, dummyMask));

          const int popCountMask = std::bitset<8>(mask).count();

          const simde__m512i newInteractionIndices = simde_mm512_maskz_compress_epi64(mask, loopIndices);
          if (numAssignedRegisters + popCountMask < 8) {
            interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, popCountMask);
            numAssignedRegisters += popCountMask;
          } else {
            interactionIndices =
                simde_mm512_alignr_epi64(newInteractionIndices, interactionIndices, 8 - numAssignedRegisters);

            SoAKernel<newton3, newton3, false>(interactionIndices, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr,
                                               type1ptr[i], type2ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc,
                                               fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr);

            interactionIndices = simde_mm512_alignr_epi64(newInteractionIndices, _zeroI, popCountMask);
            numAssignedRegisters += popCountMask - 8;
          }
        }
        if (numAssignedRegisters > 0) {
          // shift remainder to right
          interactionIndices = simde_mm512_alignr_epi64(_zeroI, interactionIndices, 8 - numAssignedRegisters);
          SoAKernel<newton3, newton3, true>(interactionIndices, x1, y1, z1, x2, y2, z2, x2ptr, y2ptr, z2ptr,
                                            type1ptr[i], type2ptr[j], type2ptr, drxij, dryij, drzij, drij2, fxiacc,
                                            fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr,
                                            numAssignedRegisters);
        }

        if constexpr (newton3) {
          fx2ptr[j] += simde_mm512_reduce_add_pd(fxjacc);
          fy2ptr[j] += simde_mm512_reduce_add_pd(fyjacc);
          fz2ptr[j] += simde_mm512_reduce_add_pd(fzjacc);
        }
      }

      fx1ptr[i] += simde_mm512_reduce_add_pd(fxiacc);
      fy1ptr[i] += simde_mm512_reduce_add_pd(fyiacc);
      fz1ptr[i] += simde_mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      autopas::utils::ExceptionHandler::exception("Globals not yet implemented.");  // TODO
    }
  }

  /**
   * Implementation function of SoAFunctorTriple(soa1, soa2, soa3, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   * @param soa3
   */
  template <bool newton3>
  void SoAFunctorTripleImplOld(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                               autopas::SoAView<SoAArraysType> soa3) {
    if (soa1.size() == 0 or soa2.size() == 0 or soa3.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x3ptr = soa3.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y3ptr = soa3.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z3ptr = soa3.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState3ptr = soa3.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx3ptr = soa3.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy3ptr = soa3.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz3ptr = soa3.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type3ptr = soa3.template begin<Particle::AttributeNames::typeId>();

    std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesK;

    if constexpr (not useMixing) {
      _nuPd = simde_mm512_set1_pd(_nu);
    }

    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      const simde__m512d x1 = simde_mm512_set1_pd(x1ptr[i]);
      const simde__m512d y1 = simde_mm512_set1_pd(y1ptr[i]);
      const simde__m512d z1 = simde_mm512_set1_pd(z1ptr[i]);

      simde__m512d fxiacc = simde_mm512_setzero_pd();
      simde__m512d fyiacc = simde_mm512_setzero_pd();
      simde__m512d fziacc = simde_mm512_setzero_pd();

      // std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesJ;
      for (size_t j = 0; j < soa2.size(); ++j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, j)) {
          continue;
        }

        // TODO scalar operations and set1 ?
        const simde__m512d x2 = simde_mm512_set1_pd(x2ptr[j]);
        const simde__m512d y2 = simde_mm512_set1_pd(y2ptr[j]);
        const simde__m512d z2 = simde_mm512_set1_pd(z2ptr[j]);

        simde__m512d fxjacc = simde_mm512_setzero_pd();
        simde__m512d fyjacc = simde_mm512_setzero_pd();
        simde__m512d fzjacc = simde_mm512_setzero_pd();

        const simde__m512d drxij = simde_mm512_sub_pd(x2, x1);
        const simde__m512d dryij = simde_mm512_sub_pd(y2, y1);
        const simde__m512d drzij = simde_mm512_sub_pd(z2, z1);

        const simde__m512d drxij2 = simde_mm512_mul_pd(drxij, drxij);
        const simde__m512d dryij2 = simde_mm512_mul_pd(dryij, dryij);
        const simde__m512d drzij2 = simde_mm512_mul_pd(drzij, drzij);

        const simde__m512d drij2PART = simde_mm512_add_pd(drxij2, dryij2);
        const simde__m512d drij2 = simde_mm512_add_pd(drij2PART, drzij2);

        indicesK.clear();
        for (size_t k = 0; k < soa3.size(); ++k) {
          if (ownedState3ptr[k] == autopas::OwnershipState::dummy or
              not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x3ptr, y3ptr, z3ptr, i, k) or
              not SoAParticlesInCutoff(x2ptr, y2ptr, z2ptr, x3ptr, y3ptr, z3ptr, j, k)) {
            continue;
          }

          indicesK.push_back(k);
        }

        size_t kIndex = 0;
        for (; kIndex < (indicesK.size() & ~(vecLength - 1)); kIndex += vecLength) {
          SoAKernel<newton3, newton3, false>(indicesK, kIndex, x1, y1, z1, x2, y2, z2, x3ptr, y3ptr, z3ptr, type1ptr[i],
                                             type2ptr[j], type3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                             fxjacc, fyjacc, fzjacc, fx3ptr, fy3ptr, fz3ptr);
        }
        if (kIndex < indicesK.size()) {
          SoAKernel<newton3, newton3, true>(indicesK, kIndex, x1, y1, z1, x2, y2, z2, x3ptr, y3ptr, z3ptr, type1ptr[i],
                                            type2ptr[j], type3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                            fxjacc, fyjacc, fzjacc, fx3ptr, fy3ptr, fz3ptr, indicesK.size() - kIndex);
        }

        if constexpr (newton3) {
          fx2ptr[j] += simde_mm512_reduce_add_pd(fxjacc);
          fy2ptr[j] += simde_mm512_reduce_add_pd(fyjacc);
          fz2ptr[j] += simde_mm512_reduce_add_pd(fzjacc);
        }
      }

      fx1ptr[i] += simde_mm512_reduce_add_pd(fxiacc);
      fy1ptr[i] += simde_mm512_reduce_add_pd(fyiacc);
      fz1ptr[i] += simde_mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      autopas::utils::ExceptionHandler::exception("Globals not yet implemented.");  // TODO
    }
  }

  template <bool newton3>
  void SoAFunctorTripleImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                            autopas::SoAView<SoAArraysType> soa3) {
    if (soa1.size() == 0 or soa2.size() == 0 or soa3.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x3ptr = soa3.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y3ptr = soa3.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z3ptr = soa3.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState3ptr = soa3.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx3ptr = soa3.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy3ptr = soa3.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz3ptr = soa3.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type3ptr = soa3.template begin<Particle::AttributeNames::typeId>();

    if constexpr (not useMixing) {
      _nuPd = simde_mm512_set1_pd(_nu);
    }

    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      const simde__m512d x1 = simde_mm512_set1_pd(x1ptr[i]);
      const simde__m512d y1 = simde_mm512_set1_pd(y1ptr[i]);
      const simde__m512d z1 = simde_mm512_set1_pd(z1ptr[i]);

      simde__m512d fxiacc = simde_mm512_setzero_pd();
      simde__m512d fyiacc = simde_mm512_setzero_pd();
      simde__m512d fziacc = simde_mm512_setzero_pd();

      // std::vector<size_t, autopas::AlignedAllocator<size_t, 64>> indicesJ;
      for (size_t j = 0; j < soa2.size(); ++j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, j)) {
          continue;
        }

        // TODO scalar operations and set1 ?
        const simde__m512d x2 = simde_mm512_set1_pd(x2ptr[j]);
        const simde__m512d y2 = simde_mm512_set1_pd(y2ptr[j]);
        const simde__m512d z2 = simde_mm512_set1_pd(z2ptr[j]);

        simde__m512d fxjacc = simde_mm512_setzero_pd();
        simde__m512d fyjacc = simde_mm512_setzero_pd();
        simde__m512d fzjacc = simde_mm512_setzero_pd();

        const simde__m512d drxij = simde_mm512_sub_pd(x2, x1);
        const simde__m512d dryij = simde_mm512_sub_pd(y2, y1);
        const simde__m512d drzij = simde_mm512_sub_pd(z2, z1);

        const simde__m512d drxij2 = simde_mm512_mul_pd(drxij, drxij);
        const simde__m512d dryij2 = simde_mm512_mul_pd(dryij, dryij);
        const simde__m512d drzij2 = simde_mm512_mul_pd(drzij, drzij);

        const simde__m512d drij2PART = simde_mm512_add_pd(drxij2, dryij2);
        const simde__m512d drij2 = simde_mm512_add_pd(drij2PART, drzij2);

        simde__m512i interactionIndices = _zeroI;
        int numAssignedRegisters = 0;

        size_t k = 0;
        for (; k < (soa3.size() & ~(vecLength - 1)); k += vecLength) {
          SoACompressAlignr<newton3, newton3, false>(interactionIndices, numAssignedRegisters, k, x1, y1, z1, x2, y2,
                                                     z2, x3ptr, y3ptr, z3ptr, ownedState3ptr, type1ptr[i], type2ptr[j],
                                                     type3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                                     fxjacc, fyjacc, fzjacc, fx3ptr, fy3ptr, fz3ptr);
        }
        if (k < soa3.size()) {
          SoACompressAlignr<newton3, newton3, false>(interactionIndices, numAssignedRegisters, k, x1, y1, z1, x2, y2,
                                                     z2, x3ptr, y3ptr, z3ptr, ownedState3ptr, type1ptr[i], type2ptr[j],
                                                     type3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc,
                                                     fxjacc, fyjacc, fzjacc, fx3ptr, fy3ptr, fz3ptr, soa3.size() - k);
        }
        if (numAssignedRegisters > 0) {
          // shift remainder to right
          interactionIndices = simde_mm512_alignr_epi64(_zeroI, interactionIndices, 8 - numAssignedRegisters);
          SoAKernel<newton3, newton3, true>(interactionIndices, x1, y1, z1, x2, y2, z2, x3ptr, y3ptr, z3ptr,
                                            type1ptr[i], type2ptr[j], type3ptr, drxij, dryij, drzij, drij2, fxiacc,
                                            fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fx3ptr, fy3ptr, fz3ptr,
                                            numAssignedRegisters);
        }

        if constexpr (newton3) {
          fx2ptr[j] += simde_mm512_reduce_add_pd(fxjacc);
          fy2ptr[j] += simde_mm512_reduce_add_pd(fyjacc);
          fz2ptr[j] += simde_mm512_reduce_add_pd(fzjacc);
        }
      }

      fx1ptr[i] += simde_mm512_reduce_add_pd(fxiacc);
      fy1ptr[i] += simde_mm512_reduce_add_pd(fyiacc);
      fz1ptr[i] += simde_mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      autopas::utils::ExceptionHandler::exception("Globals not yet implemented.");  // TODO
    }
  }

  template <bool newton3>
  void SoAFunctorSingleImplScalar(autopas::SoAView<SoAArraysType> soa) {
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorSingle() is not implemented.");
    //  TODO
    if (soa.size() == 0) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    std::vector<unsigned int> indicesK;

    for (size_t i = 0; i < soa.size(); ++i) {
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      for (size_t j = i + 1; j < soa.size(); ++j) {
        if (ownedStatePtr[j] == autopas::OwnershipState::dummy) {
          continue;
        }

        const SoAFloatPrecision drxij = xptr[j] - xptr[i];
        const SoAFloatPrecision dryij = yptr[j] - yptr[i];
        const SoAFloatPrecision drzij = zptr[j] - zptr[i];

        const SoAFloatPrecision drxij2 = drxij * drxij;
        const SoAFloatPrecision dryij2 = dryij * dryij;
        const SoAFloatPrecision drzij2 = drzij * drzij;

        const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
        if (drij2 > _cutoffSquared) {
          continue;
        }

        indicesK.clear();
        for (size_t k = j + 1; k < soa.size(); ++k) {
          if (ownedStatePtr[k] == autopas::OwnershipState::dummy) {
            continue;
          }

          const SoAFloatPrecision drxjk = xptr[k] - xptr[j];
          const SoAFloatPrecision dryjk = yptr[k] - yptr[j];
          const SoAFloatPrecision drzjk = zptr[k] - zptr[j];

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          const SoAFloatPrecision drxki = xptr[i] - xptr[k];
          const SoAFloatPrecision dryki = yptr[i] - yptr[k];
          const SoAFloatPrecision drzki = zptr[i] - zptr[k];

          const SoAFloatPrecision drxki2 = drxki * drxki;
          const SoAFloatPrecision dryki2 = dryki * dryki;
          const SoAFloatPrecision drzki2 = drzki * drzki;

          const SoAFloatPrecision drki2 = drxki2 + dryki2 + drzki2;

          if (drjk2 > _cutoffSquared or drki2 > _cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision drxi2 = drxij * drxki;
          const SoAFloatPrecision dryi2 = dryij * dryki;
          const SoAFloatPrecision drzi2 = drzij * drzki;

          const SoAFloatPrecision dri2 = drxi2 + dryi2 + drzi2;

          const SoAFloatPrecision drxj2 = drxij * drxjk;
          const SoAFloatPrecision dryj2 = dryij * dryjk;
          const SoAFloatPrecision drzj2 = drzij * drzjk;

          const SoAFloatPrecision drj2 = drxj2 + dryj2 + drzj2;

          const SoAFloatPrecision drxk2 = drxjk * drxki;
          const SoAFloatPrecision dryk2 = dryjk * dryki;
          const SoAFloatPrecision drzk2 = drzjk * drzki;

          const SoAFloatPrecision drk2 = drxk2 + dryk2 + drzk2;

          const SoAFloatPrecision drijk2 = dri2 * drj2 * drk2;

          const SoAFloatPrecision dr2 = drij2 * drjk2 * drki2;
          const SoAFloatPrecision dr5 = dr2 * dr2 * std::sqrt(dr2);

          auto nu = _nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(typeptr[i], typeptr[j], typeptr[k]);
          }
          const SoAFloatPrecision invdr5 = nu / dr5;

          const auto fxi = drxjk * dri2 * (drj2 - drk2) + drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           drxki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fxacc += fxi * 3.0 * invdr5;
          const auto fyi = dryjk * dri2 * (drj2 - drk2) + dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           dryki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fyacc += fyi * 3.0 * invdr5;
          const auto fzi = drzjk * dri2 * (drj2 - drk2) + drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           drzki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fzacc += fzi * 3.0 * invdr5;

          const auto fxj = drxki * drj2 * (drk2 - dri2) +
                           drxij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                           drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
          const auto fyj = dryki * drj2 * (drk2 - dri2) +
                           dryij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                           dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
          const auto fzj = drzki * drj2 * (drk2 - dri2) +
                           drzij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                           drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);

          fxptr[j] += fxj * 3.0 * invdr5;
          fyptr[j] += fyj * 3.0 * invdr5;
          fzptr[j] += fzj * 3.0 * invdr5;

          /*const auto fxk = drxij * drk2 * (dri2 - drj2) +
                           drxjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                           drxki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
          const auto fyk = dryij * drk2 * (dri2 - drj2) +
                           dryjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                           dryki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
          const auto fzk = drzij * drk2 * (dri2 - drj2) +
                           drzjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                           drzki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);*/
          const auto fxk = (fxi + fxj) * (-1.0);
          const auto fyk = (fyi + fyj) * (-1.0);
          const auto fzk = (fzi + fzj) * (-1.0);

          fxptr[k] += fxk * 3.0 * invdr5;
          fyptr[k] += fyk * 3.0 * invdr5;
          fzptr[k] += fzk * 3.0 * invdr5;
        }
      }
      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
  }

  template <bool newton3>
  void SoAFunctorPairImplScalar(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    // TODO: should always calculate forces for all particles in soa1, even when newton3 == false
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorPair() is not implemented.");
    // TODO
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    std::vector<unsigned int> indicesK;

    for (size_t i = 0; i < soa1.size(); ++i) {
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      // particle 2 from soa1 and 3 from soa2

      for (size_t j = i + 1; j < soa1.size(); ++j) {
        if (ownedState1ptr[j] == autopas::OwnershipState::dummy) {
          continue;
        }

        const SoAFloatPrecision drxij = x1ptr[j] - x1ptr[i];
        const SoAFloatPrecision dryij = y1ptr[j] - y1ptr[i];
        const SoAFloatPrecision drzij = z1ptr[j] - z1ptr[i];

        const SoAFloatPrecision drxij2 = drxij * drxij;
        const SoAFloatPrecision dryij2 = dryij * dryij;
        const SoAFloatPrecision drzij2 = drzij * drzij;

        const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
        if (drij2 > _cutoffSquared) {
          continue;
        }

        indicesK.clear();
        for (size_t k = 0; k < soa2.size(); ++k) {
          if (ownedState2ptr[k] == autopas::OwnershipState::dummy) {
            continue;
          }

          const SoAFloatPrecision drxjk = x2ptr[k] - x1ptr[j];
          const SoAFloatPrecision dryjk = y2ptr[k] - y1ptr[j];
          const SoAFloatPrecision drzjk = z2ptr[k] - z1ptr[j];

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          const SoAFloatPrecision drxki = x1ptr[i] - x2ptr[k];
          const SoAFloatPrecision dryki = y1ptr[i] - y2ptr[k];
          const SoAFloatPrecision drzki = z1ptr[i] - z2ptr[k];

          const SoAFloatPrecision drxki2 = drxki * drxki;
          const SoAFloatPrecision dryki2 = dryki * dryki;
          const SoAFloatPrecision drzki2 = drzki * drzki;

          const SoAFloatPrecision drki2 = drxki2 + dryki2 + drzki2;

          if (drjk2 > _cutoffSquared or drki2 > _cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision drxi2 = drxij * drxki;
          const SoAFloatPrecision dryi2 = dryij * dryki;
          const SoAFloatPrecision drzi2 = drzij * drzki;

          const SoAFloatPrecision dri2 = drxi2 + dryi2 + drzi2;

          const SoAFloatPrecision drxj2 = drxij * drxjk;
          const SoAFloatPrecision dryj2 = dryij * dryjk;
          const SoAFloatPrecision drzj2 = drzij * drzjk;

          const SoAFloatPrecision drj2 = drxj2 + dryj2 + drzj2;

          const SoAFloatPrecision drxk2 = drxjk * drxki;
          const SoAFloatPrecision dryk2 = dryjk * dryki;
          const SoAFloatPrecision drzk2 = drzjk * drzki;

          const SoAFloatPrecision drk2 = drxk2 + dryk2 + drzk2;

          const SoAFloatPrecision drijk2 = dri2 * drj2 * drk2;

          const SoAFloatPrecision dr2 = drij2 * drjk2 * drki2;
          const SoAFloatPrecision dr5 = dr2 * dr2 * std::sqrt(dr2);

          auto nu = _nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(type1ptr[i], type1ptr[j], type2ptr[k]);
          }
          const SoAFloatPrecision invdr5 = nu / dr5;

          const auto fxi = drxjk * dri2 * (drj2 - drk2) + drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           drxki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fxacc += fxi * 3.0 * invdr5;
          const auto fyi = dryjk * dri2 * (drj2 - drk2) + dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           dryki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fyacc += fyi * 3.0 * invdr5;
          const auto fzi = drzjk * dri2 * (drj2 - drk2) + drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           drzki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fzacc += fzi * 3.0 * invdr5;

          const auto fxj = drxki * drj2 * (drk2 - dri2) +
                           drxij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                           drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
          const auto fyj = dryki * drj2 * (drk2 - dri2) +
                           dryij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                           dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
          const auto fzj = drzki * drj2 * (drk2 - dri2) +
                           drzij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                           drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);

          fx1ptr[j] += fxj * 3.0 * invdr5;
          fy1ptr[j] += fyj * 3.0 * invdr5;
          fz1ptr[j] += fzj * 3.0 * invdr5;

          if (newton3) {
            /*const auto fxk = drxij * drk2 * (dri2 - drj2) +
                             drxjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                             drxki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
            const auto fyk = dryij * drk2 * (dri2 - drj2) +
                             dryjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                             dryki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
            const auto fzk = drzij * drk2 * (dri2 - drj2) +
                             drzjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                             drzki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);*/
            const auto fxk = (fxi + fxj) * (-1.0);
            const auto fyk = (fyi + fyj) * (-1.0);
            const auto fzk = (fzi + fzj) * (-1.0);

            fx2ptr[k] += fxk * 3.0 * invdr5;
            fy2ptr[k] += fyk * 3.0 * invdr5;
            fz2ptr[k] += fzk * 3.0 * invdr5;
          }
        }
      }

      // both particles 2 and 3 from soa2

      for (size_t j = 0; j < soa2.size(); ++j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy) {
          continue;
        }

        const SoAFloatPrecision drxij = x2ptr[j] - x1ptr[i];
        const SoAFloatPrecision dryij = y2ptr[j] - y1ptr[i];
        const SoAFloatPrecision drzij = z2ptr[j] - z1ptr[i];

        const SoAFloatPrecision drxij2 = drxij * drxij;
        const SoAFloatPrecision dryij2 = dryij * dryij;
        const SoAFloatPrecision drzij2 = drzij * drzij;

        const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
        if (drij2 > _cutoffSquared) {
          continue;
        }
        // particle 3 from soa 2
        {
          indicesK.clear();
          for (size_t k = j + 1; k < soa2.size(); ++k) {
            if (ownedState2ptr[k] == autopas::OwnershipState::dummy) {
              continue;
            }

            const SoAFloatPrecision drxjk = x2ptr[k] - x2ptr[j];
            const SoAFloatPrecision dryjk = y2ptr[k] - y2ptr[j];
            const SoAFloatPrecision drzjk = z2ptr[k] - z2ptr[j];

            const SoAFloatPrecision drxjk2 = drxjk * drxjk;
            const SoAFloatPrecision dryjk2 = dryjk * dryjk;
            const SoAFloatPrecision drzjk2 = drzjk * drzjk;

            const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

            const SoAFloatPrecision drxki = x1ptr[i] - x2ptr[k];
            const SoAFloatPrecision dryki = y1ptr[i] - y2ptr[k];
            const SoAFloatPrecision drzki = z1ptr[i] - z2ptr[k];

            const SoAFloatPrecision drxki2 = drxki * drxki;
            const SoAFloatPrecision dryki2 = dryki * dryki;
            const SoAFloatPrecision drzki2 = drzki * drzki;

            const SoAFloatPrecision drki2 = drxki2 + dryki2 + drzki2;

            if (drjk2 > _cutoffSquared or drki2 > _cutoffSquared) {
              continue;
            }

            const SoAFloatPrecision drxi2 = drxij * drxki;
            const SoAFloatPrecision dryi2 = dryij * dryki;
            const SoAFloatPrecision drzi2 = drzij * drzki;

            const SoAFloatPrecision dri2 = drxi2 + dryi2 + drzi2;

            const SoAFloatPrecision drxj2 = drxij * drxjk;
            const SoAFloatPrecision dryj2 = dryij * dryjk;
            const SoAFloatPrecision drzj2 = drzij * drzjk;

            const SoAFloatPrecision drj2 = drxj2 + dryj2 + drzj2;

            const SoAFloatPrecision drxk2 = drxjk * drxki;
            const SoAFloatPrecision dryk2 = dryjk * dryki;
            const SoAFloatPrecision drzk2 = drzjk * drzki;

            const SoAFloatPrecision drk2 = drxk2 + dryk2 + drzk2;

            const SoAFloatPrecision drijk2 = dri2 * drj2 * drk2;

            const SoAFloatPrecision dr2 = drij2 * drjk2 * drki2;
            const SoAFloatPrecision dr5 = dr2 * dr2 * std::sqrt(dr2);

            auto nu = _nu;
            if constexpr (useMixing) {
              nu = _PPLibrary->getMixingNu(type1ptr[i], type2ptr[j], type2ptr[k]);
            }
            const SoAFloatPrecision invdr5 = nu / dr5;

            const auto fxi = drxjk * dri2 * (drj2 - drk2) +
                             drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                             drxki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
            fxacc += fxi * 3.0 * invdr5;
            const auto fyi = dryjk * dri2 * (drj2 - drk2) +
                             dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                             dryki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
            fyacc += fyi * 3.0 * invdr5;
            const auto fzi = drzjk * dri2 * (drj2 - drk2) +
                             drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                             drzki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
            fzacc += fzi * 3.0 * invdr5;

            if (newton3) {
              const auto fxj = drxki * drj2 * (drk2 - dri2) +
                               drxij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                               drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
              const auto fyj = dryki * drj2 * (drk2 - dri2) +
                               dryij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                               dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
              const auto fzj = drzki * drj2 * (drk2 - dri2) +
                               drzij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                               drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);

              fx2ptr[j] += fxj * 3.0 * invdr5;
              fy2ptr[j] += fyj * 3.0 * invdr5;
              fz2ptr[j] += fzj * 3.0 * invdr5;

              /*const auto fxk = drxij * drk2 * (dri2 - drj2) +
                               drxjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                               drxki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
              const auto fyk = dryij * drk2 * (dri2 - drj2) +
                               dryjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                               dryki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
              const auto fzk = drzij * drk2 * (dri2 - drj2) +
                               drzjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                               drzki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);*/
              const auto fxk = (fxi + fxj) * (-1.0);
              const auto fyk = (fyi + fyj) * (-1.0);
              const auto fzk = (fzi + fzj) * (-1.0);

              fx2ptr[k] += fxk * 3.0 * invdr5;
              fy2ptr[k] += fyk * 3.0 * invdr5;
              fz2ptr[k] += fzk * 3.0 * invdr5;
            }
          }
        }
      }

      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
  }

  template <bool newton3>
  void SoAFunctorTripleImplScalar(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                  autopas::SoAView<SoAArraysType> soa3) {
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorTriple() is not implemented.");
    //   TODO
    if (soa1.size() == 0 or soa2.size() == 0 or soa3.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x3ptr = soa3.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y3ptr = soa3.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z3ptr = soa3.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState3ptr = soa3.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx3ptr = soa3.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy3ptr = soa3.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz3ptr = soa3.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type3ptr = soa3.template begin<Particle::AttributeNames::typeId>();

    std::vector<unsigned int> indicesK;

    for (size_t i = 0; i < soa1.size(); ++i) {
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      for (size_t j = 0; j < soa2.size(); ++j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy) {
          continue;
        }

        const SoAFloatPrecision drxij = x2ptr[j] - x1ptr[i];
        const SoAFloatPrecision dryij = y2ptr[j] - y1ptr[i];
        const SoAFloatPrecision drzij = z2ptr[j] - z1ptr[i];

        const SoAFloatPrecision drxij2 = drxij * drxij;
        const SoAFloatPrecision dryij2 = dryij * dryij;
        const SoAFloatPrecision drzij2 = drzij * drzij;

        const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
        if (drij2 > _cutoffSquared) {
          continue;
        }

        indicesK.clear();
        for (size_t k = 0; k < soa3.size(); ++k) {
          if (ownedState3ptr[k] == autopas::OwnershipState::dummy) {
            continue;
          }

          const SoAFloatPrecision drxjk = x3ptr[k] - x2ptr[j];
          const SoAFloatPrecision dryjk = y3ptr[k] - y2ptr[j];
          const SoAFloatPrecision drzjk = z3ptr[k] - z2ptr[j];

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          const SoAFloatPrecision drxki = x1ptr[i] - x3ptr[k];
          const SoAFloatPrecision dryki = y1ptr[i] - y3ptr[k];
          const SoAFloatPrecision drzki = z1ptr[i] - z3ptr[k];

          const SoAFloatPrecision drxki2 = drxki * drxki;
          const SoAFloatPrecision dryki2 = dryki * dryki;
          const SoAFloatPrecision drzki2 = drzki * drzki;

          const SoAFloatPrecision drki2 = drxki2 + dryki2 + drzki2;

          if (drjk2 > _cutoffSquared or drki2 > _cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision drxi2 = drxij * drxki;
          const SoAFloatPrecision dryi2 = dryij * dryki;
          const SoAFloatPrecision drzi2 = drzij * drzki;

          const SoAFloatPrecision dri2 = drxi2 + dryi2 + drzi2;

          const SoAFloatPrecision drxj2 = drxij * drxjk;
          const SoAFloatPrecision dryj2 = dryij * dryjk;
          const SoAFloatPrecision drzj2 = drzij * drzjk;

          const SoAFloatPrecision drj2 = drxj2 + dryj2 + drzj2;

          const SoAFloatPrecision drxk2 = drxjk * drxki;
          const SoAFloatPrecision dryk2 = dryjk * dryki;
          const SoAFloatPrecision drzk2 = drzjk * drzki;

          const SoAFloatPrecision drk2 = drxk2 + dryk2 + drzk2;

          const SoAFloatPrecision drijk2 = dri2 * drj2 * drk2;

          const SoAFloatPrecision dr2 = drij2 * drjk2 * drki2;
          const SoAFloatPrecision dr5 = dr2 * dr2 * std::sqrt(dr2);

          auto nu = _nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(type1ptr[i], type2ptr[j], type3ptr[k]);
          }
          const SoAFloatPrecision invdr5 = nu / dr5;

          const auto fxi = drxjk * dri2 * (drj2 - drk2) + drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           drxki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fxacc += fxi * 3.0 * invdr5;
          const auto fyi = dryjk * dri2 * (drj2 - drk2) + dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           dryki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fyacc += fyi * 3.0 * invdr5;
          const auto fzi = drzjk * dri2 * (drj2 - drk2) + drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
                           drzki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2);
          fzacc += fzi * 3.0 * invdr5;

          if (newton3) {
            const auto fxj = drxki * drj2 * (drk2 - dri2) +
                             drxij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                             drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
            const auto fyj = dryki * drj2 * (drk2 - dri2) +
                             dryij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                             dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);
            const auto fzj = drzki * drj2 * (drk2 - dri2) +
                             drzij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                             drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2);

            fx2ptr[j] += fxj * 3.0 * invdr5;
            fy2ptr[j] += fyj * 3.0 * invdr5;
            fz2ptr[j] += fzj * 3.0 * invdr5;

            /*const auto fxk = drxij * drk2 * (dri2 - drj2) +
                             drxjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                             drxki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
            const auto fyk = dryij * drk2 * (dri2 - drj2) +
                             dryjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                             dryki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);
            const auto fzk = drzij * drk2 * (dri2 - drj2) +
                             drzjk * (-dri2 * drj2 + drij2 * drki2 - 5.0 * drijk2 / drjk2) +
                             drzki * (dri2 * drj2 - drij2 * drjk2 + 5.0 * drijk2 / drki2);*/
            const auto fxk = (fxi + fxj) * (-1.0);
            const auto fyk = (fyi + fyj) * (-1.0);
            const auto fzk = (fzi + fzj) * (-1.0);

            fx3ptr[k] += fxk * 3.0 * invdr5;
            fy3ptr[k] += fyk * 3.0 * invdr5;
            fz3ptr[k] += fzk * 3.0 * invdr5;
          }
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
  }

 public:
  // clang-format off
  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorVerlet() is not implemented.");
  }

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param nu
   */
  void setParticleProperties(SoAFloatPrecision nu) { _nu = nu; }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 9>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::id,   Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

  /**
   *
   * @return useMixing
   */
  constexpr static bool getMixing() { return useMixing; }

  /**
   * Get the number of flops used per kernel call for a given particle pair. This should count the
   * floating point operations needed for two particles that lie within a cutoff radius, having already calculated the
   * distance.
   * @param molAType molecule A's type id
   * @param molBType molecule B's type id
   * @param molCType molecule C's type id
   * @param newton3 is newton3 applied.
   * @note The molecule types make no difference for AxilrodTellerFunctor, but are kept to have a consistent interface
   * for other functors where they may.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, size_t molCType, bool newton3) {
    //
    // Kernel: 56 = 18 (three dot products) + 9 (coefficients) + 29 (force calculation) sum
    // Adding to particle forces: 3
    // For Newton3: 29 (second force calculation) + 3 (adding force) + 6 (adding force to third p)
    // Total = 56 + 3 + ( 29 + 3 + 6 ) = 59 or 97
    return newton3 ? 97ul : 59ul;
  }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() final {
    _potentialEnergySum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
    }
  }

  /**
   * Accumulates global values, e.g. potential energy and virial.
   * @param newton3
   */
  void endTraversal(bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    if (calculateGlobals) {
      // Accumulate potential energy and virial values.
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
        _virialSum += _aosThreadData[i].virialSum;
      }

      _postProcessed = true;

      AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential Energy.
   * @return the potential Energy
   */
  double getPotentialEnergy() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get potential energy even though calculateGlobals is false. If you want this functor to calculate "
          "global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Cannot get potential energy, because endTraversal was not called.");
    }
    return _potentialEnergySum;
  }

  /**
   * Get the virial.
   * @return
   */
  double getVirial() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get virial even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Cannot get virial, because endTraversal was not called.");
    }
    return _virialSum[0] + _virialSum[1] + _virialSum[2];
  }

 private:
  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorVerletImpl() is not implemented.");
  }

  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    // variables
    std::array<double, 3> virialSum;
    double potentialEnergySum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  const double _cutoffSquared;

  // not const because they might be reset through PPL
  double _nu = 0.0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  const simde__m512d _zero{simde_mm512_setzero_pd()};
  const simde__m512i _zeroI{simde_mm512_setzero_si512()};
  const simde__m512d _three{simde_mm512_set1_pd(3.0)};
  const simde__m512d _five{simde_mm512_set1_pd(5.0)};
  simde__m512d _nuPd;

  static constexpr int vecLength = 8;
  const simde__mmask8 _masks[vecLength]{0b00000000, 0b00000001, 0b00000011, 0b00000111,
                                        0b00001111, 0b00011111, 0b00111111, 0b01111111};

  simde__m512i _ascendingIndices{simde_mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)};
  const simde__m512d _cutoffSquaredPd;
  const simde__m512i _ownedStateDummyEpi64{
      simde_mm512_set1_epi64(static_cast<int64_t>(autopas::OwnershipState::dummy))};
};
}  // namespace mdLib

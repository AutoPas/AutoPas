/**
 * @file LJFunctorSVE.h
 *
 * @date 15.11.2021
 * @author T. Eke
 */
#pragma once

#ifndef __ARM_FEATURE_SVE
#pragma message "Requested to compile LJFunctorSVE but SVE is not available!"
#else
#include <arm_sve.h>
#endif

#include <array>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * This Version is implemented using SVE intrinsics.
 * @tparam Particle The type of particle.
 * @tparam ParticleCell The type of particlecell.
 * @tparam applyShift Switch for the lj potential to be truncated shifted.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>
class LJFunctorSVE
    : public autopas::PairwiseFunctor<
          Particle, LJFunctorSVE<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorSVE() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit LJFunctorSVE(double cutoff, void * /*dummy*/)
#ifdef __ARM_FEATURE_SVE
      : autopas::PairwiseFunctor<
            Particle, LJFunctorSVE<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},
        _cutoffSquaredAoS(cutoff * cutoff),
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
  }
#else
      : autopas::PairwiseFunctor<
            Particle, LJFunctorSVE<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff) {
    autopas::utils::ExceptionHandler::exception("AutoPas was compiled without SVE support!");
  }
#endif
 public:
  /**
   * Constructor for Functor with mixing disabled. When using this functor it is necessary to call
   * setParticleProperties() to set internal constants because it does not use a particle properties library.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit LJFunctorSVE(double cutoff) : LJFunctorSVE(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like sigma, epsilon and shift.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit LJFunctorSVE(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorSVE(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "LJFunctorSVE"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto sigmaSquared = _sigmaSquaredAoS;
    auto epsilon24 = _epsilon24AoS;
    auto shift6 = _shift6AoS;
    if constexpr (useMixing) {
      sigmaSquared = _PPLibrary->getMixingSigmaSquared(i.getTypeId(), j.getTypeId());
      epsilon24 = _PPLibrary->getMixing24Epsilon(i.getTypeId(), j.getTypeId());
      if constexpr (applyShift) {
        shift6 = _PPLibrary->getMixingShift6(i.getTypeId(), j.getTypeId());
      }
    }
    auto dr = i.getR() - j.getR();
    double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquaredAoS) {
      return;
    }

    double invdr2 = 1. / dr2;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = dr * fac;
    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }
    if (calculateGlobals) {
      auto virial = dr * f;
      double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      const int threadnum = autopas::autopas_get_thread_num();
      if (i.isOwned()) {
        if (newton3) {
          _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergy6 * 0.5;
          _aosThreadData[threadnum].virialSumN3 += virial * 0.5;
        } else {
          // for non-newton3 the division is in the post-processing step.
          _aosThreadData[threadnum].potentialEnergySumNoN3 += potentialEnergy6;
          _aosThreadData[threadnum].virialSumNoN3 += virial;
        }
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergy6 * 0.5;
        _aosThreadData[threadnum].virialSumN3 += virial * 0.5;
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always do a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImpl<true>(soa);
    } else {
      SoAFunctorSingleImpl<false>(soa);
    }
  }

  // clang-format off
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  // clang-format on
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

 private:
  /**
   * Templatized version of SoAFunctorSingle actually doing what the latter should.
   * @tparam newton3
   * @param soa
   */
  template <bool newton3>
  void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
#ifdef __ARM_FEATURE_SVE
    if (soa.size() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

    svfloat64_t virialSumX = svdup_f64(0.0);
    svfloat64_t virialSumY = svdup_f64(0.0);
    svfloat64_t virialSumZ = svdup_f64(0.0);
    svfloat64_t potentialEnergySum = svdup_f64(0.0);

    const auto vecLength = (size_t)svlen_f64(potentialEnergySum);

    // reverse outer loop s.th. inner loop always beginns at aligned array start
    // typecast to detect underflow
    for (size_t i = soa.size() - 1; (long)i >= 0; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");

      svfloat64_t fxacc = svdup_f64(0.0);
      svfloat64_t fyacc = svdup_f64(0.0);
      svfloat64_t fzacc = svdup_f64(0.0);

      const svfloat64_t x1 = svdup_f64(xptr[i]);
      const svfloat64_t y1 = svdup_f64(yptr[i]);
      const svfloat64_t z1 = svdup_f64(zptr[i]);

      svbool_t pg_1, pg_2, pg_3, pg_4;
      size_t j = 0;
      for (; j < i; j += vecLength * 4) {
        const size_t j_2 = j + vecLength;
        const size_t j_3 = j_2 + vecLength;
        const size_t j_4 = j_3 + vecLength;
        pg_1 = svwhilelt_b64(j, i);
        pg_2 = svwhilelt_b64(j_2, i);
        pg_3 = svwhilelt_b64(j_3, i);
        pg_4 = svwhilelt_b64(j_4, i);

        SoAKernel<true, false>(j, ownedStatePtr[i] == autopas::OwnershipState::owned,
                               reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr, zptr, fxptr,
                               fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc, virialSumX, virialSumY,
                               virialSumZ, potentialEnergySum, pg_1, svundef_u64(), pg_2, svundef_u64(), pg_3,
                               svundef_u64(), pg_4, svundef_u64());
      }

      fxptr[i] += svaddv(svptrue_b64(), fxacc);
      fyptr[i] += svaddv(svptrue_b64(), fyacc);
      fzptr[i] += svaddv(svptrue_b64(), fzacc);
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
      // false, since for newton3 disabled we divide by two later on.
      if (newton3) {
        _aosThreadData[threadnum].potentialEnergySumN3 += svaddv_f64(svptrue_b64(), potentialEnergySum) * 0.5;
        _aosThreadData[threadnum].virialSumN3[0] += svaddv_f64(svptrue_b64(), virialSumX) * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += svaddv_f64(svptrue_b64(), virialSumY) * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += svaddv_f64(svptrue_b64(), virialSumZ) * 0.5;
      } else {
        _aosThreadData[threadnum].potentialEnergySumNoN3 += svaddv_f64(svptrue_b64(), potentialEnergySum);
        _aosThreadData[threadnum].virialSumNoN3[0] += svaddv_f64(svptrue_b64(), virialSumX);
        _aosThreadData[threadnum].virialSumNoN3[1] += svaddv_f64(svptrue_b64(), virialSumY);
        _aosThreadData[threadnum].virialSumNoN3[2] += svaddv_f64(svptrue_b64(), virialSumZ);
      }
    }
#endif
  }

  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
#ifdef __ARM_FEATURE_SVE
    if (soa1.size() == 0 || soa2.size() == 0) return;

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeID1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    const auto *const __restrict typeID2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    svfloat64_t virialSumX = svdup_f64(0.0);
    svfloat64_t virialSumY = svdup_f64(0.0);
    svfloat64_t virialSumZ = svdup_f64(0.0);
    svfloat64_t potentialEnergySum = svdup_f64(0.0);

    const auto vecLength = (unsigned int)svlen_f64(potentialEnergySum);

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

      svfloat64_t fxacc = svdup_f64(0.0);
      svfloat64_t fyacc = svdup_f64(0.0);
      svfloat64_t fzacc = svdup_f64(0.0);

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");

      const svfloat64_t x1 = svdup_f64(x1ptr[i]);
      const svfloat64_t y1 = svdup_f64(y1ptr[i]);
      const svfloat64_t z1 = svdup_f64(z1ptr[i]);

      svbool_t pg_1, pg_2, pg_3, pg_4;
      unsigned int j = 0;
      for (; j < soa2.size(); j += vecLength * 4) {
        const unsigned int j_2 = j + vecLength;
        const unsigned int j_3 = j_2 + vecLength;
        const unsigned int j_4 = j_3 + vecLength;
        pg_1 = svwhilelt_b64(j, (unsigned int)soa2.size());
        pg_2 = svwhilelt_b64(j_2, (unsigned int)soa2.size());
        pg_3 = svwhilelt_b64(j_3, (unsigned int)soa2.size());
        pg_4 = svwhilelt_b64(j_4, (unsigned int)soa2.size());

        SoAKernel<newton3, false>(j, ownedStatePtr1[i] == autopas::OwnershipState::owned,
                                  reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr, y2ptr, z2ptr,
                                  fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc, virialSumX,
                                  virialSumY, virialSumZ, potentialEnergySum, pg_1, svundef_u64(), pg_2, svundef_u64(),
                                  pg_3, svundef_u64(), pg_4, svundef_u64());
      }

      fx1ptr[i] += svaddv_f64(svptrue_b64(), fxacc);
      fy1ptr[i] += svaddv_f64(svptrue_b64(), fyacc);
      fz1ptr[i] += svaddv_f64(svptrue_b64(), fzacc);
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      if (newton3) {
        _aosThreadData[threadnum].potentialEnergySumN3 += svaddv_f64(svptrue_b64(), potentialEnergySum) * 0.5;
        _aosThreadData[threadnum].virialSumN3[0] += svaddv_f64(svptrue_b64(), virialSumX) * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += svaddv_f64(svptrue_b64(), virialSumY) * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += svaddv_f64(svptrue_b64(), virialSumZ) * 0.5;
      } else {
        _aosThreadData[threadnum].potentialEnergySumNoN3 += svaddv_f64(svptrue_b64(), potentialEnergySum);
        _aosThreadData[threadnum].virialSumNoN3[0] += svaddv_f64(svptrue_b64(), virialSumX);
        _aosThreadData[threadnum].virialSumNoN3[1] += svaddv_f64(svptrue_b64(), virialSumY);
        _aosThreadData[threadnum].virialSumNoN3[2] += svaddv_f64(svptrue_b64(), virialSumZ);
      }
    }
#endif
  }
#ifdef __ARM_FEATURE_SVE
  template <bool indexed>
  inline svbool_t distCalc(const size_t j, const svuint64_t &index, const svfloat64_t &x1, const svfloat64_t &y1,
                           const svfloat64_t &z1, const svbool_t &pg, const int64_t *const __restrict ownedStatePtr2,
                           const double *const __restrict x2ptr, const double *const __restrict y2ptr,
                           const double *const __restrict z2ptr, svfloat64_t &drx, svfloat64_t &dry, svfloat64_t &drz,
                           svfloat64_t &dr2, svint64_t &ownedStateJ) {
    const svfloat64_t x2 = (indexed) ? svld1_gather_index(pg, x2ptr, index) : svld1(pg, &x2ptr[j]);
    const svfloat64_t y2 = (indexed) ? svld1_gather_index(pg, y2ptr, index) : svld1(pg, &y2ptr[j]);
    const svfloat64_t z2 = (indexed) ? svld1_gather_index(pg, z2ptr, index) : svld1(pg, &z2ptr[j]);

    // having these three as not _m worsens performance
    drx = svsub_m(pg, x1, x2);
    dry = svsub_m(pg, y1, y2);
    drz = svsub_m(pg, z1, z2);

    const svfloat64_t dr2_1 = svmul_x(pg, drx, drx);
    const svfloat64_t dr2_2 = svmla_x(pg, dr2_1, dry, dry);
    dr2 = svmla_x(pg, dr2_2, drz, drz);

    const svbool_t cutoffMask = svcmple(pg, dr2, _cutoffSquared);

    ownedStateJ = (indexed) ? svld1_gather_index(pg, ownedStatePtr2, index) : svld1(pg, &ownedStatePtr2[j]);
    const svbool_t dummyMask = svcmpne(pg, ownedStateJ, (int64_t)autopas::OwnershipState::dummy);
    return svand_z(pg, cutoffMask, dummyMask);
  }

  template <bool indexed>
  inline void lennardJones(const svuint64_t &index, const size_t *const typeID1ptr, const size_t *const typeID2ptr,
                           const svbool_t &pgC, const svfloat64_t &dr2, svfloat64_t &epsilon24s, svfloat64_t &shift6s,
                           svfloat64_t &lj6, svfloat64_t &fac) {
    const svuint64_t typeIds =
        useMixing ? svmul_m(pgC, (indexed) ? svld1_gather_index(pgC, typeID2ptr, index) : svld1_u64(pgC, typeID2ptr), 3)
                  : svundef_u64();
    const auto mixingDataPtr = useMixing ? _PPLibrary->getLJMixingDataPtr(*typeID1ptr, 0) : nullptr;

    const svfloat64_t sigmaSquareds =
        useMixing ? svld1_gather_index(pgC, mixingDataPtr + 1, typeIds) : svdup_f64(_sigmaSquared);
    epsilon24s = useMixing ? svld1_gather_index(pgC, mixingDataPtr, typeIds) : svdup_f64(_epsilon24);
    shift6s = (useMixing && applyShift) ? svld1_gather_index(pgC, mixingDataPtr + 2, typeIds) : svdup_f64(_shift6);

    svfloat64_t invdr2 = svrecpe(dr2);
    invdr2 = svmul_x(pgC, invdr2, svrecps(dr2, invdr2));
    invdr2 = svmul_x(pgC, invdr2, svrecps(dr2, invdr2));
    invdr2 = svmul_x(pgC, invdr2, svrecps(dr2, invdr2));
    invdr2 = svmul_x(pgC, invdr2, svrecps(dr2, invdr2));
    const svfloat64_t lj2 = svmul_x(pgC, sigmaSquareds, invdr2);

    const svfloat64_t lj4 = svmul_x(pgC, lj2, lj2);
    lj6 = svmul_x(pgC, lj2, lj4);
    const svfloat64_t lj12 = svmul_x(pgC, lj6, lj6);

    //(2*lj12 - lj6) * epsilon24 * invdr2 = -(lj6 - 2*lj12) * epsilon24 * invdr2

    const svfloat64_t lj12m6alj12 = svnmls_x(pgC, lj6, lj12, 2);

    const svfloat64_t lj12m6alj12e = svmul_x(pgC, lj12m6alj12, epsilon24s);
    fac = svmul_x(pgC, lj12m6alj12e, invdr2);
  }

  template <bool newton3, bool indexed>
  inline void applyForces(const size_t j, const svuint64_t &index, const bool ownedStateIisOwned,
                          double *const __restrict fx2ptr, double *const __restrict fy2ptr,
                          double *const __restrict fz2ptr, svfloat64_t &fxacc, svfloat64_t &fyacc, svfloat64_t &fzacc,
                          svfloat64_t &virialSumX, svfloat64_t &virialSumY, svfloat64_t &virialSumZ,
                          svfloat64_t &potentialEnergySum, const svfloat64_t &drx, const svfloat64_t &dry,
                          const svfloat64_t &drz, const double *const __restrict x2ptr,
                          const double *const __restrict y2ptr, const double *const __restrict z2ptr,
                          const svint64_t &ownedStateJ, const svbool_t &pgC, const svfloat64_t &epsilon24s,
                          const svfloat64_t &shift6s, const svfloat64_t &lj6, const svfloat64_t &fac) {
    const svfloat64_t fx = svmul_x(pgC, drx, fac);
    const svfloat64_t fy = svmul_x(pgC, dry, fac);
    const svfloat64_t fz = svmul_x(pgC, drz, fac);

    fxacc = svadd_f64_m(pgC, fxacc, fx);
    fyacc = svadd_f64_m(pgC, fyacc, fy);
    fzacc = svadd_f64_m(pgC, fzacc, fz);

    if (newton3) {
      const svfloat64_t fx2 = (indexed) ? svld1_gather_index(pgC, fx2ptr, index) : svld1_f64(pgC, &fx2ptr[j]);
      const svfloat64_t fy2 = (indexed) ? svld1_gather_index(pgC, fy2ptr, index) : svld1_f64(pgC, &fy2ptr[j]);
      const svfloat64_t fz2 = (indexed) ? svld1_gather_index(pgC, fz2ptr, index) : svld1_f64(pgC, &fz2ptr[j]);

      const svfloat64_t fx2new = svsub_x(pgC, fx2, fx);
      const svfloat64_t fy2new = svsub_x(pgC, fy2, fy);
      const svfloat64_t fz2new = svsub_x(pgC, fz2, fz);

      if constexpr (indexed) {
        svst1_scatter_index(pgC, &fx2ptr[0], index, fx2new);
        svst1_scatter_index(pgC, &fy2ptr[0], index, fy2new);
        svst1_scatter_index(pgC, &fz2ptr[0], index, fz2new);
      } else {
        svst1(pgC, &fx2ptr[j], fx2new);
        svst1(pgC, &fy2ptr[j], fy2new);
        svst1(pgC, &fz2ptr[j], fz2new);
      }
    }

    if (calculateGlobals) {
      // Global Virial
      const svfloat64_t virialX = svmul_x(pgC, fx, drx);
      const svfloat64_t virialY = svmul_x(pgC, fy, dry);
      const svfloat64_t virialZ = svmul_x(pgC, fz, drz);

      // Global Potential
      const svfloat64_t lj12m6 = svnmls_x(pgC, lj6, lj6, lj6);
      const svfloat64_t potentialEnergy6 = svmad_x(pgC, epsilon24s, lj12m6, shift6s);
      svfloat64_t energyFactor = svdup_f64(ownedStateIisOwned ? 1.0 : 0.0);

      if constexpr (newton3) {
        svbool_t ownedMaskJ = svcmpeq(pgC, ownedStateJ, (int64_t)autopas::OwnershipState::owned);
        energyFactor = svadd_m(ownedMaskJ, energyFactor, 1.0);
      }
      potentialEnergySum = svmla_m(pgC, potentialEnergySum, energyFactor, potentialEnergy6);
      virialSumX = svmla_m(pgC, virialSumX, energyFactor, virialX);
      virialSumY = svmla_m(pgC, virialSumY, energyFactor, virialY);
      virialSumZ = svmla_m(pgC, virialSumZ, energyFactor, virialZ);
    }
  }

  template <bool newton3, bool indexed>
  // FCC needs to be forced to inline this function. Otherwise a dramatic loss in performance can be observed.
  __attribute__((always_inline)) inline void SoAKernel(
      const size_t j, const bool ownedStateIisOwned, const int64_t *const __restrict ownedStatePtr2,
      const svfloat64_t &x1, const svfloat64_t &y1, const svfloat64_t &z1, const double *const __restrict x2ptr,
      const double *const __restrict y2ptr, const double *const __restrict z2ptr, double *const __restrict fx2ptr,
      double *const __restrict fy2ptr, double *const __restrict fz2ptr, const size_t *const typeID1ptr,
      const size_t *const typeID2ptr, svfloat64_t &fxacc, svfloat64_t &fyacc, svfloat64_t &fzacc,
      svfloat64_t &virialSumX, svfloat64_t &virialSumY, svfloat64_t &virialSumZ, svfloat64_t &potentialEnergySum,

      const svbool_t &pg_1, const svuint64_t &index_1, const svbool_t &pg_2, const svuint64_t &index_2,
      const svbool_t &pg_3, const svuint64_t &index_3, const svbool_t &pg_4, const svuint64_t &index_4

  ) {
    svfloat64_t drx_1;
    svfloat64_t dry_1;
    svfloat64_t drz_1;
    svfloat64_t dr2_1;
    svint64_t ownedStateJ_1;
    const svbool_t pgC_1 = distCalc<indexed>(j, index_1, x1, y1, z1, pg_1, ownedStatePtr2, x2ptr, y2ptr, z2ptr, drx_1,
                                             dry_1, drz_1, dr2_1, ownedStateJ_1);

    svfloat64_t drx_2;
    svfloat64_t dry_2;
    svfloat64_t drz_2;
    svfloat64_t dr2_2;
    svint64_t ownedStateJ_2;
    const svbool_t pgC_2 = distCalc<indexed>(j + svlen(x1), index_2, x1, y1, z1, pg_2, ownedStatePtr2, x2ptr, y2ptr,
                                             z2ptr, drx_2, dry_2, drz_2, dr2_2, ownedStateJ_2);

    svfloat64_t drx_3;
    svfloat64_t dry_3;
    svfloat64_t drz_3;
    svfloat64_t dr2_3;
    svint64_t ownedStateJ_3;
    const svbool_t pgC_3 = distCalc<indexed>(j + svlen(x1) * 2, index_3, x1, y1, z1, pg_3, ownedStatePtr2, x2ptr, y2ptr,
                                             z2ptr, drx_3, dry_3, drz_3, dr2_3, ownedStateJ_3);

    svfloat64_t drx_4;
    svfloat64_t dry_4;
    svfloat64_t drz_4;
    svfloat64_t dr2_4;
    svint64_t ownedStateJ_4;
    const svbool_t pgC_4 = distCalc<indexed>(j + svlen(x1) * 3, index_4, x1, y1, z1, pg_4, ownedStatePtr2, x2ptr, y2ptr,
                                             z2ptr, drx_4, dry_4, drz_4, dr2_4, ownedStateJ_4);

    const bool continue_1 = svptest_any(svptrue_b64(), pgC_1);
    const bool continue_2 = svptest_any(svptrue_b64(), pgC_2);
    const bool continue_3 = svptest_any(svptrue_b64(), pgC_3);
    const bool continue_4 = svptest_any(svptrue_b64(), pgC_4);

    svfloat64_t epsilon24s_1;
    svfloat64_t shift6s_1;
    svfloat64_t lj6_1;
    svfloat64_t fac_1;
    if (continue_1)
      lennardJones<indexed>(index_1, typeID1ptr, typeID2ptr, pgC_1, dr2_1, epsilon24s_1, shift6s_1, lj6_1, fac_1);

    svfloat64_t epsilon24s_2;
    svfloat64_t shift6s_2;
    svfloat64_t lj6_2;
    svfloat64_t fac_2;
    if (continue_2)
      lennardJones<indexed>(index_2, typeID1ptr, typeID2ptr, pgC_2, dr2_2, epsilon24s_2, shift6s_2, lj6_2, fac_2);

    svfloat64_t epsilon24s_3;
    svfloat64_t shift6s_3;
    svfloat64_t lj6_3;
    svfloat64_t fac_3;
    if (continue_3)
      lennardJones<indexed>(index_3, typeID1ptr, typeID2ptr, pgC_3, dr2_3, epsilon24s_3, shift6s_3, lj6_3, fac_3);

    svfloat64_t epsilon24s_4;
    svfloat64_t shift6s_4;
    svfloat64_t lj6_4;
    svfloat64_t fac_4;
    if (continue_4)
      lennardJones<indexed>(index_4, typeID1ptr, typeID2ptr, pgC_4, dr2_4, epsilon24s_4, shift6s_4, lj6_4, fac_4);

    if (continue_1)
      applyForces<newton3, indexed>(j, index_1, ownedStateIisOwned, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc,
                                    virialSumX, virialSumY, virialSumZ, potentialEnergySum, drx_1, dry_1, drz_1, x2ptr,
                                    y2ptr, z2ptr, ownedStateJ_1, pgC_1, epsilon24s_1, shift6s_1, lj6_1, fac_1);
    if (continue_2)
      applyForces<newton3, indexed>(j + svlen(x1), index_2, ownedStateIisOwned, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc,
                                    fzacc, virialSumX, virialSumY, virialSumZ, potentialEnergySum, drx_2, dry_2, drz_2,
                                    x2ptr, y2ptr, z2ptr, ownedStateJ_2, pgC_2, epsilon24s_2, shift6s_2, lj6_2, fac_2);

    if (continue_3)
      applyForces<newton3, indexed>(j + svlen(x1) * 2, index_3, ownedStateIisOwned, fx2ptr, fy2ptr, fz2ptr, fxacc,
                                    fyacc, fzacc, virialSumX, virialSumY, virialSumZ, potentialEnergySum, drx_3, dry_3,
                                    drz_3, x2ptr, y2ptr, z2ptr, ownedStateJ_3, pgC_3, epsilon24s_3, shift6s_3, lj6_3,
                                    fac_3);

    if (continue_4)
      applyForces<newton3, indexed>(j + svlen(x1) * 3, index_4, ownedStateIisOwned, fx2ptr, fy2ptr, fz2ptr, fxacc,
                                    fyacc, fzacc, virialSumX, virialSumY, virialSumZ, potentialEnergySum, drx_4, dry_4,
                                    drz_4, x2ptr, y2ptr, z2ptr, ownedStateJ_4, pgC_4, epsilon24s_4, shift6s_4, lj6_4,
                                    fac_4);
  }
#endif

 public:
  // clang-format off
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

 private:
  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
#ifdef __ARM_FEATURE_SVE
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();
    if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
      return;
    }
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

    // accumulators
    svfloat64_t virialSumX = svdup_f64(0.0);
    svfloat64_t virialSumY = svdup_f64(0.0);
    svfloat64_t virialSumZ = svdup_f64(0.0);
    svfloat64_t potentialEnergySum = svdup_f64(0.0);

    svfloat64_t fxacc = svdup_f64(0.0);
    svfloat64_t fyacc = svdup_f64(0.0);
    svfloat64_t fzacc = svdup_f64(0.0);

    // broadcast particle 1
    const auto x1 = svdup_f64(xptr[indexFirst]);
    const auto y1 = svdup_f64(yptr[indexFirst]);
    const auto z1 = svdup_f64(zptr[indexFirst]);

    svbool_t pg_1;
    const auto *const ownedStatePtr2 = reinterpret_cast<const int64_t *>(ownedStatePtr);
    size_t j = 0;
    for (; j < neighborList.size(); j += svlen(x1)) {
      pg_1 = svwhilelt_b64(j, neighborList.size());
      const svuint64_t index_1 = svld1(pg_1, &neighborList[j]);

      svfloat64_t drx_1;
      svfloat64_t dry_1;
      svfloat64_t drz_1;
      svfloat64_t dr2_1;
      svint64_t ownedStateJ_1;
      const svbool_t pgC_1 = distCalc<true>(0, index_1, x1, y1, z1, pg_1, ownedStatePtr2, xptr, yptr, zptr, drx_1,
                                            dry_1, drz_1, dr2_1, ownedStateJ_1);

      const bool continue_1 = svptest_any(svptrue_b64(), pgC_1);
      svfloat64_t epsilon24s_1;
      svfloat64_t shift6s_1;
      svfloat64_t lj6_1;
      svfloat64_t fac_1;
      if (continue_1)
        lennardJones<true>(index_1, typeIDptr, typeIDptr, pgC_1, dr2_1, epsilon24s_1, shift6s_1, lj6_1, fac_1);

      if (continue_1)
        applyForces<newton3, true>(0, index_1, ownedStatePtr[indexFirst] == autopas::OwnershipState::owned, fxptr,
                                   fyptr, fzptr, fxacc, fyacc, fzacc, virialSumX, virialSumY, virialSumZ,
                                   potentialEnergySum, drx_1, dry_1, drz_1, xptr, yptr, zptr, ownedStateJ_1, pgC_1,
                                   epsilon24s_1, shift6s_1, lj6_1, fac_1);
    }

    fxptr[indexFirst] += svaddv_f64(svptrue_b64(), fxacc);
    fyptr[indexFirst] += svaddv_f64(svptrue_b64(), fyacc);
    fzptr[indexFirst] += svaddv_f64(svptrue_b64(), fzacc);

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      if (newton3) {
        _aosThreadData[threadnum].potentialEnergySumN3 += svaddv_f64(svptrue_b64(), potentialEnergySum) * 0.5;
        _aosThreadData[threadnum].virialSumN3[0] += svaddv_f64(svptrue_b64(), virialSumX) * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += svaddv_f64(svptrue_b64(), virialSumY) * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += svaddv_f64(svptrue_b64(), virialSumZ) * 0.5;
      } else {
        _aosThreadData[threadnum].potentialEnergySumNoN3 += svaddv_f64(svptrue_b64(), potentialEnergySum);
        _aosThreadData[threadnum].virialSumNoN3[0] += svaddv_f64(svptrue_b64(), virialSumX);
        _aosThreadData[threadnum].virialSumNoN3[1] += svaddv_f64(svptrue_b64(), virialSumY);
        _aosThreadData[threadnum].virialSumNoN3[2] += svaddv_f64(svptrue_b64(), virialSumZ);
      }
    }
#endif
  }

 public:
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
   * @param newton3 is newton3 applied.
   * @note molAType and molBType make no difference for LJFunctor, but are kept to have a consistent interface for other
   * functors where they may.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, bool newton3) {
    // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply scale) sum
    // Adding to particle forces: 6 or 3 depending newton3
    // Total = 12 + (6 or 3) = 18 or 15
    return newton3 ? 18ul : 15ul;
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
   * Accumulates global values, e.g. potentialEnergy and virial.
   * @param newton3
   */
  void endTraversal(bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    if (calculateGlobals) {
      // We distinguish between non-newton3 and newton3 functor calls. Newton3 calls are accumulated directly.
      // Non-newton3 calls are accumulated temporarily and later divided by 2.
      double potentialEnergySumNoN3Acc = 0;
      std::array<double, 3> virialSumNoN3Acc = {0, 0, 0};
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        potentialEnergySumNoN3Acc += _aosThreadData[i].potentialEnergySumNoN3;
        _potentialEnergySum += _aosThreadData[i].potentialEnergySumN3;

        virialSumNoN3Acc += _aosThreadData[i].virialSumNoN3;
        _virialSum += _aosThreadData[i].virialSumN3;
      }
      // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
      // here.
      potentialEnergySumNoN3Acc *= 0.5;
      virialSumNoN3Acc *= 0.5;

      _potentialEnergySum += potentialEnergySumNoN3Acc;
      _virialSum += virialSumNoN3Acc;

      // we have always calculated 6*potentialEnergy, so we divide by 6 here!
      _potentialEnergySum /= 6.;
      _postProcessed = true;

      AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential Energy
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
   * Get the virial
   * @return the virial
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

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param epsilon24
   * @param sigmaSquared
   */
  void setParticleProperties(double epsilon24, double sigmaSquared) {
#ifdef __ARM_FEATURE_SVE
    _epsilon24 = epsilon24;
    _sigmaSquared = sigmaSquared;
    if constexpr (applyShift) {
      _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared);
    } else {
      _shift6 = 0.0;
    }
#endif

    _epsilon24AoS = epsilon24;
    _sigmaSquaredAoS = sigmaSquared;
    if constexpr (applyShift) {
      _shift6AoS = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquaredAoS);
    } else {
      _shift6AoS = 0.;
    }
  }

 private:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData()
        : virialSumNoN3{0., 0., 0.},
          virialSumN3{0., 0., 0.},
          potentialEnergySumNoN3{0.},
          potentialEnergySumN3{0.},
          __remainingTo64{} {}
    void setZero() {
      virialSumNoN3 = {0., 0., 0.};
      virialSumN3 = {0., 0., 0.};
      potentialEnergySumNoN3 = 0.;
      potentialEnergySumN3 = 0.;
    }

    // variables
    std::array<double, 3> virialSumNoN3;
    std::array<double, 3> virialSumN3;
    double potentialEnergySumNoN3;
    double potentialEnergySumN3;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 8 * sizeof(double)) / sizeof(double)];
  };
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

#ifdef __ARM_FEATURE_SVE
  const double _cutoffSquared{};
  double _shift6{0.};
  double _epsilon24{0.};
  double _sigmaSquared{0.};
#endif

  const double _cutoffSquaredAoS;
  double _epsilon24AoS{0.}, _sigmaSquaredAoS{0.}, _shift6AoS{0.};

  ParticlePropertiesLibrary<double, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum{0.};

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum{0., 0., 0.};

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData{};

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed{false};
};
}  // namespace mdLib

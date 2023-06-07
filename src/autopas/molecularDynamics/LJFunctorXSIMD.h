//
// Created by robin_3kolqs9 on 17.05.2023.
//

#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "xsimd/xsimd.hpp"



namespace autopas {
/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * This Version is implemented using the XSIMD Wrapper.
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
          FunctorN3Modes useNewton3 = FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>

class LJFunctorXSIMD
    : public Functor<Particle,
                     LJFunctorXSIMD<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
    using SoAArraysType = typename Particle::SoAArraysType;

   public:
    /**
   * Deleted default constructor
     */
     LJFunctorXSIMD() = delete;

    private:
     /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor.
      */
     explicit LJFunctorXSIMD(double cutoff, void * /*dummy*/)
         : Functor<Particle,
                   LJFunctorXSIMD<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(cutoff),
           _cutoffsquare{cutoff * cutoff},
           _cutoffsquareAoS(cutoff * cutoff),
           _upotSum{0.},
           _virialSum{0., 0., 0.},
           _aosThreadData(),
           _postProcessed{false} {
       if (calculateGlobals) {
         _aosThreadData.resize(autopas_get_max_threads());
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
     explicit LJFunctorXSIMD(double cutoff) : LJFunctorXSIMD(cutoff, nullptr) {
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
     explicit LJFunctorXSIMD(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
         : LJFunctorXSIMD(cutoff, nullptr) {
       static_assert(useMixing,
                     "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                     "or set mixing to true.");
       _PPLibrary = &particlePropertiesLibrary;
     }

     bool isRelevantForTuning() final { return relevantForTuning; }

     bool allowsNewton3() final { return useNewton3 == FunctorN3Modes::Newton3Only or useNewton3 == FunctorN3Modes::Both; }

     bool allowsNonNewton3() final {
       return useNewton3 == FunctorN3Modes::Newton3Off or useNewton3 == FunctorN3Modes::Both;
     }

     inline void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
       using namespace autopas::utils::ArrayMath::literals;
       if (i.isDummy() or j.isDummy()) {
         return;
       }
       auto sigmasquare = _sigmaSquareAoS;
       auto epsilon24 = _epsilon24AoS;
       auto shift6 = _shift6AoS;
       if constexpr (useMixing) {
         sigmasquare = _PPLibrary->mixingSigmaSquare(i.getTypeId(), j.getTypeId());
         epsilon24 = _PPLibrary->mixing24Epsilon(i.getTypeId(), j.getTypeId());
         if constexpr (applyShift) {
           shift6 = _PPLibrary->mixingShift6(i.getTypeId(), j.getTypeId());
         }
       }
       auto dr = i.getR() - j.getR();
       double dr2 = utils::ArrayMath::dot(dr, dr);

       if (dr2 > _cutoffsquareAoS) {
         return;
       }


       double invdr2 = 1. / dr2;
       double lj6 = sigmasquare * invdr2;
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
         double upot = epsilon24 * lj12m6 + shift6;

         const int threadnum = autopas_get_thread_num();
         // for non-newton3 the division is in the post-processing step.
         if (newton3) {
           upot *= 0.5;
           virial *= (double)0.5;
         }
         if (i.isOwned()) {
           _aosThreadData[threadnum].upotSum += upot;
           _aosThreadData[threadnum].virialSum += virial;
         }
         // for non-newton3 the second particle will be considered in a separate calculation
         if (newton3 and j.isOwned()) {
           _aosThreadData[threadnum].upotSum += upot;
           _aosThreadData[threadnum].virialSum += virial;
         }
       }
     }

     /**
   * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
   * This functor will always do a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
      */
     inline void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) final {
       if (newton3) {
         SoAFunctorSingleImpl<true>(soa);
       } else {
         SoAFunctorSingleImpl<false>(soa);
       }
     }

     // clang-format off
  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
     // clang-format on
     inline void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, const bool newton3) final {
       if (newton3) {
         SoAFunctorPairImpl<true>(soa1, soa2);
       } else {
         SoAFunctorPairImpl<false>(soa1, soa2);
       }
     }



    private:
     template <bool newton3>
     inline void SoAFunctorSingleImpl(SoAView<SoAArraysType> soa) {
       if (soa.getNumberOfParticles() == 0) return;

       const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
       const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
       const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

       const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

       auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
       auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
       auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

       const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

       xsimd::batch<double> virialSumX{0};
       xsimd::batch<double> virialSumY{0};
       xsimd::batch<double> virialSumZ{0};
       xsimd::batch<double> upotSum{0};

       for (size_t i = soa.getNumberOfParticles() - 1; (long)i >= 0; --i) {
         if (ownedStatePtr[i] == OwnershipState::dummy) {
           // If the i-th particle is a dummy, skip this loop iteration.
           continue;
         }
         static_assert(std::is_same_v<std::underlying_type_t<OwnershipState>, int64_t>,
                       "OwnershipStates underlying type should be int64_t!");

         xsimd::batch<int64_t> ownedStateI{static_cast<int64_t>(ownedStatePtr[i])};

         xsimd::batch<double> fxacc{0};
         xsimd::batch<double> fyacc{0};
         xsimd::batch<double> fzacc{0};

         const xsimd::batch<double> x1 = xsimd::broadcast(&xptr[i]);
         const xsimd::batch<double> y1 = xsimd::broadcast(&yptr[i]);
         const xsimd::batch<double> z1 = xsimd::broadcast(&zptr[i]);

         size_t j = 0;

         // floor soa numParticles to multiple of vecLength
         // If b is a power of 2 the following holds:
         // a & ~(b -1) == a - (a mod b)
         for (; j < (i & ~(vecLength - 1)); j += 4) {
           SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr,
                                  yptr, zptr, fxptr, fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc,
                                  &virialSumX, &virialSumY, &virialSumZ, &upotSum, 0);
         }
         // If b is a power of 2 the following holds:
         // a & (b -1) == a mod b
         const int rest = (int)(i & (vecLength - 1));
         if (rest > 0) {
           SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr,
                                 yptr, zptr, fxptr, fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc,
                                 &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);
         }


         double sumfx = xsimd::reduce_add(fxacc);
         double sumfy = xsimd::reduce_add(fyacc);
         double sumfz = xsimd::reduce_add(fzacc);

         fxptr[i] += sumfx;
         fyptr[i] += sumfy;
         fzptr[i] += sumfz;
       }
     }


     template <bool newton3>
     inline void SoAFunctorPairImpl(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2) {
       if (soa1.getNumberOfParticles() == 0 || soa2.getNumberOfParticles() == 0) return;

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

       xsimd::batch<double> virialSumX{0};
       xsimd::batch<double> virialSumY{0};
       xsimd::batch<double> virialSumZ{0};
       xsimd::batch<double> upotSum{0};

       for (unsigned int i = 0; i < soa1.getNumberOfParticles(); ++i) {
         if (ownedStatePtr1[i] == OwnershipState::dummy) {
           // If the i-th particle is a dummy, skip this loop iteration.
           continue;
         }

         xsimd::batch<double> fxacc{0};
         xsimd::batch<double> fyacc{0};
         xsimd::batch<double> fzacc{0};

         static_assert(std::is_same_v<std::underlying_type_t<OwnershipState>, int64_t>,
                       "OwnershipStates underlying type should be int64_t!");
         // ownedStatePtr1 contains int64_t, so we broadcast these to make an __m256i.
         // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
         xsimd::batch<int64_t> ownedStateI{static_cast<int64_t>(ownedStatePtr1[i])};

         const xsimd::batch<double> x1 = xsimd::broadcast(&x1ptr[i]);
         const xsimd::batch<double> y1 = xsimd::broadcast(&y1ptr[i]);
         const xsimd::batch<double> z1 = xsimd::broadcast(&z1ptr[i]);

         // floor soa2 numParticles to multiple of vecLength
         unsigned int j = 0;
         for (; j < (soa2.getNumberOfParticles() & ~(vecLength - 1)); j += vecLength) {
           SoAKernel<newton3, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                     y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                     &virialSumX, &virialSumY, &virialSumZ, &upotSum, 0);
         }
         const int rest = (int)(soa2.getNumberOfParticles() & (vecLength - 1));
         if (rest > 0)
           SoAKernel<newton3, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                    y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                    &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);

         fx1ptr[i] += xsimd::reduce_add(fxacc);
         fy1ptr[i] += xsimd::reduce_add(fyacc);
         fz1ptr[i] += xsimd::reduce_add(fzacc);
       }

       if constexpr (calculateGlobals) {
         //TODO
//         const int threadnum = autopas_get_thread_num();
//
//         // horizontally reduce virialSumX and virialSumY
//         const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
//         const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
//         const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
//         const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);
//
//         // horizontally reduce virialSumZ and upotSum
//         const __m256d hSumVirialzUpot = _mm256_hadd_pd(virialSumZ, upotSum);
//         const __m128d hSumVirialzUpotLow = _mm256_extractf128_pd(hSumVirialzUpot, 0);
//         const __m128d hSumVirialzUpotHigh = _mm256_extractf128_pd(hSumVirialzUpot, 1);
//         const __m128d hSumVirialzUpotVec = _mm_add_pd(hSumVirialzUpotHigh, hSumVirialzUpotLow);
//
//         // globals = {virialX, virialY, virialZ, uPot}
//         double globals[4];
//         _mm_store_pd(&globals[0], hSumVirialxyVec);
//         _mm_store_pd(&globals[2], hSumVirialzUpotVec);
//
//         // we have duplicated calculations, i.e., we calculate interactions multiple times, so we have to take care
//         // that we do not add the energy multiple times!
//         double energyfactor = 1.;
//         if constexpr (newton3) {
//           energyfactor *= 0.5;  // we count the energies partly to one of the two cells!
//         }
//
//         _aosThreadData[threadnum].virialSum[0] += globals[0] * energyfactor;
//         _aosThreadData[threadnum].virialSum[1] += globals[1] * energyfactor;
//         _aosThreadData[threadnum].virialSum[2] += globals[2] * energyfactor;
//         _aosThreadData[threadnum].upotSum += globals[3] * energyfactor;
       }
     }
     /**
   * Actual inner kernel of the SoAFunctors.
   *
   * @tparam newton3
   * @tparam remainderIsMasked If false the full vector length is used. Otherwise the last entries are masked away
   * depending on the argument "rest".
   * @param j
   * @param ownedStateI
   * @param ownedStatePtr2
   * @param x1
   * @param y1
   * @param z1
   * @param x2ptr
   * @param y2ptr
   * @param z2ptr
   * @param fx2ptr
   * @param fy2ptr
   * @param fz2ptr
   * @param typeID1ptr
   * @param typeID2ptr
   * @param fxacc
   * @param fyacc
   * @param fzacc
   * @param virialSumX
   * @param virialSumY
   * @param virialSumZ
   * @param upotSum
   * @param rest
      */
     template <bool newton3, bool remainderIsMasked>
     inline void SoAKernel(const size_t j, const __m256i ownedStateI, const xsimd::batch<int64_t> *const __restrict ownedStatePtr2,
                           const xsimd::batch<float> &x1, const xsimd::batch<float> &y1, const xsimd::batch<float> &z1, const double *const __restrict x2ptr,
                           const double *const __restrict y2ptr, const double *const __restrict z2ptr,
                           double *const __restrict fx2ptr, double *const __restrict fy2ptr,
                           double *const __restrict fz2ptr, const size_t *const typeID1ptr, const size_t *const typeID2ptr,
                           xsimd::batch<float> &fxacc, xsimd::batch<float> &fyacc, xsimd::batch<float> &fzacc, xsimd::batch<float> *virialSumX, xsimd::batch<float> *virialSumY,
                           xsimd::batch<float> *virialSumZ, xsimd::batch<float> *upotSum, const unsigned int rest = 0) {
       xsimd::batch<double> epsilon24s = _epsilon24;
       xsimd::batch<double> sigmaSquares = _sigmaSquare;
       xsimd::batch<double> shift6s = _shift6;
       if (useMixing) {
         // the first argument for set lands in the last bits of the register
         epsilon24s = {
             not remainderIsMasked or rest > 3 ? _PPLibrary->mixing24Epsilon(*typeID1ptr, *(typeID2ptr + 3)) : 0,
             not remainderIsMasked or rest > 2 ? _PPLibrary->mixing24Epsilon(*typeID1ptr, *(typeID2ptr + 2)) : 0,
             not remainderIsMasked or rest > 1 ? _PPLibrary->mixing24Epsilon(*typeID1ptr, *(typeID2ptr + 1)) : 0,
             _PPLibrary->mixing24Epsilon(*typeID1ptr, *(typeID2ptr + 0))};
         sigmaSquares = {
             not remainderIsMasked or rest > 3 ? _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(typeID2ptr + 3)) : 0,
             not remainderIsMasked or rest > 2 ? _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(typeID2ptr + 2)) : 0,
             not remainderIsMasked or rest > 1 ? _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(typeID2ptr + 1)) : 0,
             _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(typeID2ptr + 0))};
         if constexpr (applyShift) {
           shift6s = {
               (not remainderIsMasked or rest > 3) ? _PPLibrary->mixingShift6(*typeID1ptr, *(typeID2ptr + 3)) : 0,
               (not remainderIsMasked or rest > 2) ? _PPLibrary->mixingShift6(*typeID1ptr, *(typeID2ptr + 2)) : 0,
               (not remainderIsMasked or rest > 1) ? _PPLibrary->mixingShift6(*typeID1ptr, *(typeID2ptr + 1)) : 0,
               _PPLibrary->mixingShift6(*typeID1ptr, *(typeID2ptr + 0))};
         }
       }
       const xsimd::batch<double> x2 = remainderIsMasked ? xsimd::bitwise_and(xsimd::load(&x2ptr[j]), _masks[rest - 1]) : xsimd::load(&x2ptr[j]);
       const xsimd::batch<double> y2 = remainderIsMasked ? xsimd::bitwise_and(xsimd::load(&y2ptr[j]), _masks[rest - 1]) : xsimd::load(&x2ptr[j]);
       const xsimd::batch<double> z2 = remainderIsMasked ? xsimd::bitwise_and(xsimd::load(&z2ptr[j]), _masks[rest - 1]) : xsimd::load(&x2ptr[j]);

       const xsimd::batch<double> drx = xsimd::sub(x1,x2);
       const xsimd::batch<double> dry = xsimd::sub(y1,y2);
       const xsimd::batch<double> drz = xsimd::sub(z1,z2);

       const xsimd::batch<double> drx2 = xsimd::mul(drx, drx);
       const xsimd::batch<double> dry2 = xsimd::mul(dry, dry);
       const xsimd::batch<double> drz2 = xsimd::mul(drz, drz);

       const xsimd::batch<double> dr2PART = xsimd::add(drx2, dry2);
       const xsimd::batch<double> dr2 = xsimd::add(dr2PART, drz2);

       // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
       // signaling = throw error if NaN is encountered
       // dr2 <= _cutoffsquare ? 0xFFFFFFFFFFFFFFFF : 0
       //const xsimd::batch<double> cutoffMask = _mm256_cmp_pd(dr2, _cutoffsquare, _CMP_LE_OS);
       const xsimd::batch_bool<double> cutoffMask = xsimd::le(dr2, _cutoffsquare);

       const xsimd::batch<double> ownedStateJ = remainderIsMasked
                                       ? xsimd::bitwise_and(
                                                xsimd::bitwise_cast<double>(ownedStatePtr2[j]), _masks[rest - 1])
                                       : xsimd::bitwise_cast<double>(ownedStatePtr2[j]);
       // This requires that dummy is the first entry in OwnershipState!
       const xsimd::batch_bool<double> dummyMask = xsimd::neq(ownedStateJ, _zero);
       const xsimd::batch_bool<double> cutoffDummyMask = xsimd::bitwise_and(cutoffMask, dummyMask);

       // if everything is masked away return from this function.
       if (none(cutoffDummyMask)) {
         return;
       }

       const xsimd::batch<double> invdr2 = xsimd::div(_one, dr2);
       const xsimd::batch<double> lj2 = xsimd::mul(sigmaSquares, invdr2);
       const xsimd::batch<double> lj4 = xsimd::mul(lj2, lj2);
       const xsimd::batch<double> lj6 = xsimd::mul(lj2, lj4);
       const xsimd::batch<double> lj12 = xsimd::mul(lj6, lj6);
       const xsimd::batch<double> lj12m6 = xsimd::sub(lj12, lj6);
       const xsimd::batch<double> lj12m6alj12 = xsimd::add(lj12m6, lj12);
       const xsimd::batch<double> lj12m6alj12e = xsimd::mul(lj12m6alj12, epsilon24s);
       const xsimd::batch<double> fac = xsimd::mul(lj12m6alj12e, invdr2);

       const xsimd::batch<double> facMasked =
           remainderIsMasked ? xsimd::bitwise_and(fac, xsimd::select(cutoffDummyMask, _masks[rest - 1], _zero))
                             : xsimd::select(cutoffDummyMask, fac, _zero);

       const xsimd::batch<double> fx = xsimd::mul(drx, facMasked);
       const xsimd::batch<double> fy = xsimd::mul(dry, facMasked);
       const xsimd::batch<double> fz = xsimd::mul(drz, facMasked);

       fxacc = xsimd::add(fxacc, fx);
       fyacc = xsimd::add(fyacc, fy);
       fzacc = xsimd::add(fzacc, fz);

       // if newton 3 is used subtract fD from particle j
       if constexpr (newton3) {
         const xsimd::batch<double> fx2 =
             remainderIsMasked ? xsimd::bitwise_and(xsimd::load(&fx2ptr[j]), _masks[rest - 1]) : xsimd::load(&fx2ptr[j]);
         const xsimd::batch<double> fy2 =
             remainderIsMasked ? xsimd::bitwise_and(xsimd::load(&fy2ptr[j]), _masks[rest - 1]) : xsimd::load(&fy2ptr[j]);
         const xsimd::batch<double> fz2 =
             remainderIsMasked ? xsimd::bitwise_and(xsimd::load(&fz2ptr[j]), _masks[rest - 1]) : xsimd::load(&fz2ptr[j]);

         const xsimd::batch<double> fx2new = xsimd::sub(fx2, fx);
         const xsimd::batch<double> fy2new = xsimd::sub(fy2, fy);
         const xsimd::batch<double> fz2new = xsimd::sub(fz2, fz);


         ;
         remainderIsMasked ? xsimd::store(&fx2ptr[j], xsimd::bitwise_and(_masks[rest - 1], fx2new))
                           : xsimd::store(&fx2ptr[j], fx2new);
         remainderIsMasked ? xsimd::store(&fy2ptr[j], xsimd::bitwise_and(_masks[rest - 1], fy2new))
                           : xsimd::store(&fy2ptr[j], fy2new);
         remainderIsMasked ? xsimd::store(&fz2ptr[j], xsimd::bitwise_and(_masks[rest - 1], fz2new))
                           : xsimd::store(&fz2ptr[j], fz2new);
       }

       if constexpr (calculateGlobals) {
         //TODO
       }


     }

    public:
     // clang-format off
  /**
   * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
     // clang-format on
     inline void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst,
                                  const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                                  bool newton3) final {
       if (soa.getNumberOfParticles() == 0 or neighborList.empty()) return;
       if (newton3) {
         SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
       } else {
         SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
       }
     }

    private:
     template <bool newton3>
     inline void SoAFunctorVerletImpl(SoAView<SoAArraysType> soa, const size_t indexFirst,
                                      const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
       const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();
       if (ownedStatePtr[indexFirst] == OwnershipState::dummy) {
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
       xsimd::batch<double> virialSumX{0};
       xsimd::batch<double> virialSumY{0};
       xsimd::batch<double> virialSumZ{0};
       xsimd::batch<double> upotSum{0};
       xsimd::batch<double> fxacc{0};
        xsimd::batch<double> fyacc{0};
         xsimd::batch<double> fzacc{0};

       // broadcast particle 1
       const xsimd::batch<double> x1 = xsimd::broadcast(&xptr[indexFirst]);
       const xsimd::batch<double> y1 = xsimd::broadcast(&yptr[indexFirst]);
       const xsimd::batch<double> z1 = xsimd::broadcast(&zptr[indexFirst]);
       // ownedStatePtr contains int64_t, so we broadcast these to make an __m256i.
       // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
       xsimd::batch<int64_t> ownedStateI{static_cast<int64_t>(ownedStatePtr[indexFirst])};

       alignas(64) std::array<double, vecLength> x2tmp{};
       alignas(64) std::array<double, vecLength> y2tmp{};
       alignas(64) std::array<double, vecLength> z2tmp{};
       alignas(64) std::array<double, vecLength> fx2tmp{};
       alignas(64) std::array<double, vecLength> fy2tmp{};
       alignas(64) std::array<double, vecLength> fz2tmp{};
       alignas(64) std::array<size_t, vecLength> typeID2tmp{};
       alignas(64) std::array<autopas::OwnershipState, vecLength> ownedStates2tmp{};

       // load 4 neighbors
       size_t j = 0;
       // Loop over all neighbors as long as we can fill full vectors
       // (until `neighborList.size() - neighborList.size() % vecLength`)
       //
       // If b is a power of 2 the following holds:
       // a & ~(b - 1) == a - (a mod b)
       for (; j < (neighborList.size() & ~(vecLength - 1)); j += vecLength) {
         // AVX2 variant:
         // create buffer for 4 interaction particles
         // and fill buffers via gathering
         //      const __m256d x2tmp = _mm256_i64gather_pd(&xptr[j], _vindex, 1);
         //      const __m256d y2tmp = _mm256_i64gather_pd(&yptr[j], _vindex, 1);
         //      const __m256d z2tmp = _mm256_i64gather_pd(&zptr[j], _vindex, 1);
         //      const __m256d fx2tmp = _mm256_i64gather_pd(&fxptr[j], _vindex, 1);
         //      const __m256d fy2tmp = _mm256_i64gather_pd(&fyptr[j], _vindex, 1);
         //      const __m256d fz2tmp = _mm256_i64gather_pd(&fzptr[j], _vindex, 1);
         //      const __m256i typeID2tmp = _mm256_i64gather_epi64(&typeIDptr[j], _vindex, 1);

         for (size_t vecIndex = 0; vecIndex < vecLength; ++vecIndex) {
           x2tmp[vecIndex] = xptr[neighborList[j + vecIndex]];
           y2tmp[vecIndex] = yptr[neighborList[j + vecIndex]];
           z2tmp[vecIndex] = zptr[neighborList[j + vecIndex]];
           if constexpr (newton3) {
             fx2tmp[vecIndex] = fxptr[neighborList[j + vecIndex]];
             fy2tmp[vecIndex] = fyptr[neighborList[j + vecIndex]];
             fz2tmp[vecIndex] = fzptr[neighborList[j + vecIndex]];
           }
           typeID2tmp[vecIndex] = typeIDptr[neighborList[j + vecIndex]];
           ownedStates2tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
         }

         SoAKernel<newton3, false>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2tmp.data()), x1, y1, z1,
                                   x2tmp.data(), y2tmp.data(), z2tmp.data(), fx2tmp.data(), fy2tmp.data(), fz2tmp.data(),
                                   &typeIDptr[indexFirst], typeID2tmp.data(), fxacc, fyacc, fzacc, &virialSumX,
                                   &virialSumY, &virialSumZ, &upotSum, 0);

         if constexpr (newton3) {
           for (size_t vecIndex = 0; vecIndex < vecLength; ++vecIndex) {
             fxptr[neighborList[j + vecIndex]] = fx2tmp[vecIndex];
             fyptr[neighborList[j + vecIndex]] = fy2tmp[vecIndex];
             fzptr[neighborList[j + vecIndex]] = fz2tmp[vecIndex];
           }
         }
       }
       // Remainder loop
       // If b is a power of 2 the following holds:
       // a & (b - 1) == a mod b
       const auto rest = static_cast<int>(neighborList.size() & (vecLength - 1));
       if (rest > 0) {
         // AVX2 variant:
         // create buffer for 4 interaction particles
         // and fill buffers via gathering
         //      TODO: use masked load because there will not be enough data left for the whole gather
         //      const __m256d x2tmp = _mm256_i64gather_pd(&xptr[j], _vindex, 1);
         //      const __m256d y2tmp = _mm256_i64gather_pd(&yptr[j], _vindex, 1);
         //      const __m256d z2tmp = _mm256_i64gather_pd(&zptr[j], _vindex, 1);
         //      const __m256d fx2tmp = _mm256_i64gather_pd(&fxptr[j], _vindex, 1);
         //      const __m256d fy2tmp = _mm256_i64gather_pd(&fyptr[j], _vindex, 1);
         //      const __m256d fz2tmp = _mm256_i64gather_pd(&fzptr[j], _vindex, 1);
         //      const __m256d typeID2tmp = _mm256_i64gather_pd(&typeIDptr[j], _vindex, 1);

         for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
           x2tmp[vecIndex] = xptr[neighborList[j + vecIndex]];
           y2tmp[vecIndex] = yptr[neighborList[j + vecIndex]];
           z2tmp[vecIndex] = zptr[neighborList[j + vecIndex]];
           // if newton3 is used we need to load f of particle j so the kernel can update it too
           if constexpr (newton3) {
             fx2tmp[vecIndex] = fxptr[neighborList[j + vecIndex]];
             fy2tmp[vecIndex] = fyptr[neighborList[j + vecIndex]];
             fz2tmp[vecIndex] = fzptr[neighborList[j + vecIndex]];
           }
           typeID2tmp[vecIndex] = typeIDptr[neighborList[j + vecIndex]];
           ownedStates2tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
         }

         SoAKernel<newton3, true>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2tmp.data()), x1, y1, z1,
                                  x2tmp.data(), y2tmp.data(), z2tmp.data(), fx2tmp.data(), fy2tmp.data(), fz2tmp.data(),
                                  &typeIDptr[indexFirst], typeID2tmp.data(), fxacc, fyacc, fzacc, &virialSumX, &virialSumY,
                                  &virialSumZ, &upotSum, rest);

         if constexpr (newton3) {
           for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
             fxptr[neighborList[j + vecIndex]] = fx2tmp[vecIndex];
             fyptr[neighborList[j + vecIndex]] = fy2tmp[vecIndex];
             fzptr[neighborList[j + vecIndex]] = fz2tmp[vecIndex];
           }
         }
       }

       fxptr[indexFirst] += xsimd::reduce_add(fxacc);
       fyptr[indexFirst] += xsimd::reduce_add(fyacc);
       fzptr[indexFirst] += xsimd::reduce_add(fzacc);

       if constexpr (calculateGlobals) {
         //TODO
       }

     }

    private:
     /**
   * Wrapper function for FMA. If FMA is not supported it executes first the multiplication then the addition.
   * @param factorA
   * @param factorB
   * @param summandC
   * @return A * B + C
      */
     inline xsimd::batch<double> wrapperFMA(const xsimd::batch<double> &factorA, const xsimd::batch<double> &factorB, const xsimd::batch<double> &summandC) {
       //TODO: if fma not supported, is it still working with xsimd?
       return xsimd::fma(factorA, factorB, summandC);
     }

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


     const xsimd::batch<double> _zero{0};
     const xsimd::batch<double> _one{1.};
     const xsimd::batch<int> _vindex{0, 1, 3, 4};
     const xsimd::batch<double> _masks[3]{
         xsimd::batch<double>(0, 0, 0, -1),
         xsimd::batch<double>(0, 0, -1, -1),
         xsimd::batch<double>(0, -1, -1, -1),
     };
     const __m256i _ownedStateDummyMM256i{0x0};
     const __m256i _ownedStateOwnedMM256i{_mm256_set1_epi64x(static_cast<int64_t>(OwnershipState::owned))};
     const xsimd::batch<double> _cutoffsquare{};
     xsimd::batch<double> _shift6{0};
     xsimd::batch<double> _epsilon24{};
     xsimd::batch<double> _sigmaSquare{};

     const double _cutoffsquareAoS = 0;
     double _epsilon24AoS, _sigmaSquareAoS, _shift6AoS = 0;

     ParticlePropertiesLibrary<double, size_t> *_PPLibrary = nullptr;

     // sum of the potential energy, only calculated if calculateGlobals is true
     double _upotSum;

     // sum of the virial, only calculated if calculateGlobals is true
     std::array<double, 3> _virialSum;

     // thread buffer for aos
     std::vector<AoSThreadData> _aosThreadData;

     // defines whether or whether not the global values are already preprocessed
     bool _postProcessed;

     // number of double values that fit into a vector register.
     // MUST be power of 2 because some optimizations make this assumption
     constexpr static size_t vecLength = xsimd::batch<double>::size;
};
}

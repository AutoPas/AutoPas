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

         xsimd::batch<double> x1 = xsimd::broadcast(xptr[i]);
         xsimd::batch<double> y1 = xsimd::broadcast(yptr[i]);
         xsimd::batch<double> z1 = xsimd::broadcast(zptr[i]);

         size_t j = 0;

         // floor soa numParticles to multiple of vecLength
         // If b is a power of 2 the following holds:
         // a & ~(b -1) == a - (a mod b)
         for (; j < (i & ~(vecLength - 1)); j += vecLength) {
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

       if constexpr (calculateGlobals) {
         const int threadnum = autopas_get_thread_num();

         double globals[] = {
             xsimd::reduce_add(virialSumX),
             xsimd::reduce_add(virialSumY),
             xsimd::reduce_add(virialSumZ),
             xsimd::reduce_add(upotSum)
         };

         double factor = 1.;
         // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
         // false, since for newton3 disabled we divide by two later on.
         factor *= newton3 ? .5 : 1.;
         // In case we have a non-cell-wise owned state, we have multiplied everything by two, so we divide it by 2 again.
         _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
         _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
         _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
         _aosThreadData[threadnum].upotSum += globals[3] * factor;
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

         const xsimd::batch<double> x1 = xsimd::broadcast(x1ptr[i]);
         const xsimd::batch<double> y1 = xsimd::broadcast(y1ptr[i]);
         const xsimd::batch<double> z1 = xsimd::broadcast(z1ptr[i]);

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
         const int threadnum = autopas_get_thread_num();

         double globals[] = {
             xsimd::reduce_add(virialSumX),
             xsimd::reduce_add(virialSumY),
             xsimd::reduce_add(virialSumZ),
             xsimd::reduce_add(upotSum)
         };

         // we have duplicated calculations, i.e., we calculate interactions multiple times, so we have to take care
         // that we do not add the energy multiple times!
         double energyfactor = 1.;
         if constexpr (newton3) {
           energyfactor *= 0.5;  // we count the energies partly to one of the two cells!
         }

         _aosThreadData[threadnum].virialSum[0] += globals[0] * energyfactor;
         _aosThreadData[threadnum].virialSum[1] += globals[1] * energyfactor;
         _aosThreadData[threadnum].virialSum[2] += globals[2] * energyfactor;
         _aosThreadData[threadnum].upotSum += globals[3] * energyfactor;
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
     inline void SoAKernel(const size_t j, const xsimd::batch<int64_t> ownedStateI, const int64_t *const __restrict ownedStatePtr2,
                           const xsimd::batch<double> &x1, const xsimd::batch<double> &y1, const xsimd::batch<double> &z1, const double *const __restrict x2ptr,
                           const double *const __restrict y2ptr, const double *const __restrict z2ptr,
                           double *const __restrict fx2ptr, double *const __restrict fy2ptr,
                           double *const __restrict fz2ptr, const size_t *const typeID1ptr, const size_t *const typeID2ptr,
                           xsimd::batch<double> &fxacc, xsimd::batch<double> &fyacc, xsimd::batch<double> &fzacc, xsimd::batch<double> *virialSumX, xsimd::batch<double> *virialSumY,
                           xsimd::batch<double> *virialSumZ, xsimd::batch<double> *upotSum, const unsigned int rest = 0) {
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
       xsimd::batch<double> x2;
       xsimd::batch<double> y2;
       xsimd::batch<double> z2;
       // load only masked values
       //TODO: doesnt work with other sized batches than 4
       if(remainderIsMasked) {
         x2 = {
           _masks[rest-1].get(0) ? x2ptr[j] : 0,
           _masks[rest-1].get(1) ? x2ptr[j + 1] : 0,
           _masks[rest-1].get(2) ? x2ptr[j + 2] : 0,
           _masks[rest-1].get(3) ? x2ptr[j + 3] : 0,
         };
         y2 = {
             _masks[rest-1].get(0) ? y2ptr[j] : 0,
             _masks[rest-1].get(1) ? y2ptr[j + 1] : 0,
             _masks[rest-1].get(2) ? y2ptr[j + 2] : 0,
             _masks[rest-1].get(3) ? y2ptr[j + 3] : 0,
         };
         z2 = {
             _masks[rest-1].get(0) ? z2ptr[j] : 0,
             _masks[rest-1].get(1) ? z2ptr[j + 1] : 0,
             _masks[rest-1].get(2) ? z2ptr[j + 2] : 0,
             _masks[rest-1].get(3) ? z2ptr[j + 3] : 0,
         };
       } else {
         x2 = xsimd::load_unaligned(&x2ptr[j]);
         y2 = xsimd::load_unaligned(&y2ptr[j]);
         z2 = xsimd::load_unaligned(&z2ptr[j]);
       }
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

       const xsimd::batch_bool<double> cutoffMask = xsimd::le(dr2, _cutoffsquare);
       const xsimd::batch<int64_t> _zeroI = xsimd::to_int(_zero);
       const xsimd::batch<int64_t> ownedStateJ = remainderIsMasked
                                                     ? xsimd::select(xsimd::batch_bool_cast<int64_t>(_masks[rest - 1]), 
                                                               xsimd::load_unaligned(&ownedStatePtr2[j]), _zeroI)
                                                     : xsimd::load_unaligned(&ownedStatePtr2[j]);
       const xsimd::batch_bool<double> dummyMask = xsimd::batch_bool_cast<double>(xsimd::neq(ownedStateJ, _zeroI));
       
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
           remainderIsMasked ? xsimd::select(xsimd::bitwise_and(cutoffDummyMask, _masks[rest - 1]), fac, _zero)
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
             remainderIsMasked ? xsimd::select(_masks[rest - 1], xsimd::load_unaligned(&fx2ptr[j]), _zero) : xsimd::load_unaligned(&fx2ptr[j]);
         const xsimd::batch<double> fy2 =
             remainderIsMasked ? xsimd::select(_masks[rest - 1], xsimd::load_unaligned(&fy2ptr[j]), _zero) : xsimd::load_unaligned(&fy2ptr[j]);
         const xsimd::batch<double> fz2 =
             remainderIsMasked ? xsimd::select(_masks[rest - 1], xsimd::load_unaligned(&fz2ptr[j]), _zero) : xsimd::load_unaligned(&fz2ptr[j]);

         const xsimd::batch<double> fx2new = xsimd::sub(fx2, fx);
         const xsimd::batch<double> fy2new = xsimd::sub(fy2, fy);
         const xsimd::batch<double> fz2new = xsimd::sub(fz2, fz);


        // store only masked values
         if(remainderIsMasked) {
           for(int i = 0; _masks[rest - 1].get(i); ++i) {
             fx2ptr[j + i] = fx2new.get(i);
             fy2ptr[j + i] = fy2new.get(i);
             fz2ptr[j + i] = fz2new.get(i);
           }
         } else {
           xsimd::store_unaligned(&fx2ptr[j], fx2new);
           xsimd::store_unaligned(&fy2ptr[j], fy2new);
           xsimd::store_unaligned(&fz2ptr[j], fz2new);
         }
       }

       if constexpr (calculateGlobals) {
         const xsimd::batch<double> virialX = xsimd::mul(fx, drx);
         const xsimd::batch<double> virialY = xsimd::mul(fy, dry);
         const xsimd::batch<double> virialZ = xsimd::mul(fz, drz);

         const xsimd::batch<double> updot = wrapperFMA(epsilon24s, lj12m6, shift6s);

         const __m256d upotMasked =
             remainderIsMasked ? xsimd::select(xsimd::bitwise_and(cutoffDummyMask, _masks[rest - 1]), updot, _zero)
                               : xsimd::select(cutoffDummyMask, updot, _zero);
         xsimd::batch_bool<double> ownedMaskI = xsimd::eq(to_float(ownedStateI), to_float(_ownedStateOwnedMM256i));
         xsimd::batch<double> energyFactor = xsimd::select(ownedMaskI, _one, _zero);
         if constexpr (newton3) {
           xsimd::batch_bool<double> ownedMaskJ =
               xsimd::eq(to_float(ownedStateJ), to_float(_ownedStateOwnedMM256i));
           energyFactor = xsimd::add(energyFactor, xsimd::select(ownedMaskJ, _one, _zero));
         }
         *upotSum = wrapperFMA(energyFactor, upotMasked, *upotSum);
         *virialSumX = wrapperFMA(energyFactor, virialX, *virialSumX);
         *virialSumY = wrapperFMA(energyFactor, virialY, *virialSumY);
         *virialSumZ = wrapperFMA(energyFactor, virialZ, *virialSumZ);
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
       const xsimd::batch<double> x1 = xsimd::broadcast(xptr[indexFirst]);
       const xsimd::batch<double> y1 = xsimd::broadcast(yptr[indexFirst]);
       const xsimd::batch<double> z1 = xsimd::broadcast(zptr[indexFirst]);
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
         const int threadnum = autopas_get_thread_num();

         double globals[] = {
             xsimd::reduce_add(virialSumX),
             xsimd::reduce_add(virialSumY),
             xsimd::reduce_add(virialSumZ),
             xsimd::reduce_add(upotSum)
         };

         double factor = 1.;
         // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
         // false, since for newton3 disabled we divide by two later on.
         factor *= newton3 ? .5 : 1.;
         // In case we have a non-cell-wise owned state, we have multiplied everything by two, so we divide it by 2 again.
         _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
         _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
         _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
         _aosThreadData[threadnum].upotSum += globals[3] * factor;
       }

     }
    public:
     /**
   * @copydoc Functor::getNeededAttr()
      */
     constexpr static auto getNeededAttr() {
       return std::array<typename Particle::AttributeNames, 9>{
           Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
           Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
           Particle::AttributeNames::forceZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
     }

     /**
   * @copydoc Functor::getNeededAttr(std::false_type)
      */
     constexpr static auto getNeededAttr(std::false_type) {
       return std::array<typename Particle::AttributeNames, 6>{
           Particle::AttributeNames::id,   Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
           Particle::AttributeNames::posZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
     }

     /**
   * @copydoc Functor::getComputedAttr()
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
     void initTraversal() final {
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
     void endTraversal(bool newton3) final {
       using namespace autopas::utils::ArrayMath::literals;

       if (_postProcessed) {
         throw utils::ExceptionHandler::AutoPasException(
             "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
       }

       if (calculateGlobals) {
         for (size_t i = 0; i < _aosThreadData.size(); ++i) {
           _upotSum += _aosThreadData[i].upotSum;
           _virialSum += _aosThreadData[i].virialSum;
         }
         if (not newton3) {
           // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
           // here.
           _upotSum *= 0.5;
           _virialSum *= 0.5;
         }
         // we have always calculated 6*upot, so we divide by 6 here!
         _upotSum /= 6.;
         _postProcessed = true;
       }
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
     /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param epsilon24
   * @param sigmaSquare
      */
     void setParticleProperties(double epsilon24, double sigmaSquare) {
       _epsilon24 = epsilon24;
       _sigmaSquare = sigmaSquare;
       if constexpr (applyShift) {
         _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffsquare.get(0));
       } else {
         _shift6 = 0;
       }
        _epsilon24AoS = epsilon24;
        _sigmaSquareAoS = sigmaSquare;
        if constexpr (applyShift) {
       _shift6AoS = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffsquareAoS);
     } else {
       _shift6AoS = 0.;
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
     const xsimd::batch<int64_t> _vindex{0, 1, 3, 4};
     const xsimd::batch_bool<double> _masks[3] {
         xsimd::batch_bool<double>(true, false, false, false),
         xsimd::batch_bool<double>(true, true, false, false),
         xsimd::batch_bool<double>(true, true, true, false)
     };
     const xsimd::batch<int64_t> _ownedStateDummyMM256i{0x0};
     const xsimd::batch<int64_t> _ownedStateOwnedMM256i{static_cast<int64_t>(OwnershipState::owned)};
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

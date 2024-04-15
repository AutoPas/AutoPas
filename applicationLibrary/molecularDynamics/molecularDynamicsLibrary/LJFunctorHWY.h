// Created by Luis Gall on 27.03.2024

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "molecularDynamicsLibrary/LJFunctorHWY.h"
// #include <hwy/foreach_target.h>
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();
namespace mdLib {
namespace HWY_NAMESPACE {

    namespace highway = hwy::HWY_NAMESPACE;

    // architecture specific information
    const highway::ScalableTag<double> tag_double;
    const highway::ScalableTag<int64_t> tag_long;
    const size_t _vecLengthDouble {highway::Lanes(tag_double)};
    using VectorDouble = decltype(highway::Zero(tag_double));
    using VectorLong = decltype(highway::Zero(tag_long));
    
    // template <T>
    // using Vector<T> = decltype(highway::Zero(highway::ScalableTag<T>));

    template <bool applyShift = false, bool useMixing = false,
        autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
        bool relevantForTuning = true>

    class LJFunctorHWY
        : public autopas::Functor<mdLib::MoleculeLJ_NoPPL,
                         LJFunctorHWY<applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
        
        using molecule = typename mdLib::MoleculeLJ_NoPPL;

        using SoAArraysType = typename molecule::SoAArraysType;

        public:
            LJFunctorHWY() = delete;

        private:
            explicit LJFunctorHWY(double cutoff, void*)
                : autopas::Functor<molecule,
                LJFunctorHWY<applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(cutoff),
                _cutoffSquared{highway::Set(tag_double, cutoff*cutoff)},
                _cutoffSquareAoS{cutoff * cutoff},
                _uPotSum{0.},
                _virialSum{0.,0.,0.},
                _aosThreadData(),
                _postProcessed{false} {
                    if (calculateGlobals) {
                        _aosThreadData.resize(autopas::autopas_get_max_threads());
                    }

                    // AutoPasLog(INFO, "Highway Wrapper initialized with a register size of ({}) with architecture {}.", _vecLengthDouble, hwy::TargetName(HWY_TARGET));
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
            explicit LJFunctorHWY(double cutoff) : LJFunctorHWY(cutoff, nullptr) {
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
            explicit LJFunctorHWY(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
                : LJFunctorHWY(cutoff, nullptr) {
                    static_assert(useMixing,
                        "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                        "or set mixing to true.");
                    _PPLibrary = &particlePropertiesLibrary;
                }

            bool isRelevantForTuning() final { return relevantForTuning; }

            bool allowsNewton3() final { return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both; }

            bool allowsNonNewton3() final {
                return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
            }


            inline void AoSFunctor(molecule &i, molecule &j, bool newton3) final {
                using namespace autopas::utils::ArrayMath::literals;
                if (i.isDummy() or j.isDummy()) {
                    return;
                }
                auto sigmasquare = _sigmaSquareAoS;
                auto epsilon24 = _epsilon24AoS;
                auto shift6 = _shift6AoS;
                if constexpr (useMixing) {
                    ;
                }
                auto dr = i.getR() - j.getR();
                double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

                if (dr2 > _cutoffSquareAoS) {
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

                    const int threadnum = autopas::autopas_get_thread_num();
                    // for non-newton3 the division is in the post-processing step.
                    if (newton3) {
                    upot *= 0.5;
                    virial *= (double)0.5;
                    }
                    if (i.isOwned()) {
                    _aosThreadData[threadnum].uPotSum += upot;
                    _aosThreadData[threadnum].virialSum += virial;
                    }
                    // for non-newton3 the second particle will be considered in a separate calculation
                    if (newton3 and j.isOwned()) {
                    _aosThreadData[threadnum].uPotSum += upot;
                    _aosThreadData[threadnum].virialSum += virial;
                    }
                }
                }

                /**
                 * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
                 * This functor will always do a newton3 like traversal of the soa.
                 * However, it still needs to know about newton3 to correctly add up the global values.
                */
                inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
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
                inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, const bool newton3) final {
                    if (newton3) {
                        SoAFunctorPairImpl<true>(soa1, soa2);
                    } else {
                        SoAFunctorPairImpl<false>(soa1, soa2);
                    }
                }

            private:
                template <bool newton3>
                inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {

                    if (soa.size() == 0) return;

                    // obtain iterators for the various values
                    const auto *const __restrict xPtr = soa.template begin<molecule::AttributeNames::posX>();
                    const auto *const __restrict yPtr = soa.template begin<molecule::AttributeNames::posY>();
                    const auto *const __restrict zPtr = soa.template begin<molecule::AttributeNames::posZ>();

                    const auto *const __restrict ownedStatePtr = soa.template begin<molecule::AttributeNames::ownershipState>();

                    auto *const __restrict fxPtr = soa.template begin<molecule::AttributeNames::forceX>();
                    auto *const __restrict fyPtr = soa.template begin<molecule::AttributeNames::forceY>();
                    auto *const __restrict fzPtr = soa.template begin<molecule::AttributeNames::forceZ>();

                    // const auto *const __restrict typeIDptr = soa.template begin<molecule::AttributeNames::typeId>();

                    // initialize vector variables
                    auto virialSumX = highway::Zero(tag_double);
                    auto virialSumY = highway::Zero(tag_double);
                    auto virialSumZ = highway::Zero(tag_double);
                    auto uPotSum = highway::Zero(tag_double);

                    // loop over list for first time
                    for (size_t i = soa.size() - 1; (long)i >= 0; --i) {
                        if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
                            continue;
                        }
                        // static_assert(std::is_same_v<std::underlying_type<autopas::OwnershipState>, int64_t>,
                        //      "OwnershipStates' underlying type should be int64_t");

                        VectorLong ownedStateI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr[i]));
                        VectorDouble x1 = highway::Set(tag_double, xPtr[i]);
                        VectorDouble y1 = highway::Set(tag_double, yPtr[i]);
                        VectorDouble z1 = highway::Set(tag_double, zPtr[i]);

                        VectorDouble fxAcc = highway::Zero(tag_double);
                        VectorDouble fyAcc = highway::Zero(tag_double);
                        VectorDouble fzAcc = highway::Zero(tag_double);

                        size_t j = 0;

                        // floor soa numParticles to multiple of vecLength
                        // If b is a power of 2 the following holds:
                        // a & ~(b -1) == a - (a mod b)
                        for (; j < (i & ~(_vecLengthDouble - 1)); j += _vecLengthDouble) {
                            SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xPtr,
                                                yPtr, zPtr, fxPtr, fyPtr, fzPtr, fxAcc, fyAcc, fzAcc,
                                                virialSumX, virialSumY, virialSumZ, uPotSum, 0);

                        }

                        // If b is a power of 2 the following holds:
                        // a & (b -1) == a mod b
                        const int rest = (int)(i & (_vecLengthDouble - 1));
                        if (rest > 0) {

                            SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xPtr,
                                                yPtr, zPtr, fxPtr, fyPtr, fzPtr, fxAcc, fyAcc, fzAcc,
                                                virialSumX, virialSumY, virialSumZ, uPotSum, rest);
                        }

                        double sumFx = highway::ReduceSum(tag_double, fxAcc);
                        double sumFy = highway::ReduceSum(tag_double, fyAcc);
                        double sumFz = highway::ReduceSum(tag_double, fzAcc);

                        fxPtr[i] += sumFx;
                        fyPtr[i] += sumFy;
                        fzPtr[i] += sumFz;
                    }

                    if constexpr (calculateGlobals) {
                        const int threadnum = autopas::autopas_get_thread_num();

                        double globals[] = {
                            highway::ReduceSum(tag_double, virialSumX),
                            highway::ReduceSum(tag_double, virialSumY),
                            highway::ReduceSum(tag_double, virialSumZ),
                            highway::ReduceSum(tag_double, uPotSum)
                        };

                        double factor = 1.;
                        factor *= newton3 ? 0.5 : 1.;
                        _aosThreadData[threadnum].viralSum[0] += globals[0] * factor;
                        _aosThreadData[threadnum].viralSum[1] += globals[1] * factor;
                        _aosThreadData[threadnum].viralSum[2] += globals[2] * factor;
                        _aosThreadData[threadnum].uPotSum += globals[3] * factor;
                    }
                }

                template <bool newton3>
                inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {

                    const auto *const __restrict x1Ptr = soa1.template begin<molecule::AttributeNames::posX>();
                    const auto *const __restrict y1Ptr = soa1.template begin<molecule::AttributeNames::posY>();
                    const auto *const __restrict z1Ptr = soa1.template begin<molecule::AttributeNames::posZ>();
                    const auto *const __restrict x2Ptr = soa2.template begin<molecule::AttributeNames::posX>();
                    const auto *const __restrict y2Ptr = soa2.template begin<molecule::AttributeNames::posY>();
                    const auto *const __restrict z2Ptr = soa2.template begin<molecule::AttributeNames::posZ>();

                    const auto *const __restrict ownedStatePtr1 = soa1.template begin<molecule::AttributeNames::ownershipState>();
                    const auto *const __restrict ownedStatePtr2 = soa2.template begin<molecule::AttributeNames::ownershipState>();

                    auto *const __restrict fx1Ptr = soa1.template begin<molecule::AttributeNames::forceX>();
                    auto *const __restrict fy1Ptr = soa1.template begin<molecule::AttributeNames::forceY>();
                    auto *const __restrict fz1Ptr = soa1.template begin<molecule::AttributeNames::forceZ>();
                    auto *const __restrict fx2Ptr = soa2.template begin<molecule::AttributeNames::forceX>();
                    auto *const __restrict fy2Ptr = soa2.template begin<molecule::AttributeNames::forceY>();
                    auto *const __restrict fz2Ptr = soa2.template begin<molecule::AttributeNames::forceZ>();

                    // const auto *const __restrict typeID1Ptr = soa1.template begin<molecule::AttributeNames::typeId>();
                    // const auto *const __restrict typeID2Ptr = soa2.template begin<molecule::AttributeNames::typeId>();

                    VectorDouble virialSumX = highway::Zero(tag_double);
                    VectorDouble virialSumY = highway::Zero(tag_double);
                    VectorDouble virialSumZ = highway::Zero(tag_double);
                    VectorDouble uPotSum = highway::Zero(tag_double);

                    for (size_t i = 0; i < soa1.size(); ++i) {
                        if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
                            continue;
                        }

                        VectorDouble fxAcc = highway::Zero(tag_double);
                        VectorDouble fyAcc = highway::Zero(tag_double);
                        VectorDouble fzAcc = highway::Zero(tag_double);

                        VectorLong ownedStateI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr1[i]));

                        const VectorDouble x1 = highway::Set(tag_double, x1Ptr[i]);
                        const VectorDouble y1 = highway::Set(tag_double, y1Ptr[i]);
                        const VectorDouble z1 = highway::Set(tag_double, z1Ptr[i]);

                        size_t j = 0;
                        for (; j < (soa2.size() & ~ (_vecLengthDouble - 1)); j += _vecLengthDouble) {
                            SoAKernel<newton3, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1,
                                    x2Ptr, y2Ptr, z2Ptr, fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc,
                                    virialSumX, virialSumY, virialSumZ, uPotSum, 0);

                        }

                        const int rest = (int) (soa2.size() & (_vecLengthDouble - 1));
                        if (rest > 0) {

                            SoAKernel<newton3, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1,
                                x2Ptr, y2Ptr, z2Ptr, fx2Ptr, fy2Ptr, fz2Ptr, fxAcc, fyAcc, fzAcc,
                                virialSumX, virialSumY, virialSumZ, uPotSum, rest);
                        }

                        fx1Ptr[i] += highway::ReduceSum(tag_double, fxAcc);
                        fy1Ptr[i] += highway::ReduceSum(tag_double, fyAcc);
                        fz1Ptr[i] += highway::ReduceSum(tag_double, fzAcc);
                    }
                    if constexpr (calculateGlobals) {        
                        const int threadnum = autopas::autopas_get_num_threads();

                        double globals [] = {
                            highway::ReduceSum(tag_double, virialSumX),
                            highway::ReduceSum(tag_double, virialSumY),
                            highway::ReduceSum(tag_double, virialSumZ),
                            highway::ReduceSum(tag_double, uPotSum)
                        };

                        double energyfactor = 1.;
                        if constexpr (newton3) {
                            energyfactor *= 0.5;
                        }

                        _aosThreadData[threadnum].virialSum[0] += globals[0] * energyfactor;
                        _aosThreadData[threadnum].virialSum[1] += globals[1] * energyfactor;
                        _aosThreadData[threadnum].virialSum[2] += globals[2] * energyfactor;
                        _aosThreadData[threadnum].upotSum += globals[3] * energyfactor;
                    }                
                }

                template <bool newton3, bool remainderIsMasked>
                inline void SoAKernel(const size_t j, VectorLong &ownedStateI,
                                    const int64_t * __restrict ownedStatePtr2,
                                    const VectorDouble& x1, const VectorDouble& y1, const VectorDouble& z1,
                                    const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
                                    const double *const __restrict z2Ptr, double *const __restrict fx2Ptr,
                                    double *const __restrict fy2Ptr, double *const __restrict fz2Ptr,
                                    VectorDouble& fxAcc, VectorDouble& fyAcc, VectorDouble& fzAcc,
                                    VectorDouble& virialSumX, VectorDouble& virialSumY, VectorDouble& virialSumZ,
                                    VectorDouble& potentialEnergySum, const unsigned int rest = 0) {

                    // masks for loading rest
                    auto restMaskDouble = highway::FirstN(tag_double, rest);
                    auto restMaskLong = highway::FirstN(tag_long, rest);

                    // prepare physical parameters (epsilon,...)
                    VectorDouble epsilon24s;
                    VectorDouble sigmaSquareds;
                    VectorDouble shift6s;

                    if (useMixing) {
                        double epsilon_buf[_vecLengthDouble] = {0};
                        double sigma_buf[_vecLengthDouble] = {0};
                        double shift_buf[_vecLengthDouble] = {0};

                        for(int i=0; (remainderIsMasked ? i < rest : i < _vecLengthDouble); ++i) {
                            epsilon_buf[i] = 0; // _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + i));
                            sigma_buf[i] = 0; // _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + i));
                            if constexpr (applyShift) {
                                shift_buf[i] = 0; // _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + i));
                            }
                        }

                        epsilon24s = highway::LoadU(tag_double, epsilon_buf);
                        sigmaSquareds = highway::LoadU(tag_double, sigma_buf);
                        shift6s = highway::LoadU(tag_double, shift_buf);
                    }
                    else {
                        epsilon24s = _epsilon24;
                        sigmaSquareds = _sigmaSquared;
                        shift6s = _shift6;
                    }

                    // load other interaction partners
                    VectorDouble x2;
                    VectorDouble y2;
                    VectorDouble z2;
                    if constexpr (remainderIsMasked) {

                        x2 = highway::MaskedLoad(restMaskDouble, tag_double, &x2Ptr[j]);
                        y2 = highway::MaskedLoad(restMaskDouble, tag_double, &y2Ptr[j]);
                        z2 = highway::MaskedLoad(restMaskDouble, tag_double, &z2Ptr[j]);

                        return;
                    }
                    else {
                        x2 = highway::LoadU(tag_double, &x2Ptr[j]);
                        y2 = highway::LoadU(tag_double, &y2Ptr[j]);
                        z2 = highway::LoadU(tag_double, &z2Ptr[j]);
                    }

                    // perform distance calculations
                    const VectorDouble drx = highway::Sub(x1, x2);
                    const VectorDouble dry = highway::Sub(y1, y2);
                    const VectorDouble drz = highway::Sub(z1, z2);

                    const VectorDouble drx2 = highway::Mul(drx, drx);
                    const VectorDouble dry2 = highway::Mul(dry, dry);
                    const VectorDouble drz2 = highway::Mul(drz, drz);

                    const VectorDouble dr2_part = highway::Add(drx2, dry2);
                    const VectorDouble dr2 = highway::Add(dr2_part, drz2);

                    // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
                    // signaling = throw error if NaN is encountered
                    // dr2 <= _cutoffsquare ? 0xFFFFFFFFFFFFFFFF : 0
                    const auto cutoffMask = highway::Le(dr2, _cutoffSquared);
                    const VectorLong _zeroI = highway::ConvertTo(tag_long, _zeroDouble);
                    const VectorLong ownedStateJ = remainderIsMasked ? highway::MaskedLoad(restMaskLong, tag_long, &ownedStatePtr2[j])
                                                                    : highway::LoadU(tag_long, &ownedStatePtr2[j]);

                    const auto dummyMask = highway::Ne(highway::ConvertTo(tag_double, ownedStateJ), highway::ConvertTo(tag_double, _zeroLong));
                    const auto cutoffDummyMask = highway::And(cutoffMask, dummyMask);

                    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
                        return;
                    }

                    const VectorDouble invDr2 = highway::Div(_oneDouble, dr2);
                    const VectorDouble lj2 = highway::Mul(sigmaSquareds, invDr2);
                    const VectorDouble lj4 = highway::Mul(lj2, lj2);
                    const VectorDouble lj6 = highway::Mul(lj2, lj4);
                    const VectorDouble lj12 = highway::Mul(lj6, lj6);
                    const VectorDouble lj12m6 = highway::Sub(lj12, lj6);
                    const VectorDouble lj12m6alj12 = highway::Add(lj12m6, lj12);
                    const VectorDouble lj12m6alj12e = highway::Mul(lj12m6alj12, epsilon24s);
                    const VectorDouble fAcc = highway::Mul(lj12m6alj12e, invDr2);

                    const VectorDouble fAccMasked =
                        remainderIsMasked ? highway::IfThenElse(highway::And(cutoffDummyMask, restMaskDouble), fAcc, _zeroDouble)
                                        : highway::IfThenElse(cutoffDummyMask, fAcc, _zeroDouble);

                    const VectorDouble fx = highway::Mul(drx, fAccMasked);
                    const VectorDouble fy = highway::Mul(dry, fAccMasked);
                    const VectorDouble fz = highway::Mul(drz, fAccMasked);

                    fxAcc = highway::Add(fxAcc, fx);
                    fyAcc = highway::Add(fyAcc, fy);
                    fzAcc = highway::Add(fzAcc, fz);

                    if constexpr (newton3) {
                        VectorDouble fx2;
                        VectorDouble fy2;
                        VectorDouble fz2;

                        if constexpr (remainderIsMasked) {

                            fx2 = highway::MaskedLoad(restMaskDouble, tag_double, &fx2Ptr[j]);
                            fy2 = highway::MaskedLoad(restMaskDouble, tag_double, &fy2Ptr[j]);
                            fz2 = highway::MaskedLoad(restMaskDouble, tag_double, &fz2Ptr[j]);

                            return;
                        }
                        else {
                            fx2 = highway::LoadU(tag_double, &fx2Ptr[j]);
                            fy2 = highway::LoadU(tag_double, &fy2Ptr[j]);
                            fz2 = highway::LoadU(tag_double, &fz2Ptr[j]);
                        }

                        const VectorDouble fx2New = highway::Sub(fx2, fx);
                        const VectorDouble fy2New = highway::Sub(fy2, fy);
                        const VectorDouble fz2New = highway::Sub(fz2, fz);

                        if constexpr (remainderIsMasked) {

                            highway::BlendedStore(fx2New, restMaskDouble, tag_double, &fx2Ptr[j]);
                            highway::BlendedStore(fy2New, restMaskDouble, tag_double, &fy2Ptr[j]);
                            highway::BlendedStore(fz2New, restMaskDouble, tag_double, &fz2Ptr[j]);

                            return;
                        }
                        else {
                            highway::StoreU(fx2New, tag_double, &fx2Ptr[j]);
                            highway::StoreU(fy2New, tag_double, &fy2Ptr[j]);
                            highway::StoreU(fz2New, tag_double, &fz2Ptr[j]);
                        }
                    }

                    if constexpr(calculateGlobals) {
                        const VectorDouble virialX = highway::Mul(fx, drx);
                        const VectorDouble virialY = highway::Mul(fy, dry);
                        const VectorDouble virialZ = highway::Mul(fz, drz);

                        const VectorDouble uPDot = highway::MulAdd(epsilon24s, lj12m6, shift6s);

                        const VectorDouble uPotMasked = 
                            remainderIsMasked ? highway::IfThenElse(highway::And(cutoffDummyMask, restMaskDouble) ,uPDot, _zeroDouble)
                                            : highway::IfThenElse(cutoffDummyMask, uPDot, _zeroDouble);

                        auto ownedMaskI = highway::Eq(ownedStateI, _ownedStateOwnedMM256i);
                        VectorDouble energyFactor = highway::ConvertTo(tag_double, highway::IfThenElse(ownedMaskI, _oneLong, _zeroLong));

                        if constexpr (newton3) {
                            auto ownedMaskJ = highway::Eq(ownedStateJ, _ownedStateOwnedMM256i);
                            energyFactor = highway::Add(energyFactor, highway::ConvertTo(tag_double,highway::IfThenElse(ownedMaskJ, _oneLong, _zeroLong)));
                        }
                        energyFactor = highway::MulAdd(energyFactor, uPotMasked, energyFactor);
                        virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
                        virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
                        virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
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
                inline void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
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
                inline void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                    const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {

                        const auto *const __restrict ownedStatePtr = soa.template begin<molecule::AttributeNames::ownershipState>();
                        if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
                            return;
                        }

                        const auto *const __restrict xPtr = soa.template begin<molecule::AttributeNames::posX>();
                        const auto *const __restrict yPtr = soa.template begin<molecule::AttributeNames::posY>();
                        const auto *const __restrict zPtr = soa.template begin<molecule::AttributeNames::posZ>();

                        auto *const __restrict fxPtr = soa.template begin<molecule::AttributeNames::forceX>();
                        auto *const __restrict fyPtr = soa.template begin<molecule::AttributeNames::forceY>();
                        auto *const __restrict fzPtr = soa.template begin<molecule::AttributeNames::forceZ>();

                        // const auto *const __restrict typeIDPtr = soa.template begin<molecule::AttributeNames::typeId>();

                        VectorDouble virialSumX = highway::Zero(tag_double);
                        VectorDouble virialSumY = highway::Zero(tag_double);
                        VectorDouble virialSumZ = highway::Zero(tag_double);
                        VectorDouble uPotSum = highway::Zero(tag_double);
                        VectorDouble fxAcc = highway::Zero(tag_double);
                        VectorDouble fyAcc = highway::Zero(tag_double);
                        VectorDouble fzAcc = highway::Zero(tag_double);

                        const VectorDouble x1 = highway::Set(tag_double, xPtr[indexFirst]);
                        const VectorDouble y1 = highway::Set(tag_double, yPtr[indexFirst]);
                        const VectorDouble z1 = highway::Set(tag_double, zPtr[indexFirst]);

                        VectorLong ownedStateI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr[indexFirst]));

                        alignas(64) std::array<double, _vecLengthDouble> x2Tmp{};
                        alignas(64) std::array<double, _vecLengthDouble> y2Tmp{};
                        alignas(64) std::array<double, _vecLengthDouble> z2Tmp{};
                        alignas(64) std::array<double, _vecLengthDouble> fx2Tmp{};
                        alignas(64) std::array<double, _vecLengthDouble> fy2Tmp{};
                        alignas(64) std::array<double, _vecLengthDouble> fz2Tmp{};
                        alignas(64) std::array<size_t, _vecLengthDouble> typeID2Tmp{};
                        alignas(64) std::array<autopas::OwnershipState, _vecLengthDouble> ownedStates2Tmp{};

                        size_t j = 0;

                        for (; j < (neighborList.size() & ~(_vecLengthDouble - 1)); j += _vecLengthDouble) {

                            for (size_t vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
                                x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
                                y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
                                z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
                                if constexpr (newton3) {
                                    fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
                                    fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
                                    fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
                                }
                                // typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
                                ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
                            }

                            SoAKernel<newton3, false>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                                    x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                                    fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0);

                            if constexpr (newton3) {
                                for (size_t vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
                                    fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
                                    fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
                                    fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
                                }
                            }
                        }

                        const size_t rest = neighborList.size() & (_vecLengthDouble - 1);

                        if (rest > 0) {
                            for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
                                x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
                                y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
                                z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
                                // if newton3 is used we need to load f of particle j so the kernel can update it too
                                if constexpr (newton3) {
                                    fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
                                    fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
                                    fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
                                }
                                // typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
                                ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
                            }

                            SoAKernel<newton3, true>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                                  x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                                  fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, rest);

                            if constexpr (newton3) {
                                for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
                                    fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
                                    fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
                                    fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
                                }
                            }
                        }

                        fxPtr[indexFirst] += highway::ReduceSum(tag_double, fxAcc);
                        fyPtr[indexFirst] += highway::ReduceSum(tag_double, fyAcc);
                        fzPtr[indexFirst] += highway::ReduceSum(tag_double, fzAcc);

                        if constexpr (calculateGlobals) {
                            const int threadnum = autopas::autopas_get_num_threads();

                            double globals[] = {
                                highway::ReduceSum(tag_double, virialSumX),
                                highway::ReduceSum(tag_double, virialSumY),
                                highway::ReduceSum(tag_double, virialSumZ),
                                highway::ReduceSum(tag_double, uPotSum)
                            };

                            double factor = 1.;
                            factor *= newton3 ? .5 : 1.;
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
                return std::array<typename molecule::AttributeNames, 9>{
                    molecule::AttributeNames::id,     molecule::AttributeNames::posX,   molecule::AttributeNames::posY,
                    molecule::AttributeNames::posZ,   molecule::AttributeNames::forceX, molecule::AttributeNames::forceY,
                    molecule::AttributeNames::forceZ, molecule::AttributeNames::ownershipState};
                }

                /**
             * @copydoc Functor::getNeededAttr(std::false_type)
                 */
                constexpr static auto getNeededAttr(std::false_type) {
                return std::array<typename molecule::AttributeNames, 6>{
                    molecule::AttributeNames::id,   molecule::AttributeNames::posX,   molecule::AttributeNames::posY,
                    molecule::AttributeNames::posZ, molecule::AttributeNames::ownershipState};
                }

                /**
             * @copydoc Functor::getComputedAttr()
                 */
                constexpr static auto getComputedAttr() {
                return std::array<typename molecule::AttributeNames, 3>{
                    molecule::AttributeNames::forceX, molecule::AttributeNames::forceY, molecule::AttributeNames::forceZ};
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
                _uPotSum = 0.;
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
                    throw autopas::utils::ExceptionHandler::AutoPasException(
                        "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
                }

                if (calculateGlobals) {
                    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
                    _uPotSum += _aosThreadData[i].uPotSum;
                    _virialSum += _aosThreadData[i].virialSum;
                    }
                    if (not newton3) {
                    // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
                    // here.
                    _uPotSum *= 0.5;
                    _virialSum *= 0.5;
                    }
                    // we have always calculated 6*upot, so we divide by 6 here!
                    _uPotSum /= 6.;
                    _postProcessed = true;
                }
                }

                /**
             * Get the potential Energy
             * @return the potential Energy
                 */
                double getUpot() {
                if (not calculateGlobals) {
                    throw autopas::utils::ExceptionHandler::AutoPasException(
                        "Trying to get upot even though calculateGlobals is false. If you want this functor to calculate global "
                        "values, please specify calculateGlobals to be true.");
                }
                if (not _postProcessed) {
                    throw autopas::utils::ExceptionHandler::AutoPasException("Cannot get upot, because endTraversal was not called.");
                }
                return _uPotSum;
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
                    throw autopas::utils::ExceptionHandler::AutoPasException("Cannot get virial, because endTraversal was not called.");
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
                    _sigmaSquared = sigmaSquare;
                    if constexpr (applyShift) {
                        _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, highway::GetLane(_cutoffSquared));
                    } else {
                        _shift6 = 0;
                    }
                    _epsilon24AoS = epsilon24;
                    _sigmaSquareAoS = sigmaSquare;
                    if constexpr (applyShift) {
                        _shift6AoS = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffSquareAoS);
                    } else {
                        _shift6AoS = 0.;
                    }
                }

            private:

                class AoSThreadData
                {
                    public:
                        AoSThreadData() : virialSum{0.,0.,0.}, uPotSum{0.} {}

                        void setZero() {
                            virialSum = {0.,0.,0.};
                            uPotSum = 0.;
                        }

                        std::array<double, 3> virialSum {};
                        double uPotSum {0};

                    private:
                        double _remainingTo64[4];
                };
                // static_assert(sizeof(AoSThreadData) & 64 == 0, "AoSThreadData has wrong size");

                // helper variables for the LJ-calculation
                const VectorDouble _zeroDouble {highway::Zero(tag_double)};
                const VectorLong _zeroLong {highway::Zero(tag_long)};
                const VectorDouble _oneDouble {highway::Set(tag_double, 1.)};
                const VectorLong _oneLong {highway::Set(tag_long, 1)};
                const VectorLong _ownedStateDummyMM256i{0x0};
                const VectorLong _ownedStateOwnedMM256i{static_cast<int64_t>(autopas::OwnershipState::owned)};
                const VectorDouble _cutoffSquared {};
                VectorDouble _shift6 {0};
                VectorDouble _epsilon24 {0};
                VectorDouble _sigmaSquared {0};

                const double _cutoffSquareAoS {0.};
                double _epsilon24AoS, _sigmaSquareAoS, _shift6AoS = 0.;
                ParticlePropertiesLibrary<double, size_t>* _PPLibrary = nullptr;
                double _uPotSum {0.};
                std::array<double, 3> _virialSum;
                std::vector<AoSThreadData> _aosThreadData;
                bool _postProcessed;
    };
} // Highway
} // mdLib
HWY_AFTER_NAMESPACE();
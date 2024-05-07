// Created by Luis Gall on 27.03.2024

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "VectorizationPatterns.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/ArrayMath.h"

// #undef HWY_TARGET_INCLUDE
// #define HWY_TARGET_INCLUDE "molecularDynamicsLibrary/LJFunctorHWY.h"
// #include <hwy/foreach_target.h>
#include <hwy/highway.h>

// HWY_BEFORE_NAMESPACE();
namespace mdLib {
//namespace HWY_NAMESPACE {

    namespace highway = hwy::HWY_NAMESPACE;

    // architecture specific information
    const highway::ScalableTag<double> tag_double;
    const highway::ScalableTag<int64_t> tag_long;
    const size_t _vecLengthDouble {highway::Lanes(tag_double)};
    using VectorDouble = decltype(highway::Zero(tag_double));
    using VectorLong = decltype(highway::Zero(tag_long));

    using MaskDouble = decltype(highway::FirstN(tag_double, 1));
    using MaskLong = decltype(highway::FirstN(tag_long, 2));

    template <class Particle, bool applyShift = false, bool useMixing = false,
        autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
        bool relevantForTuning = true, VectorizationPattern vecPattern = VectorizationPattern::p1xVec>

    class LJFunctorHWY
        : public autopas::Functor<Particle,
                         LJFunctorHWY<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
        using SoAArraysType = typename Particle::SoAArraysType;

        public:
            LJFunctorHWY() = delete;

        private:
            explicit LJFunctorHWY(double cutoff, void*)
                : autopas::Functor<Particle,
                LJFunctorHWY<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(cutoff),
                _cutoffSquared{highway::Set(tag_double, cutoff*cutoff)},
                _cutoffSquareAoS{cutoff * cutoff},
                _uPotSum{0.},
                _virialSum{0.,0.,0.},
                _aosThreadData(),
                _postProcessed{false} {
                    if (calculateGlobals) {
                        _aosThreadData.resize(autopas::autopas_get_max_threads());
                    }

                    initializeRestMasks();

                    std::cout << hwy::TargetName(HWY_TARGET); // AutoPasLog(INFO, "Highway Wrapper initialized with a register size of ({}) with architecture {}.", _vecLengthDouble, hwy::TargetName(HWY_TARGET));
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


            inline void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
                using namespace autopas::utils::ArrayMath::literals;
                if (i.isDummy() or j.isDummy()) {
                    return;
                }
                auto sigmasquare = _sigmaSquareAoS;
                auto epsilon24 = _epsilon24AoS;
                auto shift6 = _shift6AoS;
                if constexpr (useMixing) {
                    sigmasquare = _PPLibrary->getMixingSigmaSquared(i.getTypeId(), j.getTypeId());
                    epsilon24 = _PPLibrary->getMixing24Epsilon(i.getTypeId(), j.getTypeId());
                    if constexpr (applyShift) {
                    shift6 = _PPLibrary->getMixingShift6(i.getTypeId(), j.getTypeId());
                    }
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

                inline void initializeRestMasks() {

                    // TODO : handle different vectorization patterns

                    for (size_t n = 0; n <_vecLengthDouble-1;++n) {
                        restMasksDouble[n] = highway::FirstN(tag_double, n+1);
                        restMasksLong[n] = highway::FirstN(tag_long, n+1);
                    }
                }
                
                inline void decrementFirstLoop(size_t& i) {

                    if constexpr (vecPattern == VectorizationPattern::p1xVec)
                        --i;
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2)
                        i = i-2;
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2)
                        i = i - (_vecLengthDouble / 2);
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1)
                        i = i - _vecLengthDouble;
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec)
                        i = i - _vecLengthDouble;
                    else
                        throw std::runtime_error("No vectorization pattern matched, error!");
                    
                }

                inline void incrementFirstLoop(size_t& i) {

                    if constexpr (vecPattern == VectorizationPattern::p1xVec)
                        ++i;
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2)
                        i = i + 2;
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2)
                        i = i + (_vecLengthDouble / 2);
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1)
                        i = i + _vecLengthDouble;
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec)
                        i = i + _vecLengthDouble;
                    else
                        throw std::runtime_error("No vectorization pattern matched, error!");
                }

                inline void incrementSecondLoop(size_t& j) {

                    if constexpr (vecPattern == VectorizationPattern::p1xVec)
                        j = j + _vecLengthDouble;
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2)
                        j = j + _vecLengthDouble / 2;
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2)
                        j = j + 2;
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1)
                        ++j;
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec)
                        j = j + _vecLengthDouble;
                    else
                        throw std::runtime_error("No vectorization pattern matched, error!");
                }

                template <bool remainderI, bool remainderJ>
                inline void fillVectorRegisters(VectorDouble& x1, VectorDouble& y1, VectorDouble& z1,
                    VectorDouble& x2, VectorDouble& y2, VectorDouble& z2, VectorLong& ownedI, VectorLong& ownedJ,
                    VectorDouble& epsilon24s, VectorDouble& sigmaSquaredDiv2s, VectorDouble& shift6s,
                    const long unsigned *const __restrict typeID1ptr, const long unsigned *const __restrict typeID2ptr,
                    const double *const __restrict xPtr1, const double *const __restrict yPtr1, const double *const __restrict zPtr1,
                    const double *const __restrict xPtr2, const double *const __restrict yPtr2, const double *const __restrict zPtr2,
                    const int64_t *const __restrict ownedPtr1, const int64_t *const __restrict ownedPtr2, const size_t i, const size_t j,
                    const size_t restI = 0, const size_t restJ = 0) {

                    ownedI = highway::Set(tag_long, ownedPtr1[i]);
                    x1 = highway::Set(tag_double, xPtr1[i]);
                    y1 = highway::Set(tag_double, yPtr1[i]);
                    z1 = highway::Set(tag_double, zPtr1[i]);

                    if constexpr (vecPattern == VectorizationPattern::p1xVec) {

                        if constexpr(remainderJ) {
                            const auto restMaskDouble = restMasksDouble[restJ-1];
                            const auto restMaskLong = restMasksLong[restJ-1];

                            ownedJ = highway::MaskedLoad(restMaskLong, tag_long, &ownedPtr2[j]);
                            x2 = highway::MaskedLoad(restMaskDouble, tag_double, &xPtr2[j]);
                            y2 = highway::MaskedLoad(restMaskDouble, tag_double, &yPtr2[j]);
                            z2 = highway::MaskedLoad(restMaskDouble, tag_double, &zPtr2[j]);
                        }
                        else {
                            ownedJ = highway::LoadU(tag_long, &ownedPtr2[j]);
                            x2 = highway::LoadU(tag_double, &xPtr2[j]);
                            y2 = highway::LoadU(tag_double, &yPtr2[j]);
                            z2 = highway::LoadU(tag_double, &zPtr2[j]);
                        }
                    }
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
                        VectorLong tmpOwnedI;
                        VectorDouble tmpX1;
                        VectorDouble tmpY1;
                        VectorDouble tmpZ1;

                        if constexpr (remainderI) {
                            tmpOwnedI = _zeroLong;
                            tmpX1 = _zeroDouble;
                            tmpY1 = _zeroDouble;
                            tmpZ1 = _zeroDouble;
                        }
                        else {
                            // tmpOwnedI = highway::Set(tag_long, ownedPtr1[i+1]);
                            tmpX1 = highway::Set(tag_double, xPtr1[i+1]);
                            tmpY1 = highway::Set(tag_double, yPtr1[i+1]);
                            tmpZ1 = highway::Set(tag_double, zPtr1[i+1]);
                        }

                        ownedI = highway::ConcatLowerLower(tag_double, tmpOwnedI, ownedI);
                        x1 = highway::ConcatLowerLower(tag_double, tmpX1, x1);
                        y1 = highway::ConcatLowerLower(tag_double, tmpY1, y1);
                        z1 = highway::ConcatLowerLower(tag_double, tmpZ1, z1);

                        const auto maskLong = highway::FirstN(tag_long, remainderJ ? restJ : _vecLengthDouble / 2);
                        const auto maskDouble = highway::FirstN(tag_double, remainderJ ? restJ : _vecLengthDouble / 2);

                        // lower part of registers is filled
                        // TODO : figure out why LoadN does not work
                        // ownedJ = highway::MaskedLoad(maskLong, tag_long, &ownedPtr2[j]);
                        x2 = highway::MaskedLoad(maskDouble, tag_double, &xPtr2[j]);
                        y2 = highway::MaskedLoad(maskDouble, tag_double, &yPtr2[j]);
                        z2 = highway::MaskedLoad(maskDouble, tag_double, &zPtr2[j]);

                        // "broadcast" lower half to upper half
                        ownedJ = highway::ConcatLowerLower(tag_long, ownedJ, ownedJ);
                        x2 = highway::ConcatLowerLower(tag_double, x2, x2);
                        y2 = highway::ConcatLowerLower(tag_double, y2, y2);
                        z2 = highway::ConcatLowerLower(tag_double, z2, z2);
                    }
                    else if (vecPattern == VectorizationPattern::pVecDiv2x2) {
                        // TODO : implement
                    }
                    else if (vecPattern == VectorizationPattern::pVecx1) {
                        // TODO : implement
                    }
                    else if (vecPattern == VectorizationPattern::pVecxVec) {
                        // TODO : implement
                    }
                    else
                        throw std::runtime_error("No vectorization pattern matched, error!");

                    // leave this as it is -> noPPL handles this different anyways
                    if constexpr (useMixing) {
                        
                        double epsilon_buf[_vecLengthDouble] = {0.};
                        double sigma_buf[_vecLengthDouble] = {0.};
                        double shift_buf[_vecLengthDouble] = {0.};

                        for (size_t n = 0; (remainderJ ? n < restJ : n < _vecLengthDouble); ++n) {
                            epsilon_buf[n] = _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + n));
                            sigma_buf[n] =  _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + n));
                            if constexpr (applyShift) {
                                shift_buf[n] = _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + n));
                            }
                        }
                        epsilon24s = highway::LoadU(tag_double, epsilon_buf);
                        sigmaSquaredDiv2s = highway::LoadU(tag_double, sigma_buf);
                        if constexpr (applyShift) {
                            shift6s = highway::LoadU(tag_double, shift_buf);
                        }
                    }
                    else {
                        epsilon24s = _epsilon24;
                        sigmaSquaredDiv2s = _sigmaSquared;
                        if constexpr (applyShift) {
                            shift6s = _shift6;
                        }
                    }
                }

                template <bool remainder>
                inline void reduceAccumulatedForce(const VectorDouble& fxAcc, const VectorDouble& fyAcc, const VectorDouble& fzAcc,
                    double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr, const size_t i, const int rest = 0) {
                        
                    // TODO : think about more efficient way than for loop -> maybe highway provides something
                    
                    // TODO : handle case of rest
                    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
                        fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
                        fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
                        fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
                        double tmpX [_vecLengthDouble] = {0.};
                        double tmpY [_vecLengthDouble] = {0.};
                        double tmpZ [_vecLengthDouble] = {0.};
                        highway::StoreU(fxAcc, tag_double, tmpX);
                        highway::StoreU(fyAcc, tag_double, tmpY);
                        highway::StoreU(fzAcc, tag_double, tmpZ);

                        for (size_t n = 0; n < _vecLengthDouble; ++n) {
                            fxPtr[i + n%2] += tmpX[n];
                            fyPtr[i + n%2] += tmpY[n];
                            fzPtr[i + n%2] += tmpZ[n];
                        }
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
                        double tmpX [_vecLengthDouble] = {0.};
                        double tmpY [_vecLengthDouble] = {0.};
                        double tmpZ [_vecLengthDouble] = {0.};
                        highway::StoreU(fxAcc, tag_double, tmpX);
                        highway::StoreU(fyAcc, tag_double, tmpY);
                        highway::StoreU(fzAcc, tag_double, tmpZ);

                        for (size_t n = 0; n < _vecLengthDouble; ++n) {
                            fxPtr[i + n%(_vecLengthDouble/2)] += tmpX[n];
                            fyPtr[i + n%(_vecLengthDouble/2)] += tmpY[n];
                            fzPtr[i + n%(_vecLengthDouble/2)] += tmpZ[n];
                        }
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
                        double tmpX [_vecLengthDouble] = {0.};
                        double tmpY [_vecLengthDouble] = {0.};
                        double tmpZ [_vecLengthDouble] = {0.};
                        highway::StoreU(fxAcc, tag_double, tmpX);
                        highway::StoreU(fyAcc, tag_double, tmpY);
                        highway::StoreU(fzAcc, tag_double, tmpZ);

                        for (size_t n = 0; n < _vecLengthDouble; ++n) {
                            fxPtr[i + n] += tmpX[n];
                            fyPtr[i + n] += tmpY[n];
                            fzPtr[i + n] += tmpZ[n];
                        }
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
                        double tmpX [_vecLengthDouble] = {0.};
                        double tmpY [_vecLengthDouble] = {0.};
                        double tmpZ [_vecLengthDouble] = {0.};
                        highway::StoreU(fxAcc, tag_double, tmpX);
                        highway::StoreU(fyAcc, tag_double, tmpY);
                        highway::StoreU(fzAcc, tag_double, tmpZ);

                        for (size_t n = 0; n < _vecLengthDouble; ++n) {
                            fxPtr[i + n] += tmpX[n];
                            fyPtr[i + n] += tmpY[n];
                            fzPtr[i + n] += tmpZ[n];
                        }
                    }
                    else
                        throw std::runtime_error("No vectorization pattern matched, error!");
                }

                template <bool remainder>
                inline void handleNewton3Accumulation(const VectorDouble& fx, const VectorDouble& fy, const VectorDouble& fz,
                    double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr, const size_t j, const size_t rest = 0) {
                    
                    auto restMaskDouble = restMasksDouble[rest-1];
                    VectorDouble fx2 = _zeroDouble;
                    VectorDouble fy2 = _zeroDouble;
                    VectorDouble fz2 = _zeroDouble;

                    // TODO : handle different vectorization patterns

                    if constexpr (remainder) {
                        fx2 = highway::MaskedLoad(restMaskDouble, tag_double, &fxPtr[j]);
                        fy2 = highway::MaskedLoad(restMaskDouble, tag_double, &fyPtr[j]);
                        fz2 = highway::MaskedLoad(restMaskDouble, tag_double, &fzPtr[j]);
                    }
                    else {
                        fx2 = highway::LoadU(tag_double, &fxPtr[j]);
                        fy2 = highway::LoadU(tag_double, &fyPtr[j]);
                        fz2 = highway::LoadU(tag_double, &fzPtr[j]);
                    }

                    auto fx2New = fx2 - fx; // highway::Sub(fx2, fx);
                    auto fy2New = fy2 - fy; // highway::Sub(fy2, fy);
                    auto fz2New = fz2 - fz; // highway::Sub(fz2, fz);

                    if constexpr (remainder) {
                        highway::BlendedStore(fx2New, restMaskDouble, tag_double, &fxPtr[j]);                            
                        highway::BlendedStore(fy2New, restMaskDouble, tag_double, &fyPtr[j]);
                        highway::BlendedStore(fz2New, restMaskDouble, tag_double, &fzPtr[j]);
                    }
                    else {
                        highway::StoreU(fx2New, tag_double, &fxPtr[j]);
                        highway::StoreU(fy2New, tag_double, &fyPtr[j]);
                        highway::StoreU(fz2New, tag_double, &fzPtr[j]);
                    }
                }

                inline bool checkFirstLoopConditionSingle(const size_t i) {
                    // TODO : handle different vectorization patterns
                    return long(i) >=0;
                }

                inline bool checkFirstLoopConditionPair(const size_t i, const size_t size) {
                    // TODO : handle different vectorization patterns
                    return i < size;
                }

                inline bool checkSecondLoopCondition(const size_t i, const size_t j) {
                    // TODO : handle different vectorization patterns

                    // floor soa numParticles to multiple of vecLength
                    // If b is a power of 2 the following holds:
                    // a & ~(b -1) == a - (a mod b)
                    return j < (i & ~(_vecLengthDouble - 1));
                }

                inline int obtainFirstLoopRest(const size_t i) {
                    // TODO : handle different vectorization patterns
                    return 0;
                }

                inline int obtainSecondLoopRest(const size_t i, const size_t j) {
                    // TODO : handle different vectorization patterns

                    // If b is a power of 2 the following holds:
                    // a & (b -1) == a mod b
                    return (int)(i & (_vecLengthDouble - 1));
                }

                template <bool newton3>
                inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {

                    if (soa.size() == 0) return;

                    // obtain iterators for the various values
                    const auto *const __restrict xPtr = soa.template begin<Particle::AttributeNames::posX>();
                    const auto *const __restrict yPtr = soa.template begin<Particle::AttributeNames::posY>();
                    const auto *const __restrict zPtr = soa.template begin<Particle::AttributeNames::posZ>();

                    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

                    auto *const __restrict fxPtr = soa.template begin<Particle::AttributeNames::forceX>();
                    auto *const __restrict fyPtr = soa.template begin<Particle::AttributeNames::forceY>();
                    auto *const __restrict fzPtr = soa.template begin<Particle::AttributeNames::forceZ>();

                    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

                    // initialize and declare vector variables
                    auto virialSumX = highway::Zero(tag_double);
                    auto virialSumY = highway::Zero(tag_double);
                    auto virialSumZ = highway::Zero(tag_double);
                    auto uPotSum = highway::Zero(tag_double);

                    for (size_t i = soa.size() - 1; checkFirstLoopConditionSingle(i); decrementFirstLoop(i)) {

                        if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
                            continue;
                        }

                        static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                            "OwnershipStates underlying type should be int64_t!");

                        VectorDouble fxAcc = highway::Zero(tag_double);
                        VectorDouble fyAcc = highway::Zero(tag_double);
                        VectorDouble fzAcc = highway::Zero(tag_double);

                        VectorLong ownedStateI = highway::Zero(tag_long);
                        VectorDouble x1 = highway::Zero(tag_double);
                        VectorDouble y1 = highway::Zero(tag_double);
                        VectorDouble z1 = highway::Zero(tag_double);

                        VectorLong ownedStateJ = highway::Zero(tag_long);
                        VectorDouble x2 = highway::Zero(tag_double);
                        VectorDouble y2 = highway::Zero(tag_double);
                        VectorDouble z2 = highway::Zero(tag_double);

                        VectorDouble epsilon24s = highway::Zero(tag_double);
                        VectorDouble sigmaSquareds = highway::Zero(tag_double);
                        VectorDouble shift6s = highway::Zero(tag_double);

                        size_t j = 0;

                        for (; checkSecondLoopCondition(i, j); incrementSecondLoop(j)) {

                            fillVectorRegisters<false, false>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, shift6s, typeIDptr, typeIDptr, xPtr, yPtr, zPtr,
                                xPtr, yPtr, zPtr, reinterpret_cast<const int64_t *>(ownedStatePtr), reinterpret_cast<const int64_t *>(ownedStatePtr), i, j, 0, 0);
                            
                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;
                            
                            SoAKernel(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);
                            
                            handleNewton3Accumulation<false>(fx, fy, fz, fxPtr, fyPtr, fzPtr, j, 0);
                        }

                        const int restJ = obtainSecondLoopRest(i, j);
                        if (restJ > 0) {
                            fillVectorRegisters<false, true>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, shift6s, typeIDptr, typeIDptr, xPtr, yPtr, zPtr,
                                xPtr, yPtr, zPtr, reinterpret_cast<const int64_t *>(ownedStatePtr), reinterpret_cast<const int64_t *>(ownedStatePtr), i, j, 0, restJ);
                            
                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);

                            handleNewton3Accumulation<true>(fx, fy, fz, fxPtr, fyPtr, fzPtr, j, restJ);
                        }

                        reduceAccumulatedForce<false>(fxAcc, fyAcc, fzAcc, fxPtr, fyPtr, fzPtr, i);

                        if constexpr (calculateGlobals) {
                            // nothing yet
                        }
                    }
                }

                template <bool newton3>
                inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {

                    if (soa1.size() == 0 || soa2.size() == 0) {
                        return;
                    }

                    const auto *const __restrict x1Ptr = soa1.template begin<Particle::AttributeNames::posX>();
                    const auto *const __restrict y1Ptr = soa1.template begin<Particle::AttributeNames::posY>();
                    const auto *const __restrict z1Ptr = soa1.template begin<Particle::AttributeNames::posZ>();
                    const auto *const __restrict x2Ptr = soa2.template begin<Particle::AttributeNames::posX>();
                    const auto *const __restrict y2Ptr = soa2.template begin<Particle::AttributeNames::posY>();
                    const auto *const __restrict z2Ptr = soa2.template begin<Particle::AttributeNames::posZ>();

                    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
                    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

                    auto *const __restrict fx1Ptr = soa1.template begin<Particle::AttributeNames::forceX>();
                    auto *const __restrict fy1Ptr = soa1.template begin<Particle::AttributeNames::forceY>();
                    auto *const __restrict fz1Ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
                    auto *const __restrict fx2Ptr = soa2.template begin<Particle::AttributeNames::forceX>();
                    auto *const __restrict fy2Ptr = soa2.template begin<Particle::AttributeNames::forceY>();
                    auto *const __restrict fz2Ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

                    const auto *const __restrict typeID1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
                    const auto *const __restrict typeID2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

                    VectorDouble virialSumX = highway::Zero(tag_double);
                    VectorDouble virialSumY = highway::Zero(tag_double);
                    VectorDouble virialSumZ = highway::Zero(tag_double);
                    VectorDouble uPotSum = highway::Zero(tag_double);

                    for (size_t i = 0; checkFirstLoopConditionPair(i, soa1.size()); incrementFirstLoop(i)) {

                        if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
                            continue;
                        }

                        VectorDouble fxAcc = _zeroDouble;
                        VectorDouble fyAcc = _zeroDouble;
                        VectorDouble fzAcc = _zeroDouble;

                        VectorLong ownedStateI = _zeroLong;
                        VectorDouble x1 = _zeroDouble;
                        VectorDouble y1 = _zeroDouble;
                        VectorDouble z1 = _zeroDouble;

                        VectorLong ownedStateJ = _zeroLong;
                        VectorDouble x2 = _zeroDouble;
                        VectorDouble y2 = _zeroDouble;
                        VectorDouble z2 = _zeroDouble;

                        VectorDouble epsilon24s = _zeroDouble;
                        VectorDouble sigmaSquareds = _zeroDouble;
                        VectorDouble shift6s = _zeroDouble;

                        size_t j = 0;

                        for (; checkSecondLoopCondition(soa2.size(), j); incrementSecondLoop(j)) {
                            
                            fillVectorRegisters<false, false>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, shift6s, typeID1ptr, typeID2ptr, x1Ptr, y1Ptr, z1Ptr,
                                x2Ptr, y2Ptr, z2Ptr, reinterpret_cast<const int64_t *>(ownedStatePtr1), reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                i, j, 0, 0);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);

                            if constexpr (newton3) {
                                handleNewton3Accumulation<false>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j);
                            }
                        }

                        const int restJ = obtainSecondLoopRest(soa2.size(), j);
                        if (restJ > 0) {

                            fillVectorRegisters<false, true>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, shift6s, typeID1ptr, typeID2ptr, x1Ptr, y1Ptr, z1Ptr,
                                x2Ptr, y2Ptr, z2Ptr, reinterpret_cast<const int64_t *>(ownedStatePtr1), reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                i, j, 0, restJ);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);
                            
                            if constexpr (newton3) {
                                handleNewton3Accumulation<true>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
                            }
                        }

                        reduceAccumulatedForce<false>(fxAcc, fyAcc, fzAcc, fx1Ptr, fy1Ptr, fz1Ptr, i);
                    }

                    if constexpr (calculateGlobals) {        
                        // nothing yet
                    }                
                }

                inline void SoAKernel(const size_t j, const VectorLong& ownedStateI, const VectorLong& ownedStateJ,
                    const VectorDouble& x1, const VectorDouble& y1, const VectorDouble& z1, const VectorDouble& x2, const VectorDouble& y2,
                    const VectorDouble& z2, const VectorDouble& epsilon24s, const VectorDouble& sigmaSquareds,
                    VectorDouble& fxAcc, VectorDouble& fyAcc, VectorDouble& fzAcc, VectorDouble& fx, VectorDouble& fy, VectorDouble& fz,
                    VectorDouble& virialSumX, VectorDouble& virialSumY, VectorDouble& virialSumZ, VectorDouble& uPotSum) {
                    
                        // distance calculations
                    auto drX = x1 - x2; // highway::Sub(x1, x2);
                    auto drY = y1 - y2; // highway::Sub(y1, y2);
                    auto drZ = z1 - z2; // highway::Sub(z1, z2);
                    
                    auto drX2 = drX * drX; // highway::Mul(drX, drX);
                    auto drY2 = drY * drY; // highway::Mul(drY, drY);
                    auto drZ2 = drZ * drZ; // highway::Mul(drZ, drZ);

                    auto dr2 = drX2 + drY2 + drZ2;

                    auto cutoffMask = highway::Le(dr2, _cutoffSquared);

                    auto dummyMask = highway::Ne(ownedStateJ, _ownedStateDummy);
                    // convert long mask to double mask
                    VectorLong dummyMaskVec = highway::VecFromMask(dummyMask);
                    VectorDouble dummyMaskVecDouble = highway::ConvertTo(tag_double, dummyMaskVec);
                    auto dummyMaskDouble = highway::MaskFromVec(dummyMaskVecDouble);
                    auto cutoffDummyMask = highway::And(cutoffMask, dummyMaskDouble);

                    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
                        return;
                    }

                    // compute LJ Potential
                    auto invDr2 = _oneDouble / dr2; // highway::Div(_oneDouble, dr2);
                    auto lj2 = sigmaSquareds * invDr2; // highway::Mul(sigmaSquareds, invDr2);
                    auto lj4 = lj2 * lj2; // highway::Mul(lj2, lj2);
                    auto lj6 = lj2 * lj4; // highway::Mul(lj2, lj4);
                    auto lj12 = lj6 * lj6; // highway::Mul(lj6, lj6);
                    auto lj12m6 = lj12 - lj6; // highway::Sub(lj12, lj6);
                    auto lj12m6alj12 = lj12m6 - lj12; // highway::Add(lj12m6, lj12);
                    auto lj12m6alj12e = lj12m6alj12 * epsilon24s; // highway::Mul(lj12m6alj12, epsilon24s);
                    VectorDouble fac = lj12m6alj12e * invDr2; // highway::Mul(lj12m6alj12e, invDr2);

                    VectorDouble facMasked = highway::IfThenElse(cutoffDummyMask, fac, _zeroDouble);

                    fx = drX * facMasked; // highway::Mul(drX, facMasked);
                    fy = drY * facMasked; // highway::Mul(drY, facMasked);
                    fz = drZ * facMasked; // highway::Mul(drZ, facMasked);

                    fxAcc += fx; // highway::Add(fxAcc, fx);
                    fyAcc += fy; // highway::Add(fyAcc, fy);
                    fzAcc += fz; // highway::Add(fzAcc, fz);
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
                        
                        const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();
                        if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
                            return;
                        }

                        const auto *const __restrict xPtr = soa.template begin<Particle::AttributeNames::posX>();
                        const auto *const __restrict yPtr = soa.template begin<Particle::AttributeNames::posY>();
                        const auto *const __restrict zPtr = soa.template begin<Particle::AttributeNames::posZ>();

                        auto *const __restrict fxPtr = soa.template begin<Particle::AttributeNames::forceX>();
                        auto *const __restrict fyPtr = soa.template begin<Particle::AttributeNames::forceY>();
                        auto *const __restrict fzPtr = soa.template begin<Particle::AttributeNames::forceZ>();

                        const auto *const __restrict typeIDPtr = soa.template begin<Particle::AttributeNames::typeId>();

                        VectorDouble virialSumX = highway::Zero(tag_double);
                        VectorDouble virialSumY = highway::Zero(tag_double);
                        VectorDouble virialSumZ = highway::Zero(tag_double);
                        VectorDouble uPotSum = highway::Zero(tag_double);
                        VectorDouble fxAcc = highway::Zero(tag_double);
                        VectorDouble fyAcc = highway::Zero(tag_double);
                        VectorDouble fzAcc = highway::Zero(tag_double);

                        VectorDouble x1 = _zeroDouble;
                        VectorDouble y1 = _zeroDouble;
                        VectorDouble z1 = _zeroDouble;
                        VectorLong ownedStateI = _zeroLong;

                        VectorDouble x2 = _zeroDouble;
                        VectorDouble y2 = _zeroDouble;
                        VectorDouble z2 = _zeroDouble;
                        VectorLong ownedStateJ = _zeroLong;

                        VectorDouble epsilon24s = _zeroDouble;
                        VectorDouble sigmaSquareds = _zeroDouble;
                        VectorDouble shift6s = _zeroDouble;

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

                            // load neighbor particles in consecutive array
                            for (size_t vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
                                x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
                                y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
                                z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
                                typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
                                ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
                            }

                            fillVectorRegisters<false, false>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, shift6s, &typeIDPtr[indexFirst], typeID2Tmp.data(),
                                xPtr, yPtr, zPtr, x2Tmp.data(), y2Tmp.data(), z2Tmp.data(),
                                reinterpret_cast<const int64_t *>(ownedStatePtr), reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()),
                                indexFirst, j, 0, 0);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel(j, ownedStateI, ownedStateJ, x1, y1, z1, x2, y2, z2,
                                epsilon24s, sigmaSquareds, fxAcc, fyAcc, fzAcc, fx, fy, fz,
                                virialSumX, virialSumY, virialSumZ, uPotSum);

                            if constexpr (newton3) {
                                for (size_t vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
                                    fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
                                    fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
                                    fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
                                }

                                auto fx2 = highway::LoadU(tag_double, fx2Tmp.data());
                                auto fy2 = highway::LoadU(tag_double, fy2Tmp.data());
                                auto fz2 = highway::LoadU(tag_double, fz2Tmp.data());

                                auto fx2New = fx2 - fx;
                                auto fy2New = fy2 - fx;
                                auto fz2New = fz2 - fz;

                                highway::StoreU(fx2New, tag_double, fx2Tmp.data());
                                highway::StoreU(fy2New, tag_double, fy2Tmp.data());
                                highway::StoreU(fz2New, tag_double, fz2Tmp.data());

                                for (size_t vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
                                    fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
                                    fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
                                    fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
                                }
                            }
                        }

                        const size_t restJ = neighborList.size() & (_vecLengthDouble - 1);

                        if (restJ > 0) {
                            for (size_t vecIndex = 0; vecIndex < restJ; ++vecIndex) {
                                x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
                                y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
                                z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
                                typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
                                ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
                            }

                            fillVectorRegisters<false, true>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, shift6s, &typeIDPtr[indexFirst], typeID2Tmp.data(),
                                xPtr, yPtr, zPtr, x2Tmp.data(), y2Tmp.data(), z2Tmp.data(),
                                reinterpret_cast<const int64_t *>(ownedStatePtr), reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()),
                                 indexFirst, j, 0, 0);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel(j, ownedStateI, ownedStateJ, x1, y1, z1, x2, y2, z2,
                                epsilon24s, sigmaSquareds, fxAcc, fyAcc, fzAcc, fx, fy, fz,
                                virialSumX, virialSumY, virialSumZ, uPotSum);

                            if constexpr (newton3) {
                                for (size_t vecIndex = 0; vecIndex < restJ; ++vecIndex) {

                                    fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
                                    fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
                                    fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
                                }

                                auto fx2 = highway::LoadU(tag_double, fx2Tmp.data());
                                auto fy2 = highway::LoadU(tag_double, fy2Tmp.data());
                                auto fz2 = highway::LoadU(tag_double, fz2Tmp.data());

                                auto fx2New = fx2 - fx;
                                auto fy2New = fy2 - fx;
                                auto fz2New = fz2 - fz;

                                highway::StoreU(fx2New, tag_double, fx2Tmp.data());
                                highway::StoreU(fy2New, tag_double, fy2Tmp.data());
                                highway::StoreU(fz2New, tag_double, fz2Tmp.data());

                                for (size_t vecIndex = 0; vecIndex < restJ; ++vecIndex) {
                                    fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
                                    fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
                                    fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
                                }
                            }
                        }

                        reduceAccumulatedForce<false>(fxAcc, fyAcc, fzAcc, fxPtr, fyPtr, fzPtr, indexFirst);

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
                            _aosThreadData[threadnum].uPotSum += globals[3] * factor;
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
                    _epsilon24 = highway::Set(tag_double, epsilon24);
                    _sigmaSquared = highway::Set(tag_double, sigmaSquare);
                    if constexpr (applyShift) {
                        _shift6 = highway::Set(tag_double, ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, highway::GetLane(_cutoffSquared)));
                    } else {
                        _shift6 = _zeroDouble;
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
                const VectorLong _ownedStateDummy{highway::Zero(tag_long)};
                const VectorLong _ownedStateOwned{highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned))};
                const VectorDouble _cutoffSquared {};
                VectorDouble _shift6 {highway::Zero(tag_double)};
                VectorDouble _epsilon24 {highway::Zero(tag_double)};
                VectorDouble _sigmaSquared {highway::Zero(tag_double)};

                MaskDouble restMasksDouble[_vecLengthDouble-1];
                MaskLong restMasksLong[_vecLengthDouble-1];

                const double _cutoffSquareAoS {0.};
                double _epsilon24AoS, _sigmaSquareAoS, _shift6AoS = 0.;
                ParticlePropertiesLibrary<double, size_t>* _PPLibrary = nullptr;
                double _uPotSum {0.};
                std::array<double, 3> _virialSum;
                std::vector<AoSThreadData> _aosThreadData;
                bool _postProcessed;
    };
// } // Highway
} // mdLib
// HWY_AFTER_NAMESPACE();
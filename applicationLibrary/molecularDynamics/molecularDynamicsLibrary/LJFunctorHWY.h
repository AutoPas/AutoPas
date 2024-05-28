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
    const long _vecLengthDouble {highway::Lanes(tag_double)};
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
                _aosThreadData{},
                _postProcessed{false} {
                    if (calculateGlobals) {
                        _aosThreadData.resize(autopas::autopas_get_max_threads());
                    }

                    initializeRestMasks();

                    std::cout << hwy::TargetName(HWY_TARGET) << std::endl; // AutoPasLog(INFO, "Highway Wrapper initialized with a register size of ({}) with architecture {}.", _vecLengthDouble, hwy::TargetName(HWY_TARGET));
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

                    // TODO : figure out whether different patterns are really necessary -> maybe one can handle it different as well
                    for (size_t n = 0; n <_vecLengthDouble-1;++n) {
                        restMasksDouble[n] = highway::FirstN(tag_double, n+1);
                        restMasksLong[n] = highway::FirstN(tag_long, n+1);
                    }
                }
                
                inline void decrementFirstLoop(long& i) {

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

                inline void incrementFirstLoop(long& i) {

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

                inline void incrementSecondLoop(long& j) {

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

                template <bool remainder, bool reversed>
                inline void fillIRegisters(VectorDouble& x1, VectorDouble& y1, VectorDouble& z1,
                    VectorLong& ownedI, const int64_t *const __restrict ownedPtr1,
                    const double *const __restrict xPtr1, const double *const __restrict yPtr1,
                    const double *const __restrict zPtr1, const long i, const long rest) {

                        ownedI = highway::Set(tag_long, ownedPtr1[i]);
                        x1 = highway::Set(tag_double, xPtr1[i]);
                        y1 = highway::Set(tag_double, yPtr1[i]);
                        z1 = highway::Set(tag_double, zPtr1[i]);

                        if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
                            VectorLong tmpOwnedI = _zeroLong;
                            VectorDouble tmpX1 = _zeroDouble;
                            VectorDouble tmpY1 = _zeroDouble;
                            VectorDouble tmpZ1 = _zeroDouble;

                            if constexpr (!remainder) {
                                tmpOwnedI = highway::Set(tag_long, reversed ? ownedPtr1[i-1] : ownedPtr1[i+1]);
                                tmpX1 = highway::Set(tag_double, reversed ? xPtr1[i-1] : xPtr1[i+1]);
                                tmpY1 = highway::Set(tag_double, reversed ? yPtr1[i-1] : yPtr1[i+1]);
                                tmpZ1 = highway::Set(tag_double, reversed ? zPtr1[i-1] : zPtr1[i+1]);
                            }

                            ownedI = highway::ConcatLowerLower(tag_double, tmpOwnedI, ownedI);
                            x1 = highway::ConcatLowerLower(tag_double, tmpX1, x1);
                            y1 = highway::ConcatLowerLower(tag_double, tmpY1, y1);
                            z1 = highway::ConcatLowerLower(tag_double, tmpZ1, z1);
                        }
                        else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
                            const auto maskLong = restMasksLong[remainder ? rest-1 : _vecLengthDouble/2-1];
                            const auto maskDouble = restMasksDouble[remainder ? rest-1 : _vecLengthDouble/2-1];

                            long index = reversed ? (remainder ? 0 : i-_vecLengthDouble/2+1) : i;

                            ownedI = highway::MaskedLoad(maskLong, tag_long, &ownedPtr1[index]);
                            x1 = highway::MaskedLoad(maskDouble, tag_double, &xPtr1[index]);
                            y1 = highway::MaskedLoad(maskDouble, tag_double, &yPtr1[index]);
                            z1 = highway::MaskedLoad(maskDouble, tag_double, &zPtr1[index]);

                            ownedI = highway::ConcatLowerLower(tag_long, ownedI, ownedI);
                            x1 = highway::ConcatLowerLower(tag_double, x1, x1);
                            y1 = highway::ConcatLowerLower(tag_double, y1, y1);
                            z1 = highway::ConcatLowerLower(tag_double, z1, z1);
                        }
                        else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
                            long index = reversed ? (remainder ? 0 : i-_vecLengthDouble+1) : i;

                            if constexpr (remainder) {
                                auto maskDouble = restMasksDouble[rest-1];
                                auto maskLong = restMasksLong[rest-1];

                                ownedI = highway::MaskedLoad(maskLong, tag_long, &ownedPtr1[index]);
                                x1 = highway::MaskedLoad(maskDouble, tag_double, &xPtr1[index]);
                                y1 = highway::MaskedLoad(maskDouble, tag_double, &yPtr1[index]);
                                z1 = highway::MaskedLoad(maskDouble, tag_double, &zPtr1[index]);
                            }
                            else {
                                ownedI = highway::LoadU(tag_long, &ownedPtr1[index]);
                                x1 = highway::LoadU(tag_double, &xPtr1[index]);
                                y1 = highway::LoadU(tag_double, &yPtr1[index]);
                                z1 = highway::LoadU(tag_double, &zPtr1[index]);
                            }
                        }
                        else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
                            // TODO : implement
                            throw std::runtime_error("Not yet implemented!");
                        }

                }

                template <bool remainder>
                inline void fillJRegisters(VectorDouble& x2, VectorDouble& y2, VectorDouble& z2,
                    VectorLong& ownedJ, const int64_t *const __restrict ownedPtr2,
                    const double *const __restrict xPtr2, const double *const __restrict yPtr2,
                    const double *const __restrict zPtr2, const long j, const long rest) {

                    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
                        if constexpr(remainder) {
                            const auto restMaskDouble = restMasksDouble[rest-1];
                            const auto restMaskLong = restMasksLong[rest-1];

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

                        const auto maskLong = restMasksLong[remainder ? rest-1 : _vecLengthDouble/2-1];
                        const auto maskDouble = restMasksDouble[remainder ? rest-1 : _vecLengthDouble/2-1];

                        // lower part of registers is filled
                        ownedJ = highway::MaskedLoad(maskLong, tag_long, &ownedPtr2[j]);
                        x2 = highway::MaskedLoad(maskDouble, tag_double, &xPtr2[j]);
                        y2 = highway::MaskedLoad(maskDouble, tag_double, &yPtr2[j]);
                        z2 = highway::MaskedLoad(maskDouble, tag_double, &zPtr2[j]);

                        // "broadcast" lower half to upper half
                        ownedJ = highway::ConcatLowerLower(tag_long, ownedJ, ownedJ);
                        x2 = highway::ConcatLowerLower(tag_double, x2, x2);
                        y2 = highway::ConcatLowerLower(tag_double, y2, y2);
                        z2 = highway::ConcatLowerLower(tag_double, z2, z2);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {

                        ownedJ = highway::Set(tag_long, ownedPtr2[j]);
                        x2 = highway::Set(tag_double, xPtr2[j]);
                        y2 = highway::Set(tag_double, yPtr2[j]);
                        z2 = highway::Set(tag_double, zPtr2[j]);

                        VectorLong ownedJTmp =_zeroLong;
                        VectorDouble x2Tmp = _zeroDouble;
                        VectorDouble y2Tmp = _zeroDouble;
                        VectorDouble z2Tmp = _zeroDouble;

                        if constexpr (!remainder) {
                            ownedJTmp = highway::Set(tag_long, ownedPtr2[j+1]);
                            x2Tmp = highway::Set(tag_double, xPtr2[j+1]);
                            y2Tmp = highway::Set(tag_double, yPtr2[j+1]);
                            z2Tmp = highway::Set(tag_double, zPtr2[j+1]);
                        }

                        ownedJ = highway::ConcatLowerLower(tag_double, ownedJTmp, ownedJ);
                        x2 = highway::ConcatLowerLower(tag_double, x2Tmp, x2);
                        y2 = highway::ConcatLowerLower(tag_double, y2Tmp, y2);
                        z2 = highway::ConcatLowerLower(tag_double, z2Tmp, z2);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
                        ownedJ = highway::Set(tag_long, ownedPtr2[j]);
                        x2 = highway::Set(tag_double, xPtr2[j]);
                        y2 = highway::Set(tag_double, yPtr2[j]);
                        z2 = highway::Set(tag_double, zPtr2[j]);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
                        // TODO : implement
                        throw std::runtime_error("Not yet implemented!");
                    }

                }

                template <bool remainderI, bool remainderJ, bool reversed>
                inline void fillPhysicsRegisters(VectorDouble& epsilon24s, VectorDouble& sigmaSquaredDiv2s, VectorDouble& shift6s,
                        const long unsigned *const __restrict typeID1ptr, const long unsigned *const __restrict typeID2ptr,
                        const long restI, const long restJ) {

                    // TODO : handle different vectorization patterns
                    if constexpr (useMixing) {
                        
                        double epsilon_buf[_vecLengthDouble] = {0.};
                        double sigma_buf[_vecLengthDouble] = {0.};
                        double shift_buf[_vecLengthDouble] = {0.};

                        for (long n = 0; (remainderJ ? n < restJ : n < _vecLengthDouble); ++n) {
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

                template <bool remainderI, bool remainderJ, bool reversedI>
                inline void fillVectorRegisters(VectorDouble& x1, VectorDouble& y1, VectorDouble& z1,
                    VectorDouble& x2, VectorDouble& y2, VectorDouble& z2, VectorLong& ownedI, VectorLong& ownedJ,
                    VectorDouble& epsilon24s, VectorDouble& sigmaSquaredDiv2s, VectorDouble& shift6s,
                    const long unsigned *const __restrict typeID1ptr, const long unsigned *const __restrict typeID2ptr,
                    const double *const __restrict xPtr1, const double *const __restrict yPtr1, const double *const __restrict zPtr1,
                    const double *const __restrict xPtr2, const double *const __restrict yPtr2, const double *const __restrict zPtr2,
                    const int64_t *const __restrict ownedPtr1, const int64_t *const __restrict ownedPtr2, const long i, const long j,
                    const long restI = 0, const long restJ = 0) {

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
                        VectorLong tmpOwnedI = _zeroLong;
                        VectorDouble tmpX1 = _zeroDouble;
                        VectorDouble tmpY1 = _zeroDouble;
                        VectorDouble tmpZ1 = _zeroDouble;

                        if constexpr (!remainderI) {
                            tmpOwnedI = highway::Set(tag_long, reversedI ? ownedPtr1[i-1] : ownedPtr1[i+1]);
                            tmpX1 = highway::Set(tag_double, reversedI ? xPtr1[i-1] : xPtr1[i+1]);
                            tmpY1 = highway::Set(tag_double, reversedI ? yPtr1[i-1] : yPtr1[i+1]);
                            tmpZ1 = highway::Set(tag_double, reversedI ? zPtr1[i-1] : zPtr1[i+1]);
                        }

                        ownedI = highway::ConcatLowerLower(tag_double, tmpOwnedI, ownedI);
                        x1 = highway::ConcatLowerLower(tag_double, tmpX1, x1);
                        y1 = highway::ConcatLowerLower(tag_double, tmpY1, y1);
                        z1 = highway::ConcatLowerLower(tag_double, tmpZ1, z1);

                        const auto maskLong = restMasksLong[remainderJ ? restJ-1 : _vecLengthDouble/2-1];
                        const auto maskDouble = restMasksDouble[remainderJ ? restJ-1 : _vecLengthDouble/2-1];

                        // lower part of registers is filled
                        ownedJ = highway::MaskedLoad(maskLong, tag_long, &ownedPtr2[j]);
                        x2 = highway::MaskedLoad(maskDouble, tag_double, &xPtr2[j]);
                        y2 = highway::MaskedLoad(maskDouble, tag_double, &yPtr2[j]);
                        z2 = highway::MaskedLoad(maskDouble, tag_double, &zPtr2[j]);

                        // "broadcast" lower half to upper half
                        ownedJ = highway::ConcatLowerLower(tag_long, ownedJ, ownedJ);
                        x2 = highway::ConcatLowerLower(tag_double, x2, x2);
                        y2 = highway::ConcatLowerLower(tag_double, y2, y2);
                        z2 = highway::ConcatLowerLower(tag_double, z2, z2);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
                        
                        const auto maskLong = restMasksLong[remainderI ? restI-1 : _vecLengthDouble/2-1];
                        const auto maskDouble = restMasksDouble[remainderI ? restI-1 : _vecLengthDouble/2-1];

                        long index = reversedI ? (remainderI ? 0 : i-_vecLengthDouble/2+1) : i;

                        ownedI = highway::MaskedLoad(maskLong, tag_long, &ownedPtr1[index]);
                        x1 = highway::MaskedLoad(maskDouble, tag_double, &xPtr1[index]);
                        y1 = highway::MaskedLoad(maskDouble, tag_double, &yPtr1[index]);
                        z1 = highway::MaskedLoad(maskDouble, tag_double, &zPtr1[index]);

                        ownedI = highway::ConcatLowerLower(tag_long, ownedI, ownedI);
                        x1 = highway::ConcatLowerLower(tag_double, x1, x1);
                        y1 = highway::ConcatLowerLower(tag_double, y1, y1);
                        z1 = highway::ConcatLowerLower(tag_double, z1, z1);

                        ownedJ = highway::Set(tag_long, ownedPtr2[j]);
                        x2 = highway::Set(tag_double, xPtr2[j]);
                        y2 = highway::Set(tag_double, yPtr2[j]);
                        z2 = highway::Set(tag_double, zPtr2[j]);

                        VectorLong ownedJTmp =_zeroLong;
                        VectorDouble x2Tmp = _zeroDouble;
                        VectorDouble y2Tmp = _zeroDouble;
                        VectorDouble z2Tmp = _zeroDouble;

                        if constexpr (!remainderJ) {
                            ownedJTmp = highway::Set(tag_long, ownedPtr2[j+1]);
                            x2Tmp = highway::Set(tag_double, xPtr2[j+1]);
                            y2Tmp = highway::Set(tag_double, yPtr2[j+1]);
                            z2Tmp = highway::Set(tag_double, zPtr2[j+1]);
                        }

                        ownedJ = highway::ConcatLowerLower(tag_double, ownedJTmp, ownedJ);
                        x2 = highway::ConcatLowerLower(tag_double, x2Tmp, x2);
                        y2 = highway::ConcatLowerLower(tag_double, y2Tmp, y2);
                        z2 = highway::ConcatLowerLower(tag_double, z2Tmp, z2);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
                        ownedJ = highway::Set(tag_long, ownedPtr2[j]);
                        x2 = highway::Set(tag_double, xPtr2[j]);
                        y2 = highway::Set(tag_double, yPtr2[j]);
                        z2 = highway::Set(tag_double, zPtr2[j]);

                        long index = reversedI ? (remainderI ? 0 : i-_vecLengthDouble+1) : i;

                        if constexpr (remainderI) {
                            auto maskDouble = restMasksDouble[restI-1];
                            auto maskLong = restMasksLong[restI-1];

                            ownedI = highway::MaskedLoad(maskLong, tag_long, &ownedPtr1[index]);
                            x1 = highway::MaskedLoad(maskDouble, tag_double, &xPtr1[index]);
                            y1 = highway::MaskedLoad(maskDouble, tag_double, &yPtr1[index]);
                            z1 = highway::MaskedLoad(maskDouble, tag_double, &zPtr1[index]);
                        }
                        else {
                            ownedI = highway::LoadU(tag_long, &ownedPtr1[index]);
                            x1 = highway::LoadU(tag_double, &xPtr1[index]);
                            y1 = highway::LoadU(tag_double, &yPtr1[index]);
                            z1 = highway::LoadU(tag_double, &zPtr1[index]);
                        }
                    }
                    else if (vecPattern == VectorizationPattern::pVecxVec) {
                        // TODO : implement
                    }
                    else
                        throw std::runtime_error("No vectorization pattern matched, error!");

                    // TODO : handle different vectorization patterns
                    if constexpr (useMixing) {

                        // throw std::runtime_error("Should not be used!");
                        
                        double epsilon_buf[_vecLengthDouble] = {0.};
                        double sigma_buf[_vecLengthDouble] = {0.};
                        double shift_buf[_vecLengthDouble] = {0.};

                        for (long n = 0; (remainderJ ? n < restJ : n < _vecLengthDouble); ++n) {
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

                template <bool remainder, bool reversed>
                inline void reduceAccumulatedForce(const VectorDouble& fxAcc, const VectorDouble& fyAcc, const VectorDouble& fzAcc,
                    double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr, const long i, const int rest = 0) {
                        
                    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
                        fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
                        fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
                        fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
                        auto lowerFxAcc = highway::LowerHalf(fxAcc);
                        auto lowerFyAcc = highway::LowerHalf(fyAcc);
                        auto lowerFzAcc = highway::LowerHalf(fzAcc);

                        // TODO : handle this better, unfortunately cluster needs that...
                        auto extLowerFxAcc = highway::ZeroExtendVector(tag_double, lowerFxAcc);
                        auto extLowerFyAcc = highway::ZeroExtendVector(tag_double, lowerFyAcc);
                        auto extLowerFzAcc = highway::ZeroExtendVector(tag_double, lowerFzAcc);

                        fxPtr[i] += highway::ReduceSum(tag_double, extLowerFxAcc);
                        fyPtr[i] += highway::ReduceSum(tag_double, extLowerFyAcc);
                        fzPtr[i] += highway::ReduceSum(tag_double, extLowerFzAcc);

                        if constexpr (!remainder) {
                            auto higherFxAcc = highway::UpperHalf(tag_double, fxAcc);
                            auto higherFyAcc = highway::UpperHalf(tag_double, fyAcc);
                            auto higherFzAcc = highway::UpperHalf(tag_double, fzAcc);

                            // TODO : handle this better, unfortunately cluster needs that...
                            auto extHigherFxAcc = highway::ZeroExtendVector(tag_double, higherFxAcc);
                            auto extHigherFyAcc = highway::ZeroExtendVector(tag_double, higherFyAcc);
                            auto extHigherFzAcc = highway::ZeroExtendVector(tag_double, higherFzAcc);

                            long index = reversed ? i-1 : i+1;

                            fxPtr[index] += highway::ReduceSum(tag_double, extHigherFxAcc);
                            fyPtr[index] += highway::ReduceSum(tag_double, extHigherFyAcc);
                            fzPtr[index] += highway::ReduceSum(tag_double, extHigherFzAcc);
                        }
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
                        auto lowerFxAcc = highway::LowerHalf(fxAcc);
                        auto lowerFyAcc = highway::LowerHalf(fyAcc);
                        auto lowerFzAcc = highway::LowerHalf(fzAcc);

                        auto higherFxAcc = highway::UpperHalf(tag_double, fxAcc);
                        auto higherFyAcc = highway::UpperHalf(tag_double, fyAcc);
                        auto higherFzAcc = highway::UpperHalf(tag_double, fzAcc);

                        auto fxAccCombined = lowerFxAcc + higherFxAcc;
                        auto fyAccCombined = lowerFyAcc + higherFyAcc;
                        auto fzAccCombined = lowerFzAcc + higherFzAcc;

                        auto fxCombined = highway::ZeroExtendVector(tag_double, fxAccCombined);
                        auto fyCombined = highway::ZeroExtendVector(tag_double, fyAccCombined);
                        auto fzCombined = highway::ZeroExtendVector(tag_double, fzAccCombined);

                        // either load full slice of vecLength / 2 or a shortened slice of rest particles
                        // for reversed order, the index must be adapted
                        long index = reversed ? (remainder ? 0 : i-_vecLengthDouble/2+1) : i;
                        auto mask = restMasksDouble[remainder ? rest-1 : _vecLengthDouble/2-1];

                        auto oldFx = highway::MaskedLoad(mask, tag_double, &fxPtr[index]);
                        auto oldFy = highway::MaskedLoad(mask, tag_double, &fyPtr[index]);
                        auto oldFz = highway::MaskedLoad(mask, tag_double, &fzPtr[index]);

                        auto newFx = oldFx + fxCombined;
                        auto newFy = oldFy + fyCombined;
                        auto newFz = oldFz + fzCombined;

                        // TODO : get rid of blended store for non remainder case for the sake of performance
                        highway::BlendedStore(newFx, mask, tag_double, &fxPtr[index]);
                        highway::BlendedStore(newFy, mask, tag_double, &fyPtr[index]);
                        highway::BlendedStore(newFz, mask, tag_double, &fzPtr[index]);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {

                        if constexpr (reversed && remainder) {

                            fxPtr[rest-1] += highway::ExtractLane(fxAcc, rest-1);
                            fyPtr[rest-1] += highway::ExtractLane(fyAcc, rest-1);
                            fzPtr[rest-1] += highway::ExtractLane(fzAcc, rest-1);
                            return;
                        }
                        else if constexpr (reversed) {
                            fxPtr[i] += highway::ExtractLane(fxAcc, _vecLengthDouble-1);
                            fyPtr[i] += highway::ExtractLane(fyAcc, _vecLengthDouble-1);
                            fzPtr[i] += highway::ExtractLane(fzAcc, _vecLengthDouble-1);
                            return;
                        }

                        long index = reversed ? (remainder ? 0 : i-_vecLengthDouble+1) : i;
                        auto mask = restMasksDouble[remainder ? rest-1 : _vecLengthDouble-1];

                        VectorDouble oldFx = _zeroDouble;
                        VectorDouble oldFy = _zeroDouble;
                        VectorDouble oldFz = _zeroDouble;

                        oldFx = remainder ? highway::MaskedLoad(mask, tag_double, &fxPtr[index]) : highway::LoadU(tag_double, &fxPtr[index]);
                        oldFy = remainder ? highway::MaskedLoad(mask, tag_double, &fyPtr[index]) : highway::LoadU(tag_double, &fyPtr[index]);
                        oldFz = remainder ? highway::MaskedLoad(mask, tag_double, &fzPtr[index]) : highway::LoadU(tag_double, &fzPtr[index]);

                        auto newFx = oldFx + fxAcc;
                        auto newFy = oldFy + fyAcc;
                        auto newFz = oldFz + fzAcc;

                        if constexpr (remainder) {
                            highway::BlendedStore(newFx, mask, tag_double, &fxPtr[index]);
                            highway::BlendedStore(newFy, mask, tag_double, &fyPtr[index]);
                            highway::BlendedStore(newFz, mask, tag_double, &fzPtr[index]);
                        }
                        else {
                            highway::StoreU(newFx, tag_double, &fxPtr[index]);
                            highway::StoreU(newFy, tag_double, &fyPtr[index]);
                            highway::StoreU(newFz, tag_double, &fzPtr[index]);
                        }
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
                        // TODO : compare with previous approach regarding performance
                        // this approach is much more convenient to understand & maintain
                        double tmpX [_vecLengthDouble] = {0.};
                        double tmpY [_vecLengthDouble] = {0.};
                        double tmpZ [_vecLengthDouble] = {0.};
                        highway::StoreU(fxAcc, tag_double, tmpX);
                        highway::StoreU(fyAcc, tag_double, tmpY);
                        highway::StoreU(fzAcc, tag_double, tmpZ);

                        for (long n = 0; n < _vecLengthDouble; ++n) {
                            long index = reversed ? i-n : i+n;
                            fxPtr[index] += tmpX[n];
                            fyPtr[index] += tmpY[n];
                            fzPtr[index] += tmpZ[n];
                        }
                    }
                    else
                        throw std::runtime_error("No vectorization pattern matched, error!");
                }

                template <bool remainder, bool reversed>
                inline void handleNewton3Accumulation(const VectorDouble& fx, const VectorDouble& fy, const VectorDouble& fz,
                    double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr, const long j, const long rest = 0) {
                    
                    VectorDouble fx2 = _zeroDouble;
                    VectorDouble fy2 = _zeroDouble;
                    VectorDouble fz2 = _zeroDouble;
                    if constexpr (vecPattern == VectorizationPattern::p1xVec) {

                        auto mask = restMasksDouble[rest-1];

                        if constexpr (remainder) {
                            fx2 = highway::MaskedLoad(mask, tag_double, &fxPtr[j]);
                            fy2 = highway::MaskedLoad(mask, tag_double, &fyPtr[j]);
                            fz2 = highway::MaskedLoad(mask, tag_double, &fzPtr[j]);
                        }
                        else {
                            fx2 = highway::LoadU(tag_double, &fxPtr[j]);
                            fy2 = highway::LoadU(tag_double, &fyPtr[j]);
                            fz2 = highway::LoadU(tag_double, &fzPtr[j]);
                        }

                        auto fx2New = fx2 - fx;
                        auto fy2New = fy2 - fy;
                        auto fz2New = fz2 - fz;

                        if constexpr (remainder) {
                            highway::BlendedStore(fx2New, mask, tag_double, &fxPtr[j]);                            
                            highway::BlendedStore(fy2New, mask, tag_double, &fyPtr[j]);
                            highway::BlendedStore(fz2New, mask, tag_double, &fzPtr[j]);
                        }
                        else {
                            highway::StoreU(fx2New, tag_double, &fxPtr[j]);
                            highway::StoreU(fy2New, tag_double, &fyPtr[j]);
                            highway::StoreU(fz2New, tag_double, &fzPtr[j]);
                        }
                    }
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
                        auto lowerFx = highway::LowerHalf(fx);
                        auto lowerFy = highway::LowerHalf(fy);
                        auto lowerFz = highway::LowerHalf(fz);

                        auto higherFx = highway::UpperHalf(tag_double, fx);
                        auto higherFy = highway::UpperHalf(tag_double, fy);
                        auto higherFz = highway::UpperHalf(tag_double, fz);

                        auto fxCombined = lowerFx + higherFx;
                        auto fyCombined = lowerFy + higherFy;
                        auto fzCombined = lowerFz + higherFz;

                        VectorDouble fxCombined_ext = highway::ZeroExtendVector(tag_double, fxCombined);
                        VectorDouble fyCombined_ext = highway::ZeroExtendVector(tag_double, fyCombined);
                        VectorDouble fzCombined_ext = highway::ZeroExtendVector(tag_double, fzCombined);

                        auto mask = restMasksDouble[remainder ? rest-1 : _vecLengthDouble/2-1];
                        fx2 = highway::MaskedLoad(mask, tag_double, &fxPtr[j]);
                        fy2 = highway::MaskedLoad(mask, tag_double, &fyPtr[j]);
                        fz2 = highway::MaskedLoad(mask, tag_double, &fzPtr[j]);

                        auto newFx = fx2 - fxCombined_ext;
                        auto newFy = fy2 - fyCombined_ext;
                        auto newFz = fz2 - fzCombined_ext;

                        highway::BlendedStore(newFx, mask, tag_double, &fxPtr[j]);
                        highway::BlendedStore(newFy, mask, tag_double, &fyPtr[j]);
                        highway::BlendedStore(newFz, mask, tag_double, &fzPtr[j]);                        
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
                        auto lowerFx = highway::LowerHalf(fx);
                        auto lowerFy = highway::LowerHalf(fy);
                        auto lowerFz = highway::LowerHalf(fz);

                        // TODO : handle this better, unfortunately cluster needs that...
                        auto extLowerFx = highway::ZeroExtendVector(tag_double, lowerFx);
                        auto extLowerFy = highway::ZeroExtendVector(tag_double, lowerFy);
                        auto extLowerFz = highway::ZeroExtendVector(tag_double, lowerFz);

                        fxPtr[j] -= highway::ReduceSum(tag_double, extLowerFx);
                        fyPtr[j] -= highway::ReduceSum(tag_double, extLowerFy);
                        fzPtr[j] -= highway::ReduceSum(tag_double, extLowerFz);

                        if constexpr (!remainder) {
                            auto higherFx = highway::UpperHalf(tag_double, fx);
                            auto higherFy = highway::UpperHalf(tag_double, fy);
                            auto higherFz = highway::UpperHalf(tag_double, fz);

                            // TODO : handle this better, unfortunately cluster needs that...
                            auto extHigherFx = highway::ZeroExtendVector(tag_double, higherFx);
                            auto extHigherFy = highway::ZeroExtendVector(tag_double, higherFy);
                            auto extHigherFz = highway::ZeroExtendVector(tag_double, higherFz);

                            fxPtr[j+1] -= highway::ReduceSum(tag_double, extHigherFx);
                            fyPtr[j+1] -= highway::ReduceSum(tag_double, extHigherFy);
                            fzPtr[j+1] -= highway::ReduceSum(tag_double, extHigherFz);
                        }
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
                        fxPtr[j] -= highway::ReduceSum(tag_double, fx);
                        fyPtr[j] -= highway::ReduceSum(tag_double, fy);
                        fzPtr[j] -= highway::ReduceSum(tag_double, fz);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
                        // TODO : implement
                    }
                }

                template <bool reversed>
                inline bool checkFirstLoopCondition(const long i, const long size = 0) {

                    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
                        return reversed ? i >=1 : i < size;
                    }
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
                        return reversed ? i >= 2 : i < size-1;
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
                        return reversed ? i >= _vecLengthDouble/2 : i < size-(_vecLengthDouble/2-1);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
                        return reversed ? i >= _vecLengthDouble : i < size-(_vecLengthDouble-1);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
                        return reversed ? i >= _vecLengthDouble : i < size-_vecLengthDouble+1;
                    }
                }

                template <bool reversed>
                inline bool checkSecondLoopCondition(const long i, const long j) {

                    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
                        return j < (i-_vecLengthDouble+1);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
                        return j < (i-_vecLengthDouble/2+1);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
                        return j < (i-1);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
                        return j < i;
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
                        return j < (i & ~(_vecLengthDouble - 1));
                    }
                }

                template <bool reversed>
                inline int obtainFirstLoopRest(const long i, const long size = 0) {
                    return reversed ? (i < 0 ? 0 : i+1) : size-i;
                }

                template <bool reversed>
                inline int obtainSecondLoopRest(const long j) {
                    
                    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
                        return (int)(j & (_vecLengthDouble-1));
                    }
                    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
                        return (int)(j & (_vecLengthDouble/2-1));
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
                        return (int)(j & 1);
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
                        return 0;
                    }
                    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
                        return (int)(j & (_vecLengthDouble - 1));
                    }
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
                    auto virialSumX = _zeroDouble;
                    auto virialSumY = _zeroDouble;
                    auto virialSumZ = _zeroDouble;
                    auto uPotSum = _zeroDouble;

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

                    long i = soa.size() - 1;
                    for (; checkFirstLoopCondition<true>(i); decrementFirstLoop(i)) {

                        /*
                        if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {

                            // TODO : handle different vec patterns
                            // continue;
                        }
                        */

                        static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                            "OwnershipStates underlying type should be int64_t!");

                        fxAcc = _zeroDouble;
                        fyAcc = _zeroDouble;
                        fzAcc = _zeroDouble;

                        fillIRegisters<false, true>(x1, y1, z1, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr),
                            xPtr, yPtr, zPtr, i, 0);

                        long j = 0;

                        for (; checkSecondLoopCondition<true>(i, j); incrementSecondLoop(j)) {

                            fillJRegisters<false>(x2, y2, z2, ownedStateJ, reinterpret_cast<const int64_t *>(ownedStatePtr),
                                xPtr, yPtr, zPtr, j, 0);

                            fillPhysicsRegisters<false, false, true>(epsilon24s, sigmaSquareds, shift6s, typeIDptr, typeIDptr, 0, 0);
                            
                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;
                            
                            SoAKernel<true>(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);
                            
                            handleNewton3Accumulation<false, true>(fx, fy, fz, fxPtr, fyPtr, fzPtr, j);
                        }

                        const int restJ = obtainSecondLoopRest<true>(i);
                        if (restJ > 0) {
                        
                            fillJRegisters<true>(x2, y2, z2, ownedStateJ, reinterpret_cast<const int64_t *>(ownedStatePtr),
                                xPtr, yPtr, zPtr, j, restJ);

                            fillPhysicsRegisters<false, true, true>(epsilon24s, sigmaSquareds, shift6s, typeIDptr, typeIDptr, 0, restJ);
                            
                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<true>(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);

                            handleNewton3Accumulation<true, true>(fx, fy, fz, fxPtr, fyPtr, fzPtr, j, restJ);
                        }

                        reduceAccumulatedForce<false, true>(fxAcc, fyAcc, fzAcc, fxPtr, fyPtr, fzPtr, i);
                    }

                    /*
                    for (; i >= 0 && ownedStatePtr[i] == autopas::OwnershipState::dummy; --i) {
                        // Skip all dummy particles at the end of the list
                    }
                    */
                    const long restI = obtainFirstLoopRest<true>(i);
                    if (restI > 0) {
                        fxAcc = _zeroDouble;
                        fyAcc = _zeroDouble;
                        fzAcc = _zeroDouble;

                        fillIRegisters<true, true>(x1, y1, z1, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr),
                            xPtr, yPtr, zPtr, i, restI);

                        long j = 0;

                        for (; checkSecondLoopCondition<true>(i, j); incrementSecondLoop(j)) {

                            fillJRegisters<false>(x2, y2, z2, ownedStateJ, reinterpret_cast<const int64_t *>(ownedStatePtr),
                                xPtr, yPtr, zPtr, j, 0);

                            fillPhysicsRegisters<false, true, true>(epsilon24s, sigmaSquareds, shift6s, typeIDptr, typeIDptr, restI, 0);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<true>(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);

                            handleNewton3Accumulation<false, true>(fx, fy, fz, fxPtr, fyPtr, fzPtr, j);
                        }

                        const int restJ = obtainSecondLoopRest<true>(i);
                        if (restJ > 0) {

                            fillJRegisters<true>(x2, y2, z2, ownedStateJ, reinterpret_cast<const int64_t *>(ownedStatePtr),
                                xPtr, yPtr, zPtr, j, restJ);

                            fillPhysicsRegisters<true, true, true>(epsilon24s, sigmaSquareds, shift6s, typeIDptr, typeIDptr, restI, restJ);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<true>(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);
                            
                            handleNewton3Accumulation<true, true>(fx, fy, fz, fxPtr, fyPtr, fzPtr, j, restJ);
                        }

                        reduceAccumulatedForce<true, true>(fxAcc, fyAcc, fzAcc, fxPtr, fyPtr, fzPtr, i, restI);
                    }

                    if constexpr (calculateGlobals) {
                        const int threadnum = autopas::autopas_get_thread_num();

                        double globals[4] {
                            highway::ReduceSum(tag_double, virialSumX),
                            highway::ReduceSum(tag_double, virialSumY),
                            highway::ReduceSum(tag_double, virialSumZ),
                            highway::ReduceSum(tag_double, uPotSum)
                        };

                        double factor = 1.;
                        if constexpr (newton3) {
                            factor = 0.5;
                        }

                        _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
                        _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
                        _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
                        _aosThreadData[threadnum].uPotSum += globals[3] * factor;
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

                    VectorDouble virialSumX = _zeroDouble;
                    VectorDouble virialSumY = _zeroDouble;
                    VectorDouble virialSumZ = _zeroDouble;
                    VectorDouble uPotSum = _zeroDouble;

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

                    long i = 0;
                    for (; checkFirstLoopCondition<false>(i, soa1.size()); incrementFirstLoop(i)) {

                        /*
                        if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {

                            // TODO : handle different vec patterns
                            // continue;
                        }
                        */

                        fxAcc = _zeroDouble;
                        fyAcc = _zeroDouble;
                        fzAcc = _zeroDouble;

                        fillIRegisters<false, false>(x1, y1, z1, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr1),
                            x1Ptr, y1Ptr, z1Ptr, i, 0);

                        long j = 0;

                        for (; checkSecondLoopCondition<false>(soa2.size(), j); incrementSecondLoop(j)) {

                            fillJRegisters<false>(x2, y2, z2, ownedStateJ, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                x2Ptr, y2Ptr, z2Ptr, j, 0);

                            fillPhysicsRegisters<false, false, false>(epsilon24s, sigmaSquareds, shift6s, typeID1ptr, typeID2ptr, 0, 0);
                            
                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<newton3>(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);

                            if constexpr (newton3) {
                                handleNewton3Accumulation<false, false>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, 0);
                            }
                        }

                        const int restJ = obtainSecondLoopRest<false>(soa2.size());
                        if (restJ > 0) {

                            fillJRegisters<true>(x2, y2, z2, ownedStateJ, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                x2Ptr, y2Ptr, z2Ptr, j, restJ);

                            fillPhysicsRegisters<false, true, false>(epsilon24s, sigmaSquareds, shift6s, typeID1ptr, typeID2ptr, 0, restJ);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<newton3>(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);
                            
                            if constexpr (newton3) {
                                handleNewton3Accumulation<true, false>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
                            }
                        }

                        reduceAccumulatedForce<false, false>(fxAcc, fyAcc, fzAcc, fx1Ptr, fy1Ptr, fz1Ptr, i);
                    }
                    /*
                    for (; ownedStatePtr1[i] == autopas::OwnershipState::dummy && i < soa1.size()-1;) {
                        ++i;
                    }
                    */
                    const long restI = obtainFirstLoopRest<false>(i, soa1.size());
                    if (restI > 0) {
                        fxAcc = _zeroDouble;
                        fyAcc = _zeroDouble;
                        fzAcc = _zeroDouble;

                        fillIRegisters<true, false>(x1, y1, z1, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr1),
                            x1Ptr, y1Ptr, z1Ptr, i, restI);

                        long j = 0;

                        for (; checkSecondLoopCondition<false>(soa2.size(), j); incrementSecondLoop(j)) {

                            fillJRegisters<false>(x2, y2, z2, ownedStateJ, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                x2Ptr, y2Ptr, z2Ptr, j, 0);

                            fillPhysicsRegisters<false, false, false>(epsilon24s, sigmaSquareds, shift6s, typeID1ptr, typeID2ptr, restI, 0);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<newton3>(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);

                            if constexpr (newton3) {
                                handleNewton3Accumulation<false, false>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, 0);
                            }
                        }

                        const int restJ = obtainSecondLoopRest<false>(soa2.size());
                        if (restJ > 0) {

                            fillJRegisters<true>(x2, y2, z2, ownedStateJ, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                x2Ptr, y2Ptr, z2Ptr, j, restJ);

                            fillPhysicsRegisters<true, true, false>(epsilon24s, sigmaSquareds, shift6s, typeID1ptr, typeID2ptr, restI, restJ);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<newton3>(j, ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc,
                                fx, fy, fz, virialSumX, virialSumY, virialSumZ, uPotSum);
                            
                            if constexpr (newton3) {
                                handleNewton3Accumulation<true, false>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
                            }
                        }

                        reduceAccumulatedForce<true, false>(fxAcc, fyAcc, fzAcc, fx1Ptr, fy1Ptr, fz1Ptr, i, restI);
                    }

                    if constexpr (calculateGlobals) {
                        const int threadnum = autopas::autopas_get_thread_num();

                        double globals[4] {
                            highway::ReduceSum(tag_double, virialSumX),
                            highway::ReduceSum(tag_double, virialSumY),
                            highway::ReduceSum(tag_double, virialSumZ),
                            highway::ReduceSum(tag_double, uPotSum)
                        };

                        double factor = 1.;
                        if constexpr (newton3) {
                            factor = 0.5;
                        }

                        _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
                        _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
                        _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
                        _aosThreadData[threadnum].uPotSum += globals[3] * factor;
                    }
                }

                template <bool newton3>
                inline void SoAKernel(const long j, const VectorLong& ownedStateI, const VectorLong& ownedStateJ,
                    const VectorDouble& x1, const VectorDouble& y1, const VectorDouble& z1, const VectorDouble& x2, const VectorDouble& y2,
                    const VectorDouble& z2, const VectorDouble& epsilon24s, const VectorDouble& sigmaSquareds, const VectorDouble& shift6s,
                    VectorDouble& fxAcc, VectorDouble& fyAcc, VectorDouble& fzAcc, VectorDouble& fx, VectorDouble& fy, VectorDouble& fz,
                    VectorDouble& virialSumX, VectorDouble& virialSumY, VectorDouble& virialSumZ, VectorDouble& uPotSum) {
                    
                    // distance calculations
                    auto drX = x1 - x2;
                    auto drY = y1 - y2;
                    auto drZ = z1 - z2;
                    
                    auto drX2 = drX * drX;
                    auto drY2 = drY * drY;
                    auto drZ2 = drZ * drZ;

                    auto dr2 = drX2 + drY2 + drZ2;

                    auto cutoffMask = highway::Le(dr2, _cutoffSquared);
                    auto zeroMask = highway::Ne(dr2, _zeroDouble); // TODO : discuss this approach

                    auto combinedCutoffMask = highway::And(cutoffMask, zeroMask);

                    auto dummyMaskI = highway::Ne(ownedStateI, _ownedStateDummy);
                    auto dummyMaskJ = highway::Ne(ownedStateJ, _ownedStateDummy);
                    auto dummyMask = highway::And(dummyMaskI, dummyMaskJ);
                    // convert long mask to double mask
                    VectorLong dummyMaskVec = highway::VecFromMask(dummyMask);
                    VectorDouble dummyMaskVecDouble = highway::ConvertTo(tag_double, dummyMaskVec);
                    auto dummyMaskDouble = highway::MaskFromVec(dummyMaskVecDouble);
                    auto cutoffDummyMask = highway::And(combinedCutoffMask, dummyMaskDouble);

                    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
                        return;
                    }

                    // compute LJ Potential
                    auto invDr2 = _oneDouble / dr2;
                    auto lj2 = sigmaSquareds * invDr2;
                    auto lj4 = lj2 * lj2;
                    auto lj6 = lj2 * lj4;
                    auto lj12 = lj6 * lj6;
                    auto lj12m6 = lj12 - lj6;
                    auto lj12m6alj12 = lj12m6 + lj12;
                    auto lj12m6alj12e = lj12m6alj12 * epsilon24s;
                    VectorDouble fac = lj12m6alj12e * invDr2;

                    VectorDouble facMasked = highway::IfThenElseZero(cutoffDummyMask, fac);

                    fx = drX * facMasked;
                    fy = drY * facMasked;
                    fz = drZ * facMasked;

                    fxAcc = fxAcc + fx;
                    fyAcc = fyAcc + fy;
                    fzAcc = fzAcc + fz;

                    if constexpr (calculateGlobals) {

                        auto virialX = fx * drX;
                        auto virialY = fy * drY;
                        auto virialZ = fz * drZ;

                        auto uPot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
                        auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);
                        
                        auto ownedMaskI = highway::Eq(highway::ConvertTo(tag_double, ownedStateI), highway::ConvertTo(tag_double, _ownedStateOwned));
                        auto energyFactor = highway::IfThenElse(ownedMaskI, _oneDouble, _zeroDouble);

                        if constexpr (newton3) {
                            auto ownedMaskJ = highway::Eq(highway::ConvertTo(tag_double, ownedStateJ), highway::ConvertTo(tag_double, _ownedStateOwned));
                            energyFactor = energyFactor + highway::IfThenElse(ownedMaskJ, _oneDouble, _zeroDouble);
                        }

                        uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
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

                        long j = 0;

                        for (; j < (neighborList.size() & ~(_vecLengthDouble - 1)); j += _vecLengthDouble) {

                            // load neighbor particles in consecutive array
                            for (long vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
                                x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
                                y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
                                z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
                                if constexpr (newton3) {
                                    fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
                                    fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
                                    fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
                                }
                                typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
                                ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
                            }

                            fillVectorRegisters<false, false, false>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, shift6s, &typeIDPtr[indexFirst], typeID2Tmp.data(),
                                xPtr, yPtr, zPtr, x2Tmp.data(), y2Tmp.data(), z2Tmp.data(),
                                reinterpret_cast<const int64_t *>(ownedStatePtr), reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()),
                                indexFirst, j, 0, 0);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<newton3>(j, ownedStateI, ownedStateJ, x1, y1, z1, x2, y2, z2,
                                epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc, fx, fy, fz,
                                virialSumX, virialSumY, virialSumZ, uPotSum);

                            if constexpr (newton3) {

                                handleNewton3Accumulation<false, false>(fx, fy, fz, fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(), 0, 0);

                                for (long vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
                                    fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
                                    fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
                                    fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
                                }
                            }
                        }

                        const long restJ = static_cast<int>(neighborList.size() & (_vecLengthDouble - 1));

                        if (restJ > 0) {
                            for (long vecIndex = 0; vecIndex < restJ; ++vecIndex) {
                                x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
                                y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
                                z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
                                if constexpr (newton3) {
                                    fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
                                    fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
                                    fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
                                }
                                typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
                                ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
                            }

                            fillVectorRegisters<false, true, false>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, shift6s, &typeIDPtr[indexFirst], typeID2Tmp.data(),
                                xPtr, yPtr, zPtr, x2Tmp.data(), y2Tmp.data(), z2Tmp.data(),
                                reinterpret_cast<const int64_t *>(ownedStatePtr), reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()),
                                 indexFirst, j, 0, restJ);

                            VectorDouble fx = _zeroDouble;
                            VectorDouble fy = _zeroDouble;
                            VectorDouble fz = _zeroDouble;

                            SoAKernel<newton3>(j, ownedStateI, ownedStateJ, x1, y1, z1, x2, y2, z2,
                                epsilon24s, sigmaSquareds, shift6s, fxAcc, fyAcc, fzAcc, fx, fy, fz,
                                virialSumX, virialSumY, virialSumZ, uPotSum);

                            if constexpr (newton3) {

                                handleNewton3Accumulation<false, false>(fx, fy, fz, fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(), 0, restJ);

                                for (long vecIndex = 0; vecIndex < _vecLengthDouble && vecIndex < restJ; ++vecIndex) {
                                    fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
                                    fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
                                    fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
                                }
                            }
                        }

                        reduceAccumulatedForce<false, false>(fxAcc, fyAcc, fzAcc, fxPtr, fyPtr, fzPtr, indexFirst);

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
// Created by Luis Gall on 27.03.2024

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "VectorizationPatterns.h"

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
    typedef decltype(highway::Zero(tag_double)) VectorDouble;
    typedef decltype(highway::Zero(tag_long)) VectorLong;

    template <class Particle, bool applyShift = false, bool useMixing = false,
        autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
        bool relevantForTuning = true>

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
                inline void decrementFirstLoop(size_t& i) {

                    switch (_vectorizationPattern)
                    {
                    case VectorizationPattern::p1xVec:
                        --i;
                        break;
                    
                    case VectorizationPattern::p2xVecDiv2:
                        i = i-2;
                        break;
                    case VectorizationPattern::pVecDiv2x2:
                        i = i - (_vecLengthDouble / 2);
                        break;
                    case VectorizationPattern::pVecx1:
                        i = i - _vecLengthDouble;
                        break;
                    case VectorizationPattern::pVecxVec:
                        i = i - _vecLengthDouble;
                        break;
                    default:
                        throw std::runtime_error("No vectorization pattern matched, error!");
                    }
                }

                inline void incrementFirstLoop(size_t& i) {
                    switch (_vectorizationPattern)
                    {
                    case VectorizationPattern::p1xVec:
                        ++i;
                        break;
                    
                    case VectorizationPattern::p2xVecDiv2:
                        i = i + 2;
                        break;
                    case VectorizationPattern::pVecDiv2x2:
                        i = i + (_vecLengthDouble / 2);
                        break;
                    case VectorizationPattern::pVecx1:
                        i = i + _vecLengthDouble;
                        break;
                    case VectorizationPattern::pVecxVec:
                        i = i + _vecLengthDouble;
                        break;
                    default:
                        throw std::runtime_error("No vectorization pattern matched, error!");
                    }
                }

                inline void incrementSecondLoop(size_t& j) {
                    switch (_vectorizationPattern)
                    {
                    case VectorizationPattern::p1xVec:
                        j = j + _vecLengthDouble;
                        break;
                    
                    case VectorizationPattern::p2xVecDiv2:
                        j = j + _vecLengthDouble / 2;
                        break;
                    case VectorizationPattern::pVecDiv2x2:
                        j = j + 2;
                        break;
                    case VectorizationPattern::pVecx1:
                        ++j;
                        break;
                    case VectorizationPattern::pVecxVec:
                        j = j + _vecLengthDouble;
                        break;
                    default:
                        throw std::runtime_error("No vectorization pattern matched, error!");
                    }
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

                    switch (_vectorizationPattern)
                    {
                    case VectorizationPattern::p1xVec: {

                        if constexpr (remainderJ) {
                            const auto restMaskDouble = highway::FirstN(tag_double, restJ);
                            const auto restMaskLong = highway::FirstN(tag_long, restJ);

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
                        break;
                    case VectorizationPattern::p2xVecDiv2: {
                        
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
                            tmpOwnedI = highway::Set(tag_long, ownedPtr1[i+1]);
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
                        ownedJ = highway::MaskedLoad(maskLong, tag_long, &ownedPtr2[j]);
                        x2 = highway::MaskedLoad(maskDouble, tag_double, &xPtr2[j]);
                        y2 = highway::MaskedLoad(maskDouble, tag_double, &yPtr2[j]);
                        z2 = highway::MaskedLoad(maskDouble, tag_double, &zPtr2[j]);

                        // "broadcast" lower half to upper half
                        ownedJ = highway::ConcatLowerLower(tag_long, ownedJ, ownedJ);
                        x2 = highway::ConcatLowerLower(tag_double, x2, x2);
                        y2 = highway::ConcatLowerLower(tag_double, y2, y2);
                        z2 = highway::ConcatLowerLower(tag_double, z2, z2);

                        break;
                    }
                    }
                    case VectorizationPattern::pVecDiv2x2: {
                        // TODO : implement
                        break;
                    }
                    case VectorizationPattern::pVecx1: {
                        // TODO : implement
                        break;
                    }
                    case VectorizationPattern::pVecxVec: {
                        // TODO : implement
                        break;
                    }
                    default:
                        throw std::runtime_error("So the chosen vectorization pattern is not yet supported!");
                    }

                    // leave this as it is -> noPPL handles this different anyways
                    if (useMixing) {
                        double epsilon_buf[_vecLengthDouble] = {0};
                        double sigma_buf[_vecLengthDouble] = {0};
                        double shift_buf[_vecLengthDouble] = {0};

                        for(int n=0; (remainder ? n < restI : n < _vecLengthDouble); ++n) {
                            epsilon_buf[n] = _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + j));
                            sigma_buf[n] = _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + j));
                            if constexpr (applyShift) {
                                shift_buf[n] = _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + j));
                            }
                        }

                        epsilon24s = highway::LoadU(tag_double, epsilon_buf);
                        sigmaSquaredDiv2s = highway::LoadU(tag_double, sigma_buf);
                        shift6s = highway::LoadU(tag_double, epsilon_buf);
                    }
                    else {
                        epsilon24s = _epsilon24;
                        sigmaSquaredDiv2s = _sigmaSquared;
                        shift6s = _shift6;
                    }
                }

                template <bool remainder>
                inline void reduceAccumulatedForce(const VectorDouble& fxAcc, const VectorDouble& fyAcc, const VectorDouble& fzAcc,
                    double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr, const size_t i, const int rest = 0) {
                        
                    // TODO : think about more efficient way than for loop -> maybe highway provides something
                    
                    // TODO : handle case of rest

                    switch (_vectorizationPattern) {
                        case VectorizationPattern::p1xVec: {
                            fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
                            fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
                            fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
                            break;
                        }
                        case VectorizationPattern::p2xVecDiv2: {
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
                            break;
                        }
                        case VectorizationPattern::pVecDiv2x2: {
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

                            break;
                        }
                        case VectorizationPattern::pVecx1: {
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
                            break;
                        }
                        case VectorizationPattern::pVecxVec: {
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
                            break;
                        }
                    }
                }

                template <bool remainder>
                inline void handleNewton3Accumulation(const VectorDouble& fx, const VectorDouble& fy, const VectorDouble& fz,
                    double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr, const size_t j, const size_t rest = 0) {
                    
                    VectorDouble fx2;
                    VectorDouble fy2;
                    VectorDouble fz2;

                    auto mask = highway::FirstN(tag_double, rest);

                    switch (_vectorizationPattern)
                    {
                    case VectorizationPattern::p1xVec: {
                        // Load current force of j particles
                        if constexpr (!remainder) {
                            fx2 = highway::MaskedLoad(mask, tag_double, &fxPtr[j]);
                            fy2 = highway::MaskedLoad(mask, tag_double, &fyPtr[j]);
                            fz2 = highway::MaskedLoad(mask, tag_double, &fzPtr[j]);
                        }
                        else {
                            fx2 = highway::LoadU(tag_double, &fxPtr[j]);
                            fy2 = highway::LoadU(tag_double, &fyPtr[j]);
                            fz2 = highway::LoadU(tag_double, &fzPtr[j]);
                        }
                        // subtract force
                        fx2 = highway::Sub(fx2, fx);
                        fy2 = highway::Sub(fy2, fy);
                        fz2 = highway::Sub(fz2, fz);
                        // store force back into array
                        if constexpr (remainder) {
                            highway::BlendedStore(fx2, mask, tag_double, &fxPtr[j]);
                            highway::BlendedStore(fy2, mask, tag_double, &fyPtr[j]);
                            highway::BlendedStore(fz2, mask, tag_double, &fzPtr[j]);
                        }
                        else {
                            highway::StoreU(fx2, tag_double, &fxPtr[j]);
                            highway::StoreU(fy2, tag_double, &fyPtr[j]);
                            highway::StoreU(fz2, tag_double, &fzPtr[j]);
                        }
                        break;
                    }
                    case VectorizationPattern::p2xVecDiv2: {
                        // Load current force of j particles
                        mask = highway::And(mask, highway::FirstN(tag_double, _vecLengthDouble / 2));
                        fx2 = highway::MaskedLoad(mask, tag_double, &fxPtr[j]);
                        fy2 = highway::MaskedLoad(mask, tag_double, &fyPtr[j]);
                        fz2 = highway::MaskedLoad(mask, tag_double, &fzPtr[j]);

                        auto fx2_low = highway::LowerHalf(fx2);
                        auto fy2_low = highway::LowerHalf(fy2);
                        auto fz2_low = highway::LowerHalf(fz2);

                        auto fx_high = highway::UpperHalf(tag_double, fx);
                        auto fx_low = highway::LowerHalf(tag_double, fx);
                        auto fy_high = highway::UpperHalf(tag_double, fy);
                        auto fy_low = highway::LowerHalf(tag_double, fy);
                        auto fz_high = highway::UpperHalf(tag_double, fz);
                        auto fz_low = highway::LowerHalf(tag_double, fz);

                        auto fx_tmp = highway::Add(fx_high, fx_low);
                        auto fy_tmp = highway::Add(fy_high, fy_low);
                        auto fz_tmp = highway::Add(fz_high, fz_low);

                        // force substraction
                        fx2_low = highway::Sub(fx2_low, fx_tmp);
                        fy2_low = highway::Sub(fy2_low, fy_tmp);
                        fz2_low = highway::Sub(fz2_low, fz_tmp);

                        // TODO : figure out how to do this without combine and BlendedStore
                        fx2 = highway::Combine(tag_double, fx2_low, fx2_low);
                        fy2 = highway::Combine(tag_double, fy2_low, fy2_low);
                        fz2 = highway::Combine(tag_double, fz2_low, fz2_low);
                        
                        highway::BlendedStore(fx2, mask, tag_double, &fxPtr[j]);
                        highway::BlendedStore(fy2, mask, tag_double, &fyPtr[j]);
                        highway::BlendedStore(fz2, mask, tag_double, &fzPtr[j]);
                        
                        break;
                    }
                    case VectorizationPattern::pVecDiv2x2: {
                        // Handle j particle
                        auto fx_low = highway::LowerHalf(fx);
                        auto fy_low = highway::LowerHalf(fy);
                        auto fz_low = highway::LowerHalf(fz);

                        fxPtr[j] -= highway::ReduceSum(tag_double, fx_low);
                        fyPtr[j] -= highway::ReduceSum(tag_double, fy_low);
                        fzPtr[j] -= highway::ReduceSum(tag_double, fz_low);

                        // If no remainder, handle j+1 particle
                        if constexpr (!remainder) {
                            auto fx_high = highway::UpperHalf(tag_double, fx);
                            auto fy_high = highway::UpperHalf(tag_double, fy);
                            auto fz_high = highway::UpperHalf(tag_double, fz);

                            fxPtr[j+1] -= highway::ReduceSum(tag_double, fx_high);
                            fyPtr[j+1] -= highway::ReduceSum(tag_double, fy_high);
                            fzPtr[j+1] -= highway::ReduceSum(tag_double, fz_high);
                        }        
                        break;
                    }
                    case VectorizationPattern::pVecx1: {
                        
                        fxPtr[j] -= highway::ReduceSum(tag_double, fx);
                        fyPtr[j] -= highway::ReduceSum(tag_double, fy);
                        fzPtr[j] -= highway::ReduceSum(tag_double, fz);

                        break;
                    }
                    case VectorizationPattern::pVecxVec: {
                        // TODO : implement
                        break;
                    }
                    default: {
                        throw std::runtime_error("Unsupported Vectorization pattern!");
                    }
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

                    // loop over list for first time
                    for (size_t i = soa.size() - 1; (long)i >= 0; --i) {

                        if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
                            continue;
                        }

                        static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                            "OwnershipStates underlying type should be int64_t!");

                        VectorLong ownedStateI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr[i]));

                        VectorDouble fxAcc = highway::Zero(tag_double);
                        VectorDouble fyAcc = highway::Zero(tag_double);
                        VectorDouble fzAcc = highway::Zero(tag_double);

                        VectorDouble x1 = highway::Set(tag_double, xPtr[i]);
                        VectorDouble y1 = highway::Set(tag_double, yPtr[i]);
                        VectorDouble z1 = highway::Set(tag_double, zPtr[i]);

                        size_t j = 0;

                        for (; j < (i & ~(_vecLengthDouble - 1)); j+=_vecLengthDouble) {

                            SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1,
                                xPtr, yPtr, zPtr, fxPtr, fyPtr, fzPtr, &typeIDptr[i], typeIDptr, fxAcc, fyAcc, fzAcc,
                                virialSumX, virialSumY, virialSumZ, uPotSum, 0);
                        }

                        const int restJ = (int) (i & (_vecLengthDouble - 1));
                        if (restJ > 0) {
                            SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const int64_t*>(ownedStatePtr), x1, y1, z1,
                                xPtr, yPtr, zPtr, fxPtr, fyPtr, fzPtr, &typeIDptr[i], typeIDptr, fxAcc, fyAcc, fzAcc,
                                virialSumX, virialSumY, virialSumZ, uPotSum, restJ);
                        }

                        double sumFx = highway::ReduceSum(tag_double, fxAcc);
                        double sumFy = highway::ReduceSum(tag_double, fyAcc);
                        double sumFz = highway::ReduceSum(tag_double, fzAcc);

                        fxPtr[i] += sumFx;
                        fyPtr[i] += sumFy;
                        fzPtr[i] += sumFz;

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

                    for (unsigned int i = 0; i < soa1.size(); ++i) {

                        if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
                            continue;
                        }

                        VectorDouble fxAcc = _zeroDouble;
                        VectorDouble fyAcc = _zeroDouble;
                        VectorDouble fzAcc = _zeroDouble;

                        VectorLong ownedStateI = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr1[i]));

                        VectorDouble x1 = highway::Set(tag_double, x1Ptr[i]);
                        VectorDouble y1 = highway::Set(tag_double, y1Ptr[i]);
                        VectorDouble z1 = highway::Set(tag_double, z1Ptr[i]);

                        unsigned int j = 0;

                        for (; j < (soa2.size() & ~(_vecLengthDouble - 1)); j+=_vecLengthDouble) {
                            SoAKernel<newton3, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1,
                                x2Ptr, y2Ptr, z2Ptr, fx2Ptr, fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, fxAcc, fyAcc, fzAcc,
                                virialSumX, virialSumY, virialSumZ, uPotSum, 0);
                        }

                        const int restJ = (int) (soa2.size() & (_vecLengthDouble - 1));
                        if (restJ > 0) {
                            SoAKernel<newton3, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1,
                                x2Ptr, y2Ptr, z2Ptr, fx2Ptr, fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, fxAcc, fyAcc, fzAcc,
                                virialSumX, virialSumY, virialSumZ, uPotSum, restJ);
                        }

                        fx1Ptr[i] += highway::ReduceSum(tag_double, fxAcc);
                        fy1Ptr[i] += highway::ReduceSum(tag_double, fyAcc);
                        fz1Ptr[i] += highway::ReduceSum(tag_double, fzAcc);
                    }

                    if constexpr (calculateGlobals) {        
                        // nothing yet
                    }                
                }

                template <bool newton3, bool remainderIsMasked>
                inline void SoAKernel(const size_t j, const VectorLong& ownedStateI, const int64_t *const __restrict ownedStatePtr2,
                    const VectorDouble& x1, const VectorDouble& y1, const VectorDouble& z1, const double *const __restrict x2Ptr,
                    const double *const __restrict y2Ptr, const double *const __restrict z2Ptr, double *const __restrict fx2Ptr,
                    double *const __restrict fy2Ptr, double *const __restrict fz2Ptr, const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                    VectorDouble& fxAcc, VectorDouble& fyAcc, VectorDouble& fzAcc, VectorDouble& virialSumX, VectorDouble& virialSumY,
                    VectorDouble& virialSumZ, VectorDouble& uPotSum, const unsigned int rest = 0) {

                        // Prepare physical parameters
                        VectorDouble epsilon24s = _epsilon24;
                        VectorDouble sigmaSquareds = _sigmaSquared;
                        VectorDouble shift6s = _shift6;

                        if constexpr (useMixing) {
                            double epsilons [_vecLengthDouble];
                            double sigmas [_vecLengthDouble];
                            double shifts [_vecLengthDouble];

                            size_t n = 0;
                            for (; remainderIsMasked? n < rest : n < _vecLengthDouble; ++n) {
                                epsilons[n] = _PPLibrary->getMixing24Epsilon(*typeID1Ptr, *(typeID2Ptr + n));
                                sigmas[n] = _PPLibrary->getMixingSigmaSquared(*typeID1Ptr, *(typeID2Ptr + n));
                                if constexpr (applyShift) {
                                    shifts[n] = _PPLibrary->getMixingShift6(*typeID1Ptr, *(typeID2Ptr + n));
                                }
                            }

                            // Zero out rest
                            if constexpr (remainderIsMasked) {
                                for (; n < _vecLengthDouble; ++n) {
                                    epsilons[n] = 0.;
                                    sigmas[n] = 0.;
                                    if constexpr (applyShift) {
                                        shifts[n] = 0.;
                                    }
                                }
                            }

                            epsilon24s = highway::LoadU(tag_double, epsilons);
                            sigmaSquareds = highway::LoadU(tag_double, sigmas);
                            if constexpr (applyShift) {
                                shift6s = highway::LoadU(tag_double, shifts);
                            }
                        }

                        auto restMaskDouble = highway::FirstN(tag_double, rest);
                        auto restMaskLong = highway::FirstN(tag_long, rest);

                        // Load j particls
                        VectorDouble x2 = remainderIsMasked ? highway::MaskedLoad(restMaskDouble, tag_double, &x2Ptr[j]) : highway::LoadU(tag_double, &x2Ptr[j]);
                        VectorDouble y2 = remainderIsMasked ? highway::MaskedLoad(restMaskDouble, tag_double, &y2Ptr[j]) : highway::LoadU(tag_double, &y2Ptr[j]);
                        VectorDouble z2 = remainderIsMasked ? highway::MaskedLoad(restMaskDouble, tag_double, &z2Ptr[j]) : highway::LoadU(tag_double, &z2Ptr[j]);
                    
                        // distance calculations
                        auto drX = highway::Sub(x1, x2);
                        auto drY = highway::Sub(y1, y2);
                        auto drZ = highway::Sub(z1, z2);
                    
                        auto drX2 = highway::Mul(drX, drX);
                        auto drY2 = highway::Mul(drY, drY);
                        auto drZ2 = highway::Mul(drZ, drZ);

                        auto dr2 = highway::Add(drX2, highway::Add(drY2, drZ2));

                        auto cutoffMask = highway::Le(dr2, _cutoffSquared);

                        VectorLong ownedStateJ = remainderIsMasked ?
                            highway::MaskedLoad(restMaskLong, tag_long, &ownedStatePtr2[j])
                            : highway::LoadU(tag_long, &ownedStatePtr2[j]);

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
                        auto invDr2 = highway::Div(_oneDouble, dr2);
                        auto lj2 = highway::Mul(sigmaSquareds, invDr2);
                        auto lj4 = highway::Mul(lj2, lj2);
                        auto lj6 = highway::Mul(lj2, lj4);
                        auto lj12 = highway::Mul(lj6, lj6);
                        auto lj12m6 = highway::Sub(lj12, lj6);
                        auto lj12m6alj12 = highway::Add(lj12m6, lj12);
                        auto lj12m6alj12e = highway::Mul(lj12m6alj12, epsilon24s);
                        VectorDouble fac = highway::Mul(lj12m6alj12e, invDr2);

                        VectorDouble facMasked = remainderIsMasked ?
                            highway::IfThenElse(highway::And(restMaskDouble, cutoffDummyMask), fac, _zeroDouble)
                            : highway::IfThenElse(cutoffDummyMask, fac, _zeroDouble);

                        auto fx = highway::Mul(drX, facMasked);
                        auto fy = highway::Mul(drY, facMasked);
                        auto fz = highway::Mul(drZ, facMasked);

                        fxAcc = highway::Add(fxAcc, fx);
                        fyAcc = highway::Add(fyAcc, fy);
                        fzAcc = highway::Add(fzAcc, fz);

                        if constexpr (newton3) {

                            VectorDouble fx2 = remainderIsMasked ?
                                highway::MaskedLoad(restMaskDouble, tag_double, &fx2Ptr[j])
                                : highway::LoadU(tag_double, &fx2Ptr[j]);

                            VectorDouble fy2 = remainderIsMasked ?
                                highway::MaskedLoad(restMaskDouble, tag_double, &fy2Ptr[j])
                                : highway::LoadU(tag_double, &fy2Ptr[j]);

                            VectorDouble fz2 = remainderIsMasked ?
                                highway::MaskedLoad(restMaskDouble, tag_double, &fz2Ptr[j])
                                : highway::LoadU(tag_double, &fz2Ptr[j]);

                            auto fx2New = highway::Sub(fx2, fx);
                            auto fy2New = highway::Sub(fy2, fy);
                            auto fz2New = highway::Sub(fz2, fz);

                            if constexpr (remainderIsMasked) {
                                highway::BlendedStore(fx2New, restMaskDouble, tag_double, &fx2Ptr[j]);
                                highway::BlendedStore(fy2New, restMaskDouble, tag_double, &fy2Ptr[j]);
                                highway::BlendedStore(fz2New, restMaskDouble, tag_double, &fz2Ptr[j]);
                            }
                            else {
                                highway::StoreU(fx2New, tag_double, &fx2Ptr[j]);
                                highway::StoreU(fy2New, tag_double, &fy2Ptr[j]);
                                highway::StoreU(fz2New, tag_double, &fz2Ptr[j]);
                            }

                            if constexpr (calculateGlobals) {
                                // nothing yet xD
                            }
                        }
                    }

                template <bool newton3>
                inline std::tuple<VectorDouble, VectorDouble, VectorDouble> SoAKernel(const VectorLong& ownedStateI, const VectorLong& ownedStateJ,
                                        const VectorDouble& x1, const VectorDouble& y1, const VectorDouble& z1,
                                        const VectorDouble& x2, const VectorDouble& y2, const VectorDouble& z2,
                                        const VectorDouble& sigmaSquareds, const VectorDouble& epsilon24s,
                                        VectorDouble& virialSumX, VectorDouble& virialSumY, VectorDouble& virialSumZ,
                                        VectorDouble& potentialEnergySum) {

                    // determine distances
                    const VectorDouble drX = highway::Sub(x1, x2);
                    const VectorDouble drY = highway::Sub(y1, y2);
                    const VectorDouble drZ = highway::Sub(z1, z2);

                    const VectorDouble drX2 = highway::Mul(drX, drX);
                    const VectorDouble drY2 = highway::Mul(drY, drY);
                    const VectorDouble drZ2 = highway::Mul(drZ, drZ);
                    const VectorDouble distanceSquared = highway::Add(drX2, highway::Add(drY2, drZ2));

                    const auto cutoffMask = highway::Le(distanceSquared, _cutoffSquared);
                    const auto ownershipMask = highway::Ne(highway::ConvertTo(tag_double, ownedStateJ), highway::ConvertTo(tag_double, _ownedStateDummy));
                    const auto cutoffOwnershipMask = highway::And(cutoffMask, ownershipMask);

                    if (highway::AllFalse(tag_double, cutoffOwnershipMask)) {
                        return std::make_tuple(_zeroDouble, _zeroDouble, _zeroDouble);
                    }

                    // actual force computation
                    const VectorDouble invDistance2 = highway::Div(_oneDouble, distanceSquared);
                    const VectorDouble lj2 = highway::Mul(sigmaSquareds, invDistance2);
                    const VectorDouble lj4 = highway::Mul(lj2, lj2);
                    const VectorDouble lj6 = highway::Mul(lj2, lj4);
                    const VectorDouble lj12 = highway::Mul(lj6, lj6);
                    const VectorDouble lj12m6 = highway::Sub(lj12, lj6);
                    const VectorDouble lj12m6alj12 = highway::Add(lj12m6, lj12);
                    const VectorDouble lj12m6alj12e = highway::Mul(lj12m6alj12, epsilon24s);

                    const VectorDouble fAcc = highway::Mul(lj12m6alj12e, invDistance2);
                    const VectorDouble fAccMasked = highway::IfThenElse(cutoffOwnershipMask, fAcc, _zeroDouble);

                    const VectorDouble fx = highway::Mul(drX, fAccMasked);
                    const VectorDouble fy = highway::Mul(drY, fAccMasked);
                    const VectorDouble fz = highway::Mul(drZ, fAccMasked);

                    if constexpr (calculateGlobals) {
                        const VectorDouble virialX = highway::Mul(fx, drX);
                        const VectorDouble virialY = highway::Mul(fy, drY);
                        const VectorDouble virialZ = highway::Mul(fz, drZ);

                        const auto sigmaDivCutoffPow2 = highway::Div(sigmaSquareds, _cutoffSquared);
                        const auto sigmaDivCutoffPow4 = highway::Mul(sigmaDivCutoffPow2, sigmaDivCutoffPow2);
                        const auto sigmaDivCutoffPow6 = highway::Mul(sigmaDivCutoffPow2, sigmaDivCutoffPow4);
                        const auto sigmaDivCutoffPow12 = highway::Mul(sigmaDivCutoffPow6, sigmaDivCutoffPow6);
                        const auto sigmaDivCutoffPow6SubPow12 = highway::Sub(sigmaDivCutoffPow6, sigmaDivCutoffPow12);
                        const auto shift6s = applyShift and useMixing ? highway::Mul(epsilon24s, sigmaDivCutoffPow6SubPow12) : highway::Set(tag_double, _shift6AoS);

                        const VectorDouble uPotDot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
                        const VectorDouble uPotMasked = highway::IfThenElse(cutoffOwnershipMask, uPotDot, _zeroDouble);

                        auto ownedMaskI = highway::Eq(ownedStateI, _ownedStateOwned);
                        VectorDouble energyFactor = highway::ConvertTo(tag_double, highway::IfThenElse(ownedMaskI, _oneLong, _zeroLong));

                        if constexpr (newton3) {
                            auto ownedMaskJ = highway::Eq(ownedStateJ, _ownedStateOwned);
                            energyFactor = highway::Add(energyFactor, highway::ConvertTo(tag_double, highway::IfThenElse(ownedMaskJ, _oneLong, _zeroLong)));
                        }

                        energyFactor = highway::MulAdd(energyFactor, uPotMasked, energyFactor);
                        virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
                        virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
                        virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
                    }

                    return std::make_tuple(fx, fy, fz);

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
                                typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
                                ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
                            }

                            SoAKernel<newton3, false>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                                    x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                                    &typeIDPtr[indexFirst], typeID2Tmp.data(), fxAcc, fyAcc, fzAcc, virialSumX,
                                    virialSumY, virialSumZ, uPotSum, 0);

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
                                typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
                                ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
                            }

                            SoAKernel<newton3, true>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                                  x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                                  &typeIDPtr[indexFirst], typeID2Tmp.data(), fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                  virialSumZ, uPotSum, rest);

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

                const double _cutoffSquareAoS {0.};
                double _epsilon24AoS, _sigmaSquareAoS, _shift6AoS = 0.;
                ParticlePropertiesLibrary<double, size_t>* _PPLibrary = nullptr;
                double _uPotSum {0.};
                std::array<double, 3> _virialSum;
                std::vector<AoSThreadData> _aosThreadData;
                bool _postProcessed;

                VectorizationPattern _vectorizationPattern {VectorizationPattern::p1xVec};
    };
} // Highway
} // mdLib
HWY_AFTER_NAMESPACE();
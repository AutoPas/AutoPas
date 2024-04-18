// Created by Luis Gall on 27.03.2024

#pragma once

#include "MoleculeLJ_NoPPL.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/WrapOpenMP.h"
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
            explicit LJFunctorHWY(double cutoff)
                : autopas::Functor<molecule,
                LJFunctorHWY<applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(cutoff),
                _cutoffSquared{highway::Set(tag_double, cutoff*cutoff)},
                _cutoffSquareAoS{cutoff * cutoff},
                _uPotSum{0.},
                _virialSum{0.,0.,0.},
                _aosThreadData(),
                _postProcessed{false}
                {
                    if (calculateGlobals) {
                        _aosThreadData.resize(autopas::autopas_get_max_threads());
                    }

                    // AutoPasLog(INFO, "Highway Wrapper initialized with a register size of ({}) with architecture {}.", _vecLengthDouble, hwy::TargetName(HWY_TARGET));
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

                const auto sigmaMixed = useMixing ? i.getSigmaDiv2() + j.getSigmaDiv2() : 0;
                const auto sigmaSquared = useMixing ? sigmaMixed * sigmaMixed : _sigmaSquared;
                const auto epsilon24 = useMixing ? 24. * i.getSquareRootEpsilon() * j.getSquareRootEpsilon() : _epsilon24;

                auto dr = i.getR() - j.getR();
                double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

                if (dr2 > _cutoffSquareAoS) {
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
                    double uPot6 = epsilon24 * lj12m6;

                    if constexpr (applyShift) {
                        uPot6 += _shift6;
                    }

                    const int threadnum = autopas::autopas_get_thread_num();
                    // for non-newton3 the division is in the post-processing step.
                    if (newton3) {
                        uPot6 *= 0.5;
                        virial *= (double)0.5;
                    }
                    if (i.isOwned()) {
                        _aosThreadData[threadnum].uPotSum += uPot6;
                        _aosThreadData[threadnum].virialSum += virial;
                    }
                    // for non-newton3 the second particle will be considered in a separate calculation
                    if (newton3 and j.isOwned()) {
                        _aosThreadData[threadnum].uPotSum += uPot6;
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
                        i = i + _vecLengthDouble;
                        break;
                    
                    case VectorizationPattern::p2xVecDiv2:
                        i = i + _vecLengthDouble / 2;
                        break;
                    case VectorizationPattern::pVecDiv2x2:
                        i = i + 2;
                        break;
                    case VectorizationPattern::pVecx1:
                        ++i;
                        break;
                    case VectorizationPattern::pVecxVec:
                        i = i + _vecLengthDouble;
                        break;
                    default:
                        throw std::runtime_error("No vectorization pattern matched, error!");
                    }
                }

                template <bool remainder>
                inline void fillVectorRegisters(VectorDouble& x1, VectorDouble& y1, VectorDouble& z1,
                    VectorDouble& x2, VectorDouble& y2, VectorDouble& z2, VectorLong& ownedI, VectorLong& ownedJ,
                    VectorDouble& epsilon24s, VectorDouble& sigmaSquaredDiv2s, const double *const __restrict sigmaDiv2Ptr1,
                    const double *const __restrict sigmaDiv2Ptr2, const double *const __restrict sqrtEpsilonPtr1, const double *const __restrict sqrtEpsilonPtr2,
                    const double *const __restrict xPtr1, const double *const __restrict yPtr1, const double *const __restrict zPtr1,
                    const double *const __restrict xPtr2, const double *const __restrict yPtr2, const double *const __restrict zPtr2,
                    const int64_t *const __restrict ownedPtr1, const int64_t *const __restrict ownedPtr2, const size_t i, const size_t j, const size_t rest = 0) {

                    // TODO : handle various vectorization patterns with switch case
                    if constexpr (remainder) {
                        ownedI = highway::Set(tag_long, ownedPtr1[i]);
                        x1 = highway::Set(tag_double, xPtr1[i]);
                        y1 = highway::Set(tag_double, yPtr1[i]);
                        z1 = highway::Set(tag_double, zPtr1[i]);

                        const auto restMaskDouble = highway::FirstN(tag_double, rest);
                        const auto restMaskLong = highway::FirstN(tag_long, rest);

                        ownedJ = highway::MaskedLoad(restMaskLong, tag_long, &ownedPtr2[j]);
                        x2 = highway::MaskedLoad(restMaskDouble, tag_double, &xPtr2[j]);
                        y2 = highway::MaskedLoad(restMaskDouble, tag_double, &yPtr2[j]);
                        z2 = highway::MaskedLoad(restMaskDouble, tag_double, &zPtr2[j]);
                    }
                    else {
                        x2 = highway::LoadU(tag_double, &xPtr2[j]);
                        y2 = highway::LoadU(tag_double, &yPtr2[j]);
                        z2 = highway::LoadU(tag_double, &zPtr2[j]);
                    }

                    const auto sigmaDiv21 = highway::Set(tag_double, sigmaDiv2Ptr1[i]);
                    const auto sigmaDiv22 = highway::LoadU(tag_double, &sigmaDiv2Ptr2[j]);
                    const auto sigmaMixed = highway::Add(sigmaDiv21, sigmaDiv22);
                    sigmaSquaredDiv2s = useMixing ? highway::Mul(sigmaMixed, sigmaMixed) : highway::Set(tag_double, _sigmaSquared);

                    const auto sqrtEpsilon1Mul24 = highway::Set(tag_double, sqrtEpsilonPtr1[i] * 24.);
                    const auto sqrtEpsilon2 = highway::LoadU(tag_double, &sqrtEpsilonPtr2[j]);
                    epsilon24s = useMixing ? highway::Mul(sqrtEpsilon1Mul24, sqrtEpsilon2) : highway::Set(tag_double, _epsilon24);
                }

                inline void reduceAccumulatedForce(const VectorDouble& fxAcc, const VectorDouble& fyAcc, const VectorDouble& fzAcc,
                    double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr, const size_t i) {
                        
                    switch (_vectorizationPattern) {
                        case VectorizationPattern::p1xVec:
                            fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
                            fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
                            fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
                            break;
                        case VectorizationPattern::p2xVecDiv2:
                            double* tmpX [_vecLengthDouble] = {0.};
                            double* tmpY [_vecLengthDouble] = {0.};
                            double* tmpZ [_vecLengthDouble] = {0.};
                            highway::StoreU(fxAcc, tag_double, tmpX);
                            highway::StoreU(fyAcc, tag_double, tmpY);
                            highway::StoreU(fzAcc, tag_double, tmpZ);

                            for (size_t n = 0; n < _vecLengthDouble; ++n) {
                                fxPtr[i + n%2] += tmpX[n];
                                fyPtr[i + n%2] += tmpY[n];
                                fzPtr[i + n%2] += tmpZ[n];
                            }
                            break;
                        case VectorizationPattern::pVecDiv2x2:
                            double* tmpX [_vecLengthDouble] = {0.};
                            double* tmpY [_vecLengthDouble] = {0.};
                            double* tmpZ [_vecLengthDouble] = {0.};
                            highway::StoreU(fxAcc, tag_double, tmpX);
                            highway::StoreU(fyAcc, tag_double, tmpY);
                            highway::StoreU(fzAcc, tag_double, tmpZ);

                            for (size_t n = 0; n < _vecLengthDouble; ++n) {
                                fxPtr[i + n%(_vecLengthDouble/2)] += tmpX[n];
                                fyPtr[i + n%(_vecLengthDouble/2)] += tmpY[n];
                                fzPtr[i + n%(_vecLengthDouble/2)] += tmpZ[n];
                            }

                            break;
                        case VectorizationPattern::pVecx1:
                            double* tmpX [_vecLengthDouble] = {0.};
                            double* tmpY [_vecLengthDouble] = {0.};
                            double* tmpZ [_vecLengthDouble] = {0.};
                            highway::StoreU(fxAcc, tag_double, tmpX);
                            highway::StoreU(fyAcc, tag_double, tmpY);
                            highway::StoreU(fzAcc, tag_double, tmpZ);

                            for (size_t n = 0; n < _vecLengthDouble; ++n) {
                                fxPtr[i + n] += tmpX[n];
                                fyPtr[i + n] += tmpY[n];
                                fzPtr[i + n] += tmpZ[n];
                            }
                            break;
                        case VectorizationPattern::pVecxVec:
                            double* tmpX [_vecLengthDouble] = {0.};
                            double* tmpY [_vecLengthDouble] = {0.};
                            double* tmpZ [_vecLengthDouble] = {0.};
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

                template <bool remainder>
                inline void handleNewton3Accumulation(const VectorDouble& fx, const VectorDouble& fy, const VectorDouble& fz,
                    double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr, const size_t j, const size_t rest = 0) {
                        
                    // TODO : handle different vectorization patterns
                    VectorDouble fx2;
                    VectorDouble fy2;
                    VectorDouble fz2;

                    const auto restMask = highway::FirstN(tag_double, rest);

                    if constexpr (remainder) {
                        fx2 = highway::MaskedLoad(restMask, tag_double, &fxPtr[j]);
                        fy2 = highway::MaskedLoad(restMask, tag_double, &fyPtr[j]);
                        fz2 = highway::MaskedLoad(restMask, tag_double, &fzPtr[j]);
                    }
                    else {
                        fx2 = highway::LoadU(tag_double, &fxPtr[j]);
                        fy2 = highway::LoadU(tag_double, &fyPtr[j]);
                        fz2 = highway::LoadU(tag_double, &fzPtr[j]);
                    }

                    fx2 = highway::Sub(fx2, fx);
                    fy2 = highway::Sub(fy2, fy);
                    fz2 = highway::Sub(fz2, fz);

                    if constexpr (remainder) {
                        highway::BlendedStore(fx2, restMask, tag_double, &fxPtr[j]);
                        highway::BlendedStore(fy2, restMask, tag_double, &fyPtr[j]);
                        highway::BlendedStore(fz2, restMask, tag_double, &fzPtr[j]);
                    }
                    else {
                        highway::StoreU(fx2, tag_double, &fxPtr[j]);
                        highway::StoreU(fy2, tag_double, &fyPtr[j]);
                        highway::StoreU(fz2, tag_double, &fzPtr[j]);
                    }
                }

                template <bool newton3>
                inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {

                    if (soa.size() == 0) return;

                    // obtain iterators for the various values
                    const auto *const __restrict xPtr = soa.template begin<molecule::AttributeNames::posX>();
                    const auto *const __restrict yPtr = soa.template begin<molecule::AttributeNames::posY>();
                    const auto *const __restrict zPtr = soa.template begin<molecule::AttributeNames::posZ>();

                    const auto *const __restrict ownedStatePtr = reinterpret_cast<const int64_t *>(soa.template begin<molecule::AttributeNames::ownershipState>());

                    auto *const __restrict fxPtr = soa.template begin<molecule::AttributeNames::forceX>();
                    auto *const __restrict fyPtr = soa.template begin<molecule::AttributeNames::forceY>();
                    auto *const __restrict fzPtr = soa.template begin<molecule::AttributeNames::forceZ>();

                    const auto *const __restrict sigmaDiv2Ptr = soa.template begin<molecule::AttributeNames::sigmaDiv2>();
                    const auto *const __restrict sqrtEpsilonPtr = soa.template begin<molecule::AttributeNames::squareRootEpsilon>();

                    // initialize and declare vector variables
                    auto virialSumX = highway::Zero(tag_double);
                    auto virialSumY = highway::Zero(tag_double);
                    auto virialSumZ = highway::Zero(tag_double);
                    auto uPotSum = highway::Zero(tag_double);

                    VectorDouble x1;
                    VectorDouble y1;
                    VectorDouble z1;
                    VectorDouble x2;
                    VectorDouble y2;
                    VectorDouble z2;
                    VectorLong ownedStateI;
                    VectorLong ownedStateJ;

                    VectorDouble epsilon24s;
                    VectorDouble sigmaSquareds;

                    // loop over list for first time
                    for (size_t i = soa.size() - 1; (long)i >= 0; decrementFirstLoop(i)) {

                        VectorDouble fxAcc = highway::Zero(tag_double);
                        VectorDouble fyAcc = highway::Zero(tag_double);
                        VectorDouble fzAcc = highway::Zero(tag_double);

                        size_t j = 0;

                        // floor soa numParticles to multiple of vecLength
                        // If b is a power of 2 the following holds:
                        // a & ~(b -1) == a - (a mod b)
                        for (; j < (i & ~(_vecLengthDouble - 1)); incrementSecondLoop(j)) {

                            // load interaction partners
                            fillVectorRegisters<false>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, sigmaDiv2Ptr, sigmaDiv2Ptr, sqrtEpsilonPtr, sqrtEpsilonPtr,
                                xPtr, yPtr, zPtr, xPtr, yPtr, zPtr, ownedStatePtr, ownedStatePtr, i, j);

                            // compute forces
                            auto [fx, fy, fz] = SoAKernel<newton3>(ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, sigmaSquareds, epsilon24s, virialSumX, virialSumY, virialSumZ, uPotSum);
                            
                            // accumulate forces
                            fxAcc = highway::Add(fxAcc, fx);
                            fyAcc = highway::Add(fyAcc, fy);
                            fzAcc = highway::Add(fzAcc, fz);

                            if constexpr (newton3) {
                                handleNewton3Accumulation<false>(fx, fy, fz, fxPtr, fyPtr, fzPtr, j);
                            }
                        }

                        // If b is a power of 2 the following holds:
                        // a & (b -1) == a mod b
                        const int rest = (int)(i & (_vecLengthDouble - 1));
                        if (rest > 0) {

                            // load interaction partners
                            fillVectorRegisters<true>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ,
                                epsilon24s, sigmaSquareds, sigmaDiv2Ptr, sigmaDiv2Ptr, sqrtEpsilonPtr, sqrtEpsilonPtr,
                                xPtr, yPtr, zPtr, xPtr, yPtr, zPtr, ownedStatePtr, ownedStatePtr, i, j, rest);

                            // calculate forces
                            auto [fx, fy, fz] = SoAKernel<newton3>(ownedStateI, ownedStateJ, x1, y1, z1,
                                x2, y2, z2, sigmaSquareds, epsilon24s, virialSumX, virialSumY, virialSumZ, uPotSum);

                            // accumulate forces
                            fxAcc = highway::Add(fxAcc, fx);
                            fyAcc = highway::Add(fyAcc, fy);
                            fzAcc = highway::Add(fzAcc, fz);

                            if constexpr (newton3) {
                                handleNewton3Accumulation<true>(fx, fy, fz, fxPtr, fyPtr, fzPtr, j, rest);
                            }
                        }

                        reduceAccumulatedForce(fxAcc, fyAcc, fzAcc, fxPtr, fyPtr, fzPtr, i);
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

                    const auto *const __restrict sigmaDiv2Ptr1 = soa1.template begin<molecule::AttributeNames::sigmaDiv2>();
                    const auto *const __restrict sqrtEpsilonPtr1 = soa1.template begin<molecule::AttributeNames::squareRootEpsilon>();

                    const auto *const __restrict sigmaDiv2Ptr2 = soa2.template begin<molecule::AttributeNames::sigmaDiv2>();
                    const auto *const __restrict sqrtEpsilonPtr2 = soa2.template begin<molecule::AttributeNames::squareRootEpsilon>();

                    VectorDouble virialSumX = highway::Zero(tag_double);
                    VectorDouble virialSumY = highway::Zero(tag_double);
                    VectorDouble virialSumZ = highway::Zero(tag_double);
                    VectorDouble uPotSum = highway::Zero(tag_double);

                    VectorDouble x1;
                    VectorDouble y1;
                    VectorDouble z1;
                    VectorDouble x2;
                    VectorDouble y2;
                    VectorDouble z2;
                    VectorLong ownedStateI;
                    VectorLong ownedStateJ;

                    VectorDouble epsilon24s;
                    VectorDouble sigmaSquareds;

                    for (size_t i = 0; i < soa1.size(); incrementFirstLoop(i)) {
                        if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
                            continue;
                        }

                        VectorDouble fxAcc = highway::Zero(tag_double);
                        VectorDouble fyAcc = highway::Zero(tag_double);
                        VectorDouble fzAcc = highway::Zero(tag_double);

                        size_t j = 0;
                        for (; j < (soa2.size() & ~ (_vecLengthDouble - 1)); incrementSecondLoop(j)) {

                            fillVectorRegisters<false>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ, epsilon24s, sigmaSquareds,
                                sigmaDiv2Ptr1, sigmaDiv2Ptr2, sqrtEpsilonPtr1, sqrtEpsilonPtr2, x1Ptr, y1Ptr, z1Ptr, x2Ptr, y2Ptr, z2Ptr,
                                reinterpret_cast<const int64_t *>(ownedStatePtr1), reinterpret_cast<const int64_t *>(ownedStatePtr2), i, j);
                            
                            auto [fx, fy, fz] = SoAKernel<newton3>(ownedStateI, ownedStateJ, x1, y1, z1, x2, y2, z2,
                                sigmaSquareds, epsilon24s, virialSumX, virialSumY, virialSumZ, uPotSum);

                            fxAcc = highway::Add(fxAcc, fx);
                            fyAcc = highway::Add(fyAcc, fy);
                            fzAcc = highway::Add(fzAcc, fz);

                            if constexpr (newton3) {
                                handleNewton3Accumulation<false>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j);
                            }
                        }

                        const int rest = (int) (soa2.size() & (_vecLengthDouble - 1));
                        if (rest > 0) {

                            fillVectorRegisters<true>(x1, y1, z1, x2, y2, z2, ownedStateI, ownedStateJ, epsilon24s, sigmaSquareds,
                                sigmaDiv2Ptr1, sigmaDiv2Ptr2, sqrtEpsilonPtr1, sqrtEpsilonPtr2, x1Ptr, y1Ptr, z1Ptr, x2Ptr, y2Ptr, z2Ptr,
                                reinterpret_cast<const int64_t *>(ownedStatePtr1), reinterpret_cast<const int64_t *>(ownedStatePtr2), i, j, rest);
                            
                            auto [fx, fy, fz] = SoAKernel<newton3>(ownedStateI, ownedStateJ, x1, y1, z1, x2, y2, z2,
                                sigmaSquareds, epsilon24s, virialSumX, virialSumY, virialSumZ, uPotSum);

                            fxAcc = highway::Add(fxAcc, fx);
                            fyAcc = highway::Add(fyAcc, fy);
                            fzAcc = highway::Add(fzAcc, fz);

                            if constexpr (newton3) {
                                handleNewton3Accumulation<true>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, rest);
                            }
                        }

                        fx1Ptr[i] += highway::ReduceSum(tag_double, fxAcc);
                        fy1Ptr[i] += highway::ReduceSum(tag_double, fyAcc);
                        fz1Ptr[i] += highway::ReduceSum(tag_double, fzAcc);
                    }
                    if constexpr (calculateGlobals) {        
                        const int threadnum = autopas::autopas_get_thread_num();

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
                    if (ownedStatePtr[indexFirst] == (autopas::OwnershipState::dummy)) {
                        return;
                    }

                    const auto *const __restrict xPtr = soa.template begin<molecule::AttributeNames::posX>();
                    const auto *const __restrict yPtr = soa.template begin<molecule::AttributeNames::posY>();
                    const auto *const __restrict zPtr = soa.template begin<molecule::AttributeNames::posZ>();

                    auto *const __restrict fxPtr = soa.template begin<molecule::AttributeNames::forceX>();
                    auto *const __restrict fyPtr = soa.template begin<molecule::AttributeNames::forceY>();
                    auto *const __restrict fzPtr = soa.template begin<molecule::AttributeNames::forceZ>();

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

                        // TODO : reuse the generalized Kernel
                        // SoAKernel<newton3, false>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                        //         x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                        //         &typeIDPtr[indexFirst], typeID2Tmp.data(), fxAcc, fyAcc, fzAcc, virialSumX,
                        //         virialSumY, virialSumZ, uPotSum, 0);

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

                        // TODO : reuse the generalized Kernel
                        // SoAKernel<newton3, true>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                        //         x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                        //         &typeIDPtr[indexFirst], typeID2Tmp.data(), fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                        //         virialSumZ, uPotSum, rest);

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
                        const int threadnum = autopas::autopas_get_thread_num();

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
                    const auto ownershipMask = highway::Ne(highway::ConvertTo(tag_double, ownedStateJ), _ownedStateDummy);
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
                        const auto shift6s = applyShift and useMixing ? highway::Mul(epsilon24s, sigmaDivCutoffPow6SubPow12) : highway::Set(tag_double, _shift6);

                        const VectorDouble uPotDot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
                        const VectorDouble uPotMasked = highway::IfThenElse(cutoffOwnershipMask, uPotDot, _zeroDouble);

                        auto ownedMaskI = highway::Eq(ownedStateI, highway::ConvertTo(tag_long, _ownedStateOwned));
                        VectorDouble energyFactor = highway::ConvertTo(tag_double, highway::IfThenElse(ownedMaskI, _oneLong, _zeroLong));

                        if constexpr (newton3) {
                            auto ownedMaskJ = highway::Eq(ownedStateJ, highway::ConvertTo(tag_long, _ownedStateOwned));
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
                        const auto sigmaDivCutoffPow2 = sigmaSquare / _cutoffSquareAoS;
                        const auto sigmaDivCutoffPow6 = sigmaDivCutoffPow2 * sigmaDivCutoffPow2 * sigmaDivCutoffPow2;
                        _shift6 = epsilon24 * (sigmaDivCutoffPow6 - sigmaDivCutoffPow6 * sigmaDivCutoffPow6);
                    } else {
                        _shift6 = 0;
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
                const VectorDouble _ownedStateDummy{highway::ConvertTo(tag_double, highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::dummy)))};
                const VectorDouble _ownedStateOwned{highway::ConvertTo(tag_double, highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned)))};
                const VectorDouble _cutoffSquared {};

                const double _cutoffSquareAoS {0.};
                double _epsilon24, _sigmaSquared, _shift6 = 0.;
                double _uPotSum {0.};
                std::array<double, 3> _virialSum;
                std::vector<AoSThreadData> _aosThreadData;
                bool _postProcessed;

                VectorizationPattern _vectorizationPattern {VectorizationPattern::p1xVec};
    };
} // Highway
} // mdLib
HWY_AFTER_NAMESPACE();
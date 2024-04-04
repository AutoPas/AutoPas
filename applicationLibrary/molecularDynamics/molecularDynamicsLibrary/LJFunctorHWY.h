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

    namespace hn = hwy::HWY_NAMESPACE;

    // architecture specific information
    const hn::ScalableTag<double> tag_double;
    const hn::ScalableTag<int64_t> tag_long;
    const size_t _vecLengthDouble {hn::Lanes(tag_double)};
    typedef decltype(hn::Zero(tag_double)) VectorDouble;
    typedef decltype(hn::Zero(tag_long)) VectorLong;

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
                _cutoffSquared{cutoff * cutoff},
                _cutoffSquareAoS{cutoff * cutoff},
                _uPotSum{0.},
                _virialSum{0.,0.,0.},
                _aosThreadData(),
                _postProcessed{false} {
                    if (calculateGlobals) {
                        _aosThreadData.resize(autopas::autopas_get_max_threads());
                    }

                    AutoPasLog(INFO, "Highway Wrapper initialized with a register size of ({}) with architecture {}.", _vecLengthDouble, hwy::TargetName(HWY_TARGET));
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

                    // initialize vector variables
                    auto virialSumX = hn::Zero(tag_double);
                    auto virialSumY = hn::Zero(tag_double);
                    auto virialSumZ = hn::Zero(tag_double);
                    auto uPotSum = hn::Zero(tag_double);

                    // loop over list for first time
                    for (size_t i = soa.size() - 1; (long)i >= 0; --i) {
                        if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
                            continue;
                        }
                        // static_assert(std::is_same_v<std::underlying_type<autopas::OwnershipState>, int64_t>,
                        //      "OwnershipStates' underlying type should be int64_t");

                        VectorLong ownedStateI = hn::Set(tag_long, static_cast<int64_t>(ownedStatePtr[i]));
                        VectorDouble x1 = hn::Set(tag_double, xPtr[i]);
                        VectorDouble y1 = hn::Set(tag_double, yPtr[i]);
                        VectorDouble z1 = hn::Set(tag_double, yPtr[i]);

                        VectorDouble fxAcc = hn::Zero(tag_double);
                        VectorDouble fyAcc = hn::Zero(tag_double);
                        VectorDouble fzAcc = hn::Zero(tag_double);

                        size_t j = 0;

                        // floor soa numParticles to multiple of vecLength
                        // If b is a power of 2 the following holds:
                        // a & ~(b -1) == a - (a mod b)
                        for (; j < (i & ~(_vecLengthDouble - 1)); j += _vecLengthDouble) {
                            SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xPtr,
                                                yPtr, zPtr, fxPtr, fyPtr, fzPtr, &typeIDptr[i], typeIDptr, fxAcc, fyAcc, fzAcc,
                                                virialSumX, virialSumY, virialSumZ, uPotSum, 0);

                        }

                        // If b is a power of 2 the following holds:
                        // a & (b -1) == a mod b
                        const int rest = (int)(i & (_vecLengthDouble - 1));
                        if (rest > 0) {
                            SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xPtr,
                                                yPtr, zPtr, fxPtr, fyPtr, fzPtr, &typeIDptr[i], typeIDptr, fxAcc, fyAcc, fzAcc,
                                                virialSumX, virialSumY, virialSumZ, uPotSum, rest);
                        }

                        double sumFx = hn::ReduceSum(tag_double, fxAcc);
                        double sumFy = hn::ReduceSum(tag_double, fyAcc);
                        double sumFz = hn::ReduceSum(tag_double, fzAcc);

                        fxPtr[i] += sumFx;
                        fyPtr[i] += sumFy;
                        fzPtr[i] += sumFz;
                    }

                    if constexpr (calculateGlobals) {
                        const int threadnum = autopas::autopas_get_thread_num();

                        double globals[] = {
                            hn::ReduceSum(tag_double, virialSumX),
                            hn::ReduceSum(tag_double, virialSumY),
                            hn::ReduceSum(tag_double, virialSumZ),
                            hn::ReduceSum(tag_double, uPotSum)
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

                    // TODO : implement

                }

                template <bool newton3, bool remainderIsMasked>
                inline void SoAKernel(const size_t j, VectorLong &ownedStateI,
                                    const int64_t * __restrict owendStatePtr2,
                                    const VectorDouble& x1, const VectorDouble& y1, const VectorDouble& z1,
                                    const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
                                    const double *const __restrict z2Ptr, double *const __restrict fx2Ptr,
                                    double *const __restrict fy2Ptr, double *const __restrict fz2Ptr, const size_t *const __restrict typeID1ptr,
                                    const size_t *const __restrict typeID2ptr, VectorDouble& fxAcc, VectorDouble& fyAcc,
                                    VectorDouble& fzAcc, VectorDouble& virialSumX, VectorDouble& virialSumY, VectorDouble& virialSumZ,
                                    VectorDouble& potentialEnergySum, const unsigned int rest = 0) {

                    // prepare physical parameters (epsilon,...)
                    VectorDouble epsilon24s;
                    VectorDouble sigmaSquareds;
                    VectorDouble shift6s;

                    if (useMixing) {
                        double epsilon_buf[_vecLengthDouble] = {0};
                        double sigma_buf[_vecLengthDouble] = {0};
                        double shift_buf[_vecLengthDouble] = {0};

                        for(int i=0; (remainderIsMasked ? 0 /* TODO : figure out what */ : i < _vecLengthDouble); ++i) {
                            epsilon_buf[i] = _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + i));
                            sigma_buf[i] = _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + i));
                            if constexpr (applyShift) {
                                shift_buf[i] = _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + i));
                            }
                        }

                        epsilon24s = hn::LoadU(tag_double, epsilon_buf);
                        sigmaSquareds = hn::LoadU(tag_double, sigma_buf);
                        shift6s = hn::LoadU(tag_double, shift_buf);
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
                    if (remainderIsMasked) {
                        // TODO : implement
                    }
                    else {
                        x2 = hn::LoadU(tag_double, &x2Ptr[j]);
                        y2 = hn::LoadU(tag_double, &y2Ptr[j]);
                        z2 = hn::LoadU(tag_double, &z2Ptr[j]);
                    }

                    // perform distance calculations
                    const VectorDouble drx = hn::Sub(x1, x2);
                    const VectorDouble dry = hn::Sub(y1, y2);
                    const VectorDouble drz = hn::Sub(z1, z2);

                    const VectorDouble drx2 = hn::Mul(drx, drx);
                    const VectorDouble dry2 = hn::Mul(dry, dry);
                    const VectorDouble drz2 = hn::Mul(drz, drz);

                    const VectorDouble dr2_part = hn::Add(drx2, dry2);
                    const VectorDouble dr2 = hn::Add(dr2_part, drz2);

                    // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
                    // signaling = throw error if NaN is encountered
                    // dr2 <= _cutoffsquare ? 0xFFFFFFFFFFFFFFFF : 0
                    const auto cutoffMask = hn::Le(dr2, _cutoffSquared);
                    const VectorLong _zeroI = hn::ConvertTo(tag_long, _zero);

                    // TODO : continue with calculations
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
                    
                        AutoPasLog(INFO, "Running Highway with vector length ({}) on target {}", _vecLengthDouble, hwy::TargetName(HWY_TARGET));

                        // TODO : implement

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
                    _epsilon24 = epsilon24;
                    _sigmaSquared = sigmaSquare;
                    if constexpr (applyShift) {
                        _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, hn::GetLane(_cutoffSquared));
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
                inline VectorDouble wrapperFMA(const VectorDouble &factorA, const VectorDouble &factorB, const VectorDouble &summandC) {
                    return hn::MulAdd(factorA, factorB, summandC);
                }

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
                const VectorDouble _zero {hn::Zero(tag_double)};
                const VectorDouble _one  {hn::Set(tag_double, 1.)};
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
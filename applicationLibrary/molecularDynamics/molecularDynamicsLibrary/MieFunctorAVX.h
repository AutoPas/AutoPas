#ifndef AUTOPAS_MIEFUNCTORAVX2_H
#define AUTOPAS_MIEFUNCTORAVX2_H

#endif  // AUTOPAS_MIEFUNCTORAVX2_H
/**
  * @file MieFunctor
  *
  * @date 15 Sep 2023
  * @author K. Cole
  */
#pragma once
#ifndef __AVX__
#pragma message "Requested to compile MieFunctorAVX but AVX is not available!"
#else
#include <immintrin.h>
#endif


#include <array>
#include <cmath>

#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

    namespace mdLib {

    /**
         * A functor to handle mie interactions between two particles (molecules).
         * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
         * scheme.
         * This Version is implemented using AVX intrinsics.
         * @tparam Particle The type of particle.
         * @tparam ParticleCell The type of particlecell.
         * @tparam applyShift Switch for the Mie potential to be truncated shifted.
         * @tparam useMixing Switch for the functor to be used with multiple particle types.
         * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
         * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
         * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
         * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
     */
    template <class Particle, bool applyShift = false, bool useMixing = false,
              autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
              bool relevantForTuning = true>
    class MieFunctorAVX
        : public autopas::Functor<
              Particle, MieFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
      using SoAArraysType = typename Particle::SoAArraysType;

     private:
      /**
       * Computes and sets the DoubbleAdditionChain for th eexponents
       */
       void computeDoubleAdditionChain(){
         uint16_t  a,b;

         if(_mexpAoS>_nexpAoS){
           a = _nexpAoS;
           b = _mexpAoS;
         }
         else{
           a = _mexpAoS;
           b = _nexpAoS;
         }
         int chain = 0;
         int pointer = 1;
         int i = 0;
         while(a || b != 1) {
           i++;
           if (a <= b / 2) {
             i++;
             if (b & 1) {
               chain |= pointer;
             }
             pointer <<= 1;
             b = b / 2;
           } else {
             auto diff = b - a;
             b = a;
             a = diff;
             chain |= pointer;
           }
           pointer <<= 1;
         }

         int rev_chain = 0;

         for (int k = 0; k < i; k++,chain >>= 1) {
           rev_chain = (rev_chain<<1) | (chain & 1);
         }
         doubleAdditionChain = rev_chain;
         chain_len = i;
       }
     public:
      /**
       * Deleted default constructor
       */
      MieFunctorAVX() = delete;

      /**
       * Sets the particle properties constants for this functor.
       *
       * This is only necessary if no particlePropertiesLibrary is used.
       *
       * @param epsilon
       * @param sigmaSquared
       */
      void setParticleProperties(double epsilon, double sigmaSquared) {
         double c = static_cast<double>(_nexpAoS) / static_cast<double>(_nexpAoS - _mexpAoS);
         c = c * pow(static_cast<double>(_nexpAoS) / static_cast<double>(_mexpAoS), static_cast<double>(_mexpAoS) / (_nexpAoS - _mexpAoS));

         double cepsilon = c * epsilon;
         double cnepsilon = _nexpAoS * cepsilon;
         double cmepsilon = _mexpAoS * cepsilon;

#ifdef __AVX__
        _cepsilon = _mm256_set1_pd(cepsilon);
        _cnepsilon = _mm256_set1_pd(cnepsilon);
        _cmepsilon = _mm256_set1_pd(cmepsilon);
        _sigmaSquared = _mm256_set1_pd(sigmaSquared);


        if constexpr (applyShift) {
          _shift6 = _mm256_set1_pd(
              ParticlePropertiesLibrary<double, size_t>::calcShiftMie(cepsilon, sigmaSquared, _cutoffSquared[0], _nexpAoS,_mexpAoS));
        } else {
          _shift6 = _mm256_setzero_pd();
        }
#endif

        _cepsilonAoS = cepsilon;
        _cnepsilonAoS = cnepsilon;
        _cmepsilonAoS = cmepsilon;
        _sigmaSquaredAoS = sigmaSquared;

            if constexpr (applyShift) {
          _shift6AoS = ParticlePropertiesLibrary<double, size_t>::calcShiftMie(cepsilon, sigmaSquared, _cutoffSquaredAoS, _nexpAoS, _mexpAoS);
        } else {
          _shift6AoS = 0.0;
        }
      }

     private:
      /**
       * Internal, actual constructor.
       * @param cutoff
       * @param n_exp
       * @param m_exp
       * @note param dummy unused, only there to make the signature different from the public constructor.
       */
      explicit MieFunctorAVX(double cutoff, uint16_t n_exp, uint16_t m_exp, void * /*dummy*/)
    #ifdef __AVX__
          : autopas::Functor<
                Particle, MieFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
                cutoff),
            _cutoffSquared{_mm256_set1_pd(cutoff * cutoff)},
            _cutoffSquaredAoS(cutoff * cutoff),
            _nexp(_mm256_set1_pd(n_exp)),
            _nexpAoS(n_exp),
            _mexp(_mm256_set1_pd(m_exp)),
            _mexpAoS(m_exp),

            _potentialEnergySum{0.},
            _virialSum{0., 0., 0.},
            _aosThreadData(),
            _postProcessed{false} {
        if (calculateGlobals) {
          _aosThreadData.resize(autopas::autopas_get_max_threads());
        }
        computeDoubleAdditionChain();


      }

    #else
          : autopas::Functor<
                Particle, MieFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
                cutoff) {
        autopas::utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
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
       * @param n_exp
       * @param m_exp
       */
      explicit MieFunctorAVX(double cutoff, double n_exp, double m_exp) : MieFunctorAVX(cutoff, n_exp, m_exp, nullptr) {
        static_assert(not useMixing,
                      "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                      "mixing to false.");
      }

      /**
       * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
       * properties like sigma, epsilon and shift.
       * @param cutoff
       * @param n_exp
       * @param m_exp
       * @param particlePropertiesLibrary
       */
      explicit MieFunctorAVX(double cutoff, double n_exp, double m_exp,
                             ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
          : MieFunctorAVX(cutoff, n_exp, m_exp, nullptr) {
        static_assert(useMixing,
                      "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                      "or set mixing to true.");
        _PPLibrary = &particlePropertiesLibrary;
      }

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

      bool isRelevantForTuning() final { return relevantForTuning; }

      bool allowsNewton3() final {
        return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
      }

      bool allowsNonNewton3() final {
        return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
      }

      inline void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
        using namespace autopas::utils::ArrayMath::literals;
        if (i.isDummy() or j.isDummy()) {
          return;
        }
        auto sigmaSquared = _sigmaSquaredAoS;
        auto cepsilon = _cepsilonAoS;
        auto cnepsilon = _cnepsilonAoS;
        auto cmepsilon = _cmepsilonAoS;

        auto shift6 = _shift6AoS;
        auto m_exp = _mexpAoS;
        auto n_exp = _nexpAoS;

        if constexpr (useMixing) {
          sigmaSquared = _PPLibrary->getMixingSigmaSquaredMie(i.getTypeId(), j.getTypeId());
          cepsilon = _PPLibrary->getMixingCEpsilonMie(i.getTypeId(), j.getTypeId());
          cnepsilon = _PPLibrary->getMixingCNEpsilonMie(i.getTypeId(), j.getTypeId());
          cmepsilon = _PPLibrary->getMixingCMEpsilonMie(i.getTypeId(), j.getTypeId());

          if constexpr (applyShift) {
            shift6 = _PPLibrary->getMixingShift6Mie(i.getTypeId(), j.getTypeId());
          }
        }
        auto dr = i.getR() - j.getR();
        double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

        if (dr2 > _cutoffSquaredAoS) {
          return;
        }
        // for now inefficient power function



        double invdr2 = 1. / dr2;
        double fract = sigmaSquared * invdr2;

        double Mie_m = 1;
        double Mie_n = 1;

        if(mode==0) {

          Mie_m = m_exp & 1 ? sqrt(fract) : 1;
          Mie_n = n_exp & 1 ? sqrt(fract) : 1;

          for (size_t k = 1; k < m_exp; k += 2) {
            Mie_m = Mie_m * fract;
          }

          for (size_t k = 1; k < n_exp; k += 2) {
            Mie_n = Mie_n * fract;
          }
        }

        else if(mode==1) {

           if(m_exp&1) {
            Mie_m = sqrt(fract);
           }
           if(n_exp&1) {
            Mie_n = sqrt(fract);
           }

           m_exp>>=1;
           n_exp>>=1;

           while(m_exp||n_exp) {
            if(m_exp&1)
              Mie_m*=fract;
            if(n_exp&1)
              Mie_n*=fract;
            fract*=fract;
            m_exp>>=1;
            n_exp>>=1;
           }
        }

        else if(mode==2) {
           auto chain = doubleAdditionChain;
           auto base = sqrt(fract);
           double a = 1.0, b = base;
           for(uint16_t k=0;k<chain_len;k++,chain>>=1){

            if(chain&1){
              auto tmp = b*a;
              a = b;
              b = tmp;
            }
            else{
              chain>>=1;
              k++;
              b *=b;
              if(chain&1){
                b *=base;
              }
            }
           }
        Mie_m = a;
        Mie_n = b;

        }
        else if(mode==3) {
        }
           

       double Mie = (cnepsilon * Mie_n - cmepsilon * Mie_m);

        auto f = dr * (invdr2 * Mie);
        i.addF(f);
        if (newton3) {
          // only if we use newton 3 here, we want to
          j.subF(f);
        }
        if (calculateGlobals) {
          auto virial = dr * f;
          // Here we calculate either the potential energy * 6 or the potential energy * 12.
          // For newton3, this potential energy contribution is distributed evenly to the two molecules.
          // For non-newton3, the full potential energy is added to the one molecule.
          // The division by 6 is handled in endTraversal, as well as the division by two needed if newton3 is not used.
          double potentialEnergy6 = 6 * cepsilon * (Mie_n - Mie_m) + shift6;

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
       * @copydoc autopas::Functor::SoAFunctorSingle()
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
       * @copydoc autopas::Functor::SoAFunctorPair()
       */
      // clang-format on
      inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
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
      inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
    #ifdef __AVX__
        if (soa.getNumberOfParticles() == 0) return;

        const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
        const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
        const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

        const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

        auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
        auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
        auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

        const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

        __m256d virialSumX = _mm256_setzero_pd();
        __m256d virialSumY = _mm256_setzero_pd();
        __m256d virialSumZ = _mm256_setzero_pd();
        __m256d potentialEnergySum = _mm256_setzero_pd();

        // reverse outer loop s.th. inner loop always beginns at aligned array start
        // typecast to detect underflow
        for (size_t i = soa.getNumberOfParticles() - 1; (long)i >= 0; --i) {
          if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
            // If the i-th particle is a dummy, skip this loop iteration.
            continue;
          }

          static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                        "OwnershipStates underlying type should be int64_t!");
          // ownedStatePtr contains int64_t, so we broadcast these to make an __m256i.
          // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
          __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[i]));

          __m256d fxacc = _mm256_setzero_pd();
          __m256d fyacc = _mm256_setzero_pd();
          __m256d fzacc = _mm256_setzero_pd();

          const __m256d x1 = _mm256_broadcast_sd(&xptr[i]);
          const __m256d y1 = _mm256_broadcast_sd(&yptr[i]);
          const __m256d z1 = _mm256_broadcast_sd(&zptr[i]);

          size_t j = 0;
          // floor soa numParticles to multiple of vecLength
          // If b is a power of 2 the following holds:
          // a & ~(b -1) == a - (a mod b)
          for (; j < (i & ~(vecLength - 1)); j += 4) {
            SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                                   zptr, fxptr, fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc, &virialSumX,
                                   &virialSumY, &virialSumZ, &potentialEnergySum, 0);
          }
          // If b is a power of 2 the following holds:
          // a & (b -1) == a mod b
          const int rest = (int)(i & (vecLength - 1));
          if (rest > 0) {
            SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                                  zptr, fxptr, fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc, &virialSumX,
                                  &virialSumY, &virialSumZ, &potentialEnergySum, rest);
          }

          // horizontally reduce fDacc to sumfD
          const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
          const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

          const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
          const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

          const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
          const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

          const __m128d sumfxfyVEC = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
          const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

          const double sumfx = sumfxfyVEC[0];
          const double sumfy = sumfxfyVEC[1];
          const double sumfz = _mm_cvtsd_f64(sumfzVEC);

          fxptr[i] += sumfx;
          fyptr[i] += sumfy;
          fzptr[i] += sumfz;
        }

        if constexpr (calculateGlobals) {
          const int threadnum = autopas::autopas_get_thread_num();

          // horizontally reduce virialSumX and virialSumY
          const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
          const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
          const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
          const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

          // horizontally reduce virialSumZ and potentialEnergySum
          const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
          const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
          const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
          const __m128d hSumVirialzPotentialEnergyVec =
              _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

          // globals = {virialX, virialY, virialZ, potentialEnergy}
          double globals[4];
          _mm_store_pd(&globals[0], hSumVirialxyVec);
          _mm_store_pd(&globals[2], hSumVirialzPotentialEnergyVec);

          // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
          // false, since for newton3 disabled we divide by two later on.
          if constexpr (newton3) {
            _aosThreadData[threadnum].virialSumN3[0] += globals[0] * 0.5;
            _aosThreadData[threadnum].virialSumN3[1] += globals[1] * 0.5;
            _aosThreadData[threadnum].virialSumN3[2] += globals[2] * 0.5;
            _aosThreadData[threadnum].potentialEnergySumN3 += globals[3] * 0.5;
          } else {
            _aosThreadData[threadnum].virialSumNoN3[0] += globals[0];
            _aosThreadData[threadnum].virialSumNoN3[1] += globals[1];
            _aosThreadData[threadnum].virialSumNoN3[2] += globals[2];
            _aosThreadData[threadnum].potentialEnergySumNoN3 += globals[3];
          }
        }
    #endif
      }

      template <bool newton3>
      inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    #ifdef __AVX__
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

        __m256d virialSumX = _mm256_setzero_pd();
        __m256d virialSumY = _mm256_setzero_pd();
        __m256d virialSumZ = _mm256_setzero_pd();
        __m256d potentialEnergySum = _mm256_setzero_pd();

        for (unsigned int i = 0; i < soa1.getNumberOfParticles(); ++i) {
          if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
            // If the i-th particle is a dummy, skip this loop iteration.
            continue;
          }

          __m256d fxacc = _mm256_setzero_pd();
          __m256d fyacc = _mm256_setzero_pd();
          __m256d fzacc = _mm256_setzero_pd();

          static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                        "OwnershipStates underlying type should be int64_t!");
          // ownedStatePtr1 contains int64_t, so we broadcast these to make an __m256i.
          // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
          __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr1[i]));

          const __m256d x1 = _mm256_broadcast_sd(&x1ptr[i]);
          const __m256d y1 = _mm256_broadcast_sd(&y1ptr[i]);
          const __m256d z1 = _mm256_broadcast_sd(&z1ptr[i]);

          // floor soa2 numParticles to multiple of vecLength
          unsigned int j = 0;
          for (; j < (soa2.getNumberOfParticles() & ~(vecLength - 1)); j += 4) {
            SoAKernel<newton3, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                      y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                      &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, 0);
          }
          const int rest = (int)(soa2.getNumberOfParticles() & (vecLength - 1));
          if (rest > 0)
            SoAKernel<newton3, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                     y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                     &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, rest);

          // horizontally reduce fDacc to sumfD
          const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
          const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

          const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
          const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

          const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
          const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

          const __m128d sumfxfyVEC = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
          const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

          const double sumfx = sumfxfyVEC[0];
          const double sumfy = sumfxfyVEC[1];
          const double sumfz = _mm_cvtsd_f64(sumfzVEC);

          fx1ptr[i] += sumfx;
          fy1ptr[i] += sumfy;
          fz1ptr[i] += sumfz;
        }

        if constexpr (calculateGlobals) {
          const int threadnum = autopas::autopas_get_thread_num();

          // horizontally reduce virialSumX and virialSumY
          const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
          const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
          const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
          const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

          // horizontally reduce virialSumZ and potentialEnergySum
          const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
          const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
          const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
          const __m128d hSumVirialzPotentialEnergyVec =
              _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

          // globals = {virialX, virialY, virialZ, potentialEnergy}
          double globals[4];
          _mm_store_pd(&globals[0], hSumVirialxyVec);
          _mm_store_pd(&globals[2], hSumVirialzPotentialEnergyVec);

          // we have duplicated calculations, i.e., we calculate interactions multiple times, so we have to take care
          // that we do not add the energy multiple times!
          if constexpr (newton3) {
            // we count the energies partly to one of the two cells!
            _aosThreadData[threadnum].virialSumN3[0] += globals[0] * 0.5;
            _aosThreadData[threadnum].virialSumN3[1] += globals[1] * 0.5;
            _aosThreadData[threadnum].virialSumN3[2] += globals[2] * 0.5;
            _aosThreadData[threadnum].potentialEnergySumN3 += globals[3] * 0.5;
          } else {
            _aosThreadData[threadnum].virialSumNoN3[0] += globals[0];
            _aosThreadData[threadnum].virialSumNoN3[1] += globals[1];
            _aosThreadData[threadnum].virialSumNoN3[2] += globals[2];
            _aosThreadData[threadnum].potentialEnergySumNoN3 += globals[3];
          }
        }
    #endif
      }

    #ifdef __AVX__
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
       * @param potentialEnergySum
       * @param rest
       */
      template <bool newton3, bool remainderIsMasked>
      inline void SoAKernel(const size_t j, const __m256i ownedStateI, const int64_t *const __restrict ownedStatePtr2,
                            const __m256d &x1, const __m256d &y1, const __m256d &z1, const double *const __restrict x2ptr,
                            const double *const __restrict y2ptr, const double *const __restrict z2ptr,
                            double *const __restrict fx2ptr, double *const __restrict fy2ptr,
                            double *const __restrict fz2ptr, const size_t *const typeID1ptr, const size_t *const typeID2ptr,
                            __m256d &fxacc, __m256d &fyacc, __m256d &fzacc, __m256d *virialSumX, __m256d *virialSumY,
                            __m256d *virialSumZ, __m256d *potentialEnergySum, const unsigned int rest = 0) {
        __m256d cepsilons = _cepsilon;
        __m256d cnepsilons = _cnepsilon;
        __m256d cmepsilons = _cmepsilon;
        __m256d sigmaSquareds = _sigmaSquared;
        __m256d shift6s = _shift6;
        if (useMixing) {
          // the first argument for set lands in the last bits of the register
          cepsilons = _mm256_set_pd(
              not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingCEpsilonMie(*typeID1ptr, *(typeID2ptr + 3)) : 0,
              not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingCEpsilonMie(*typeID1ptr, *(typeID2ptr + 2)) : 0,
              not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingCEpsilonMie(*typeID1ptr, *(typeID2ptr + 1)) : 0,
              _PPLibrary->getMixingCEpsilonMie(*typeID1ptr, *(typeID2ptr + 0)));
          cnepsilons = _mm256_set_pd(
              not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingCNEpsilonMie(*typeID1ptr, *(typeID2ptr + 3)) : 0,
              not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingCNEpsilonMie(*typeID1ptr, *(typeID2ptr + 2)) : 0,
              not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingCNEpsilonMie(*typeID1ptr, *(typeID2ptr + 1)) : 0,
              _PPLibrary->getMixingCNEpsilonMie(*typeID1ptr, *(typeID2ptr + 0)));
          cmepsilons = _mm256_set_pd(
              not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingCMEpsilonMie(*typeID1ptr, *(typeID2ptr + 3)) : 0,
              not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingCMEpsilonMie(*typeID1ptr, *(typeID2ptr + 2)) : 0,
              not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingCMEpsilonMie(*typeID1ptr, *(typeID2ptr + 1)) : 0,
              _PPLibrary->getMixingCMEpsilonMie(*typeID1ptr, *(typeID2ptr + 0)));
          sigmaSquareds = _mm256_set_pd(
              not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingSigmaSquaredMie(*typeID1ptr, *(typeID2ptr + 3)) : 0,
              not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingSigmaSquaredMie(*typeID1ptr, *(typeID2ptr + 2)) : 0,
              not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingSigmaSquaredMie(*typeID1ptr, *(typeID2ptr + 1)) : 0,
              _PPLibrary->getMixingSigmaSquaredMie(*typeID1ptr, *(typeID2ptr + 0)));
          if constexpr (applyShift) {
            shift6s = _mm256_set_pd(
                (not remainderIsMasked or rest > 3) ? _PPLibrary->getMixingShift6Mie(*typeID1ptr, *(typeID2ptr + 3)) : 0,
                (not remainderIsMasked or rest > 2) ? _PPLibrary->getMixingShift6Mie(*typeID1ptr, *(typeID2ptr + 2)) : 0,
                (not remainderIsMasked or rest > 1) ? _PPLibrary->getMixingShift6Mie(*typeID1ptr, *(typeID2ptr + 1)) : 0,
                _PPLibrary->getMixingShift6Mie(*typeID1ptr, *(typeID2ptr + 0)));
          }
        }

        const __m256d x2 = remainderIsMasked ? _mm256_maskload_pd(&x2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&x2ptr[j]);
        const __m256d y2 = remainderIsMasked ? _mm256_maskload_pd(&y2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&y2ptr[j]);
        const __m256d z2 = remainderIsMasked ? _mm256_maskload_pd(&z2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&z2ptr[j]);

        const __m256d drx = _mm256_sub_pd(x1, x2);
        const __m256d dry = _mm256_sub_pd(y1, y2);
        const __m256d drz = _mm256_sub_pd(z1, z2);

        const __m256d drx2 = _mm256_mul_pd(drx, drx);
        const __m256d dry2 = _mm256_mul_pd(dry, dry);
        const __m256d drz2 = _mm256_mul_pd(drz, drz);

        const __m256d dr2PART = _mm256_add_pd(drx2, dry2);
        const __m256d dr2 = _mm256_add_pd(dr2PART, drz2);

        // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
        // signaling = throw error if NaN is encountered
        // dr2 <= _cutoffSquared ? 0xFFFFFFFFFFFFFFFF : 0
        const __m256d cutoffMask = _mm256_cmp_pd(dr2, _cutoffSquared, _CMP_LE_OS);

        // This requires that dummy is zero (otherwise when loading using a mask the owned state will not be zero)
        const __m256i ownedStateJ = remainderIsMasked
                                        ? _mm256_castpd_si256(_mm256_maskload_pd(
                                              reinterpret_cast<double const *>(&ownedStatePtr2[j]), _masks[rest - 1]))
                                        : _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&ownedStatePtr2[j]));
        // This requires that dummy is the first entry in OwnershipState!
        const __m256d dummyMask = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateJ), _zero, _CMP_NEQ_UQ);
        const __m256d cutoffDummyMask = _mm256_and_pd(cutoffMask, dummyMask);

        // if everything is masked away return from this function.
        if (_mm256_movemask_pd(cutoffDummyMask) == 0) {
          return;
        }
        uint16_t m_exp = _mexpAoS;
        uint16_t n_exp = _nexpAoS;
        const __m256d invdr2 = _mm256_div_pd(_one, dr2);
        __m256d miem, mien;
        const __m256d mie2 = _mm256_mul_pd(sigmaSquareds, invdr2);

        if(mode==0) {
          miem = m_exp & 1 ? _mm256_sqrt_pd(mie2) : _one;
          mien = n_exp & 1 ? _mm256_sqrt_pd(mie2) : _one;

          for (size_t k = 1; k < m_exp; k += 2) {
            miem = _mm256_mul_pd(miem, mie2);
          }
          for (size_t k = 1; k < n_exp; k += 2) {
            mien = _mm256_mul_pd(mien, mie2);
          }
        }
        else if(mode==1) {

          miem = m_exp & 1 ? _mm256_sqrt_pd(mie2) : _one;
          mien = n_exp & 1 ? _mm256_sqrt_pd(mie2) : _one;

          m_exp>>=1;
          n_exp>>=1;
          __m256d fact = mie2;
          while (m_exp||n_exp) {
            if (m_exp & 1) {
              miem = _mm256_mul_pd(miem,fact);
            }
            if (n_exp & 1) {
              mien = _mm256_mul_pd(mien,fact);
            }
            fact = _mm256_mul_pd(fact,fact);
            m_exp >>= 1;
            n_exp >>= 1;
          }
        }
        else if(mode==2) {
          auto chain = doubleAdditionChain;
          __m256d base =  _mm256_sqrt_pd(mie2);
          __m256d a = _one, b = base;
          for(size_t k=0;k<chain_len;k++,chain>>=1){

            if(chain&1){
              auto tmp = b*a;
              a = b;
              b = tmp;
            }
            else{
              chain>>=1;
              k++;
              b *=b;
              if(chain&1){
                b *=base;
              }
            }
          }
          miem = a;
          mien = b;

        }
        const __m256d nmien = _mm256_mul_pd(_cnepsilon, mien);
        const __m256d mmiem = _mm256_mul_pd(_cmepsilon, miem);
        const __m256d miensubmiem = _mm256_sub_pd(nmien, mmiem);
        const __m256d fac = _mm256_mul_pd(miensubmiem, invdr2);

        const __m256d facMasked =
            remainderIsMasked ? _mm256_and_pd(fac, _mm256_and_pd(cutoffDummyMask, _mm256_castsi256_pd(_masks[rest - 1])))
                              : _mm256_and_pd(fac, cutoffDummyMask);

        const __m256d fx = _mm256_mul_pd(drx, facMasked);
        const __m256d fy = _mm256_mul_pd(dry, facMasked);
        const __m256d fz = _mm256_mul_pd(drz, facMasked);

        fxacc = _mm256_add_pd(fxacc, fx);
        fyacc = _mm256_add_pd(fyacc, fy);
        fzacc = _mm256_add_pd(fzacc, fz);

        // if newton 3 is used subtract fD from particle j
        if constexpr (newton3) {
          const __m256d fx2 =
              remainderIsMasked ? _mm256_maskload_pd(&fx2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fx2ptr[j]);
          const __m256d fy2 =
              remainderIsMasked ? _mm256_maskload_pd(&fy2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fy2ptr[j]);
          const __m256d fz2 =
              remainderIsMasked ? _mm256_maskload_pd(&fz2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fz2ptr[j]);

          const __m256d fx2new = _mm256_sub_pd(fx2, fx);
          const __m256d fy2new = _mm256_sub_pd(fy2, fy);
          const __m256d fz2new = _mm256_sub_pd(fz2, fz);

          remainderIsMasked ? _mm256_maskstore_pd(&fx2ptr[j], _masks[rest - 1], fx2new)
                            : _mm256_storeu_pd(&fx2ptr[j], fx2new);
          remainderIsMasked ? _mm256_maskstore_pd(&fy2ptr[j], _masks[rest - 1], fy2new)
                            : _mm256_storeu_pd(&fy2ptr[j], fy2new);
          remainderIsMasked ? _mm256_maskstore_pd(&fz2ptr[j], _masks[rest - 1], fz2new)
                            : _mm256_storeu_pd(&fz2ptr[j], fz2new);
        }

        if constexpr (calculateGlobals) {
          // Global Virial
          const __m256d virialX = _mm256_mul_pd(fx, drx);
          const __m256d virialY = _mm256_mul_pd(fy, dry);
          const __m256d virialZ = _mm256_mul_pd(fz, drz);
          // Global Potential
          //ToDo: Maybe/store _cepsilon*_six
          const __m256d potentialEnergy = wrapperFMA(_mm256_mul_pd(_cepsilon, _six), _mm256_sub_pd(mien, miem), shift6s);

          const __m256d potentialEnergyMasked =
              remainderIsMasked
                  ? _mm256_and_pd(potentialEnergy, _mm256_and_pd(cutoffDummyMask, _mm256_castsi256_pd(_masks[rest - 1])))
                  : _mm256_and_pd(potentialEnergy, cutoffDummyMask);

          __m256d ownedMaskI =
              _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateI), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
          __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskI);
          if constexpr (newton3) {
            __m256d ownedMaskJ =
                _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateJ), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
            energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, ownedMaskJ));
          }
          *potentialEnergySum = wrapperFMA(energyFactor, potentialEnergyMasked, *potentialEnergySum);
          *virialSumX = wrapperFMA(energyFactor, virialX, *virialSumX);
          *virialSumY = wrapperFMA(energyFactor, virialY, *virialSumY);
          *virialSumZ = wrapperFMA(energyFactor, virialZ, *virialSumZ);
        }
      }
    #endif

           public:
            // clang-format off
        /**
         * @copydoc autopas::Functor::SoAFunctorVerlet()
         * @note If you want to parallelize this by openmp, please ensure that there
         * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
         */
            // clang-format on
            inline void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
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
            inline void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                                             const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
      #ifdef __AVX__
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
              __m256d virialSumX = _mm256_setzero_pd();
              __m256d virialSumY = _mm256_setzero_pd();
              __m256d virialSumZ = _mm256_setzero_pd();
              __m256d potentialEnergySum = _mm256_setzero_pd();
              __m256d fxacc = _mm256_setzero_pd();
              __m256d fyacc = _mm256_setzero_pd();
              __m256d fzacc = _mm256_setzero_pd();

              // broadcast particle 1
              const auto x1 = _mm256_broadcast_sd(&xptr[indexFirst]);
              const auto y1 = _mm256_broadcast_sd(&yptr[indexFirst]);
              const auto z1 = _mm256_broadcast_sd(&zptr[indexFirst]);
              // ownedStatePtr contains int64_t, so we broadcast these to make an __m256i.
              // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
              __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[indexFirst]));

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
                                          &virialSumY, &virialSumZ, &potentialEnergySum, 0);

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
                                         &virialSumZ, &potentialEnergySum, rest);

                if constexpr (newton3) {
                  for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
                    fxptr[neighborList[j + vecIndex]] = fx2tmp[vecIndex];
                    fyptr[neighborList[j + vecIndex]] = fy2tmp[vecIndex];
                    fzptr[neighborList[j + vecIndex]] = fz2tmp[vecIndex];
                  }
                }
              }

              // horizontally reduce fDacc to sumfD
              const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
              const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

              const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
              const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

              const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
              const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

              const __m128d sumfxfyVEC = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
              const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

              const double sumfx = sumfxfyVEC[0];
              const double sumfy = sumfxfyVEC[1];
              const double sumfz = _mm_cvtsd_f64(sumfzVEC);

              fxptr[indexFirst] += sumfx;
              fyptr[indexFirst] += sumfy;
              fzptr[indexFirst] += sumfz;

              if constexpr (calculateGlobals) {
                const int threadnum = autopas::autopas_get_thread_num();

                // horizontally reduce virialSumX and virialSumY
                const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
                const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
                const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
                const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

                // horizontally reduce virialSumZ and potentialEnergySum
                const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
                const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
                const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
                const __m128d hSumVirialzPotentialEnergyVec =
                    _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

                // globals = {virialX, virialY, virialZ, potentialEnergy}
                double globals[4];
                _mm_store_pd(&globals[0], hSumVirialxyVec);
                _mm_store_pd(&globals[2], hSumVirialzPotentialEnergyVec);

                // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
                // false, since for newton3 disabled we divide by two later on.
                if (newton3) {
                  _aosThreadData[threadnum].virialSumN3[0] += globals[0] * 0.5;
                  _aosThreadData[threadnum].virialSumN3[1] += globals[1] * 0.5;
                  _aosThreadData[threadnum].virialSumN3[2] += globals[2] * 0.5;
                  _aosThreadData[threadnum].potentialEnergySumN3 += globals[3] * 0.5;
                } else {
                  _aosThreadData[threadnum].virialSumNoN3[0] += globals[0];
                  _aosThreadData[threadnum].virialSumNoN3[1] += globals[1];
                  _aosThreadData[threadnum].virialSumNoN3[2] += globals[2];
                  _aosThreadData[threadnum].potentialEnergySumNoN3 += globals[3];
                }
              }
              // interact with i with 4 neighbors
      #endif  // __AVX__
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
   * @note molAType and molBType make no difference for MieFunctor, but are kept to have a consistent interface for other
   * functors where they may.
   * @return the number of floating point operations
       */
      static unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, bool newton3) {
        // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply scale) sum
        // Adding to particle forces: 6 or 3 depending newton3
        // Total = 12 + (6 or 3) = 18 or 15
        return newton3 ? 18ul : 15ul;
      }
     public:
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
   * Sets the Functor properties constants for this functor. Updates the epsilon value
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param n_exp
   * @param m_exp
       */
      void setFunctorProperties(double n_exp, double m_exp) {
#ifdef __AVX__
        _nexp = _mm256_set1_pd(n_exp);
        _mexp = _mm256_set1_pd(m_exp);

        if constexpr (applyShift) {
        } else {
          _shift6 = _mm256_setzero_pd();
        }
#endif
        _nexpAoS = n_exp;
        _mexpAoS = m_exp;
        if constexpr (applyShift) {
        } else {
          _shift6AoS = 0.;
        }
      }

     private:
#ifdef __AVX__
      /**
   * Wrapper function for FMA. If FMA is not supported it executes first the multiplication then the addition.
   * @param factorA
   * @param factorB
   * @param summandC
   * @return A * B + C
       */
      inline __m256d wrapperFMA(const __m256d &factorA, const __m256d &factorB, const __m256d &summandC) {
#ifdef __FMA__
        return _mm256_fmadd_pd(factorA, factorB, summandC);
#elif __AVX__
        const __m256d tmp = _mm256_mul_pd(factorA, factorB);
        return _mm256_add_pd(summandC, tmp);
#else
        // dummy return. If no vectorization is available this whole class is pointless anyways.
        return __m256d();
#endif
      }
#endif


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

  #ifdef __AVX__
    const __m256d _zero{_mm256_set1_pd(0.)};
    const __m256d _one{_mm256_set1_pd(1.)};
    const __m256d _six{_mm256_set1_pd(6.)};
    const __m256i _vindex = _mm256_set_epi64x(0, 1, 3, 4);
    const __m256i _masks[3]{
        _mm256_set_epi64x(0, 0, 0, -1),
        _mm256_set_epi64x(0, 0, -1, -1),
        _mm256_set_epi64x(0, -1, -1, -1),
    };
    const __m256i _ownedStateDummyMM256i{0x0};
    const __m256i _ownedStateOwnedMM256i{_mm256_set1_epi64x(static_cast<int64_t>(autopas::OwnershipState::owned))};
    const __m256d _cutoffSquared{};
    __m256d _nexp{};
    __m256d _mexp{};
    __m256d _shift6 = _mm256_setzero_pd();
    __m256d _cepsilon{};
    __m256d _cnepsilon{};
    __m256d _cmepsilon{};
    __m256d _sigmaSquared{};


  #endif
    const uint8_t mode = 1;
    const double _cutoffSquaredAoS = 0;
    double _cepsilonAoS, _cnepsilonAoS, _cmepsilonAoS, _sigmaSquaredAoS, _shift6AoS = 0;
    uint16_t _nexpAoS,_mexpAoS = 0;
    int doubleAdditionChain;
    uint16_t chain_len;

    ParticlePropertiesLibrary<double, size_t> *_PPLibrary = nullptr;

    // sum of the potential energy, only calculated if calculateGlobals is true
    double _potentialEnergySum;

    // sum of the virial, only calculated if calculateGlobals is true
    std::array<double, 3> _virialSum;

    // thread buffer for aos
    std::vector<AoSThreadData> _aosThreadData;

    // defines whether or whether not the global values are already preprocessed
    bool _postProcessed;

    // number of double values that fit into a vector register.
    // MUST be power of 2 because some optimizations make this assumption
    constexpr static size_t vecLength = 4;
    };
  }
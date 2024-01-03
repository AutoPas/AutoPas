#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
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
class MieFunctorFixed
    : public autopas::Functor<
          Particle, MieFunctorFixed<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  MieFunctorFixed() = delete;

 private:
  /**
       * Computes and sets the DoubbleAdditionChain for the exponents
   */
  void initDoubleAdditionChain(){
    uint16_t  a,b;
    a = (_mexp > _nexp) ? _nexp : _mexp;
    b = (_mexp > _nexp) ? _mexp : _nexp;

    int chain = 0;
    int pointer = 1;
    int i = 0;

    for(;a || b != 1;i++,pointer <<= 1) {
      if (a <= b / 2) {
        i++;
        if (b & 1) {
          chain |= pointer;
        }
        pointer <<= 1;
        b >>= 1 ;
      } else {

        auto diff = b - a;
        b = a;
        a = diff;
        chain |= pointer;
      }
    }

    int rev_chain = 0;
    for (int k = 0; k < i; k++,chain >>= 1) {
      rev_chain = (rev_chain<<1) | (chain & 1);
    }
    doubleAdditionChain = rev_chain;
    chain_len = i;
  }
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit MieFunctorFixed(double cutoff, void * /*dummy*/)
      : autopas::Functor<Particle,
                         MieFunctorFixed<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},

        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (mode==2) {
      initDoubleAdditionChain();
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
  explicit MieFunctorFixed(double cutoff) : MieFunctorFixed(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like sigma, epsilon and shift.
   * properties like sigma, epsilon and shift.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit MieFunctorFixed(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : MieFunctorFixed(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
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
    auto sigmaSquared = _sigmaSquared;
    auto cepsilon = _cepsilon;
    auto cnepsilon = _cnepsilon;
    auto cmepsilon = _cmepsilon;

    auto shift6 = _shift6;
    auto m_exp = _mexp;
    auto n_exp = _nexp;

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

    if (dr2 > _cutoffSquared) {
      return;
    }

    double invdr2 = 1. / dr2;
    double fract = sigmaSquared * invdr2;

    double Mie_m = 1;
    double Mie_n = 1;

      if constexpr(mode==0) {
      Mie_m = m_exp & 1 ? sqrt(fract) : 1;
      Mie_n = n_exp & 1 ? sqrt(fract) : 1;

      for (size_t k = 1; k < m_exp; k += 2) {
        Mie_m = Mie_m * fract;
      }
      for (size_t k = 1; k < n_exp; k += 2) {
        Mie_n = Mie_n * fract;
      }
    }
      else if constexpr (mode==1) {
      if(m_exp&1) {
        Mie_m = sqrt(fract);
      }
      if(n_exp&1) {
        Mie_n = sqrt(fract);
      }

      m_exp>>=1;
      n_exp>>=1;
      while(m_exp || n_exp) {
        if(m_exp&1)
          Mie_m*=fract;
        if(n_exp&1)
          Mie_n*=fract;
        fract*=fract;
        m_exp>>=1;
        n_exp>>=1;
      }
      }

      else if constexpr (mode==2) {
      auto chain = doubleAdditionChain;
      auto base = sqrt(fract);
      double a = 1.0, b = base;

      for(uint16_t k=0; k < chain_len; k++,chain>>=1){
        if(chain&1){
          std::swap(a,b);
          b *= a;
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
   * This functor will always use a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.getNumberOfParticles() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();
    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum =
        0.;  // Note: this is not the potential energy, but some fixed multiple of it.
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> cepsilons;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> cnepsilons;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> cmepsilons;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;
    if constexpr (useMixing) {
      // Preload all sigma and epsilons for next vectorized region.
      // Not preloading and directly using the values, will produce worse results.
      sigmaSquareds.resize(soa.getNumberOfParticles());
      cepsilons.resize(soa.getNumberOfParticles());
      cnepsilons.resize(soa.getNumberOfParticles());
      cmepsilons.resize(soa.getNumberOfParticles());
      // if no mixing or mixing but no shift shift6 is constant therefore we do not need this vector.
      if constexpr (applyShift) {
        shift6s.resize(soa.getNumberOfParticles());
      }
    }

    const SoAFloatPrecision const_shift6 = _shift6;
    const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    const SoAFloatPrecision const_cepsilon = _cepsilon;
    const SoAFloatPrecision const_cnepsilon = _cnepsilon;
    const SoAFloatPrecision const_cmepsilon = _cmepsilon;


    for (unsigned int i = 0; i < soa.getNumberOfParticles(); ++i) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      if constexpr (useMixing) {
        for (unsigned int j = 0; j < soa.getNumberOfParticles(); ++j) {
          auto mixingData = _PPLibrary->getMixingDataMie(typeptr[i], typeptr[j]);
          sigmaSquareds[j] = mixingData.sigmaSquared;
          cepsilons[j] = mixingData.cepsilon;
          cnepsilons[j] = mixingData.cnepsilon;
          cmepsilons[j] = mixingData.cmepsilon;

          if constexpr (applyShift) {
            shift6s[j] = mixingData.shift6;
          }
        }
      }

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
      for (unsigned int j = i + 1; j < soa.getNumberOfParticles(); ++j) {
        SoAFloatPrecision shift6 = const_shift6;
        SoAFloatPrecision sigmaSquared = const_sigmaSquared;
        SoAFloatPrecision cepsilon = const_cepsilon;
        SoAFloatPrecision cnepsilon = const_cnepsilon;
        SoAFloatPrecision cmepsilon = const_cmepsilon;
        if constexpr (useMixing) {
          sigmaSquared = sigmaSquareds[j];
          cepsilon = cepsilons[j];
          cnepsilon = cnepsilons[j];
          cmepsilon = cmepsilons[j];
          if constexpr (applyShift) {
            shift6 = shift6s[j];
          }
        }

        const auto ownedStateJ = ownedStatePtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

        // Mask away if distance is too large or any particle is a dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

        const SoAFloatPrecision invdr2 = 1. / dr2;

        SoAFloatPrecision fract = sigmaSquared * invdr2;

        SoAFloatPrecision miem = 1;
        SoAFloatPrecision mien = 1;
        auto m_exp = _mexp;
        auto n_exp = _nexp;

        if constexpr (mode==0) {
          miem = m_exp & 1 ? sqrt(fract) : 1;
          mien = n_exp & 1 ? sqrt(fract) : 1;

          for (size_t k = 1; k < m_exp; k += 2) {
            miem = miem * fract;
          }
          for (size_t k = 1; k < n_exp; k += 2) {
            mien = mien * fract;
          }
        }
        else if constexpr (mode==1){
          miem = m_exp & 1 ? sqrt(fract) : 1;
          mien = n_exp & 1 ? sqrt(fract) : 1;

          m_exp>>=1;
          n_exp>>=1;

          while (m_exp||n_exp) {
            if (m_exp & 1) {
              miem = miem * fract;
            }
            if (n_exp & 1) {
              mien = mien * fract;
            }
            fract = fract * fract;
            m_exp >>= 1;
            n_exp >>= 1;
          }
        }
        else if constexpr (mode==2){
          auto chain = doubleAdditionChain;
          SoAFloatPrecision base =  sqrt(fract);
          SoAFloatPrecision a = 1.0, b = base;
          for(size_t k=0;k<chain_len;k++,chain>>=1){
            if(chain&1){
              SoAFloatPrecision tmp = b*a;
              a = b;
              b = tmp;
            }
            else{
              chain>>=1;
              k++;
              b =  b *b;
              if(chain&1){
                b = b *base;
              }
            }
          }
          miem = a;
          mien = b;
        }

        const SoAFloatPrecision mie = (cnepsilon * mien - cmepsilon * miem);
        const SoAFloatPrecision fac = mask ?  (invdr2 * mie) : 0;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        // newton 3
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;

        if (calculateGlobals) {
          const SoAFloatPrecision virialx = drx * fx;
          const SoAFloatPrecision virialy = dry * fy;
          const SoAFloatPrecision virialz = drz * fz;
          const SoAFloatPrecision potentialEnergy6 = mask ? (cepsilon * (miem-mien) + shift6) : 0.;

          // Add to the potential energy sum for each particle which is owned.
          // This results in obtaining 12 * the potential energy for the SoA.
          SoAFloatPrecision energyFactor = (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
                                           (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.);
          potentialEnergySum += potentialEnergy6 * energyFactor;

          virialSumX += virialx * energyFactor;
          virialSumY += virialy * energyFactor;
          virialSumZ += virialz * energyFactor;
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
    if (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // SoAFunctorSingle obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      if (newton3) {
        _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergySum * 0.5;
        _aosThreadData[threadnum].virialSumN3[0] += virialSumX * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += virialSumY * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += virialSumZ * 0.5;
      } else {
        _aosThreadData[threadnum].potentialEnergySumNoN3 += potentialEnergySum;
        _aosThreadData[threadnum].virialSumNoN3[0] += virialSumX;
        _aosThreadData[threadnum].virialSumNoN3[1] += virialSumY;
        _aosThreadData[threadnum].virialSumNoN3[2] += virialSumZ;
      }
    }
  }

  /**
   * @copydoc autopas::Functor::SoAFunctorPair()
   */
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
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
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
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle::AttributeNames::typeId>();

    // Checks whether the cells are halo cells.
    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmaSquared = _sigmaSquared;
    SoAFloatPrecision cepsilon = _cepsilon;
    SoAFloatPrecision cnepsilon = _cnepsilon;
    SoAFloatPrecision cmepsilon = _cmepsilon;


    // preload all sigma and epsilons for next vectorized region
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> cepsilons;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> cnepsilons;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> cmepsilons;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;
    if constexpr (useMixing) {
      sigmaSquareds.resize(soa2.getNumberOfParticles());
      cepsilons.resize(soa2.getNumberOfParticles());
      cnepsilons.resize(soa2.getNumberOfParticles());
      cmepsilons.resize(soa2.getNumberOfParticles());
      // if no mixing or mixing but no shift shift6 is constant therefore we do not need this vector.
      if constexpr (applyShift) {
        shift6s.resize(soa2.getNumberOfParticles());
      }
    }

    for (unsigned int i = 0; i < soa1.getNumberOfParticles(); ++i) {
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      // preload all sigma and epsilons for next vectorized region
      if constexpr (useMixing) {
        for (unsigned int j = 0; j < soa2.getNumberOfParticles(); ++j) {
          sigmaSquareds[j] = _PPLibrary->getMixingSigmaSquaredMie(typeptr1[i], typeptr2[j]);
          cepsilons[j] = _PPLibrary->getMixingCEpsilonMie(typeptr1[i], typeptr2[j]);
          cnepsilons[j] = _PPLibrary->getMixingCNEpsilonMie(typeptr1[i], typeptr2[j]);
          cmepsilons[j] = _PPLibrary->getMixingCMEpsilonMie(typeptr1[i], typeptr2[j]);
          if constexpr (applyShift) {
            shift6s[j] = _PPLibrary->getMixingShift6Mie(typeptr1[i], typeptr2[j]);
          }
        }
      }

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
      for (unsigned int j = 0; j < soa2.getNumberOfParticles(); ++j) {
        if constexpr (useMixing) {
          sigmaSquared = sigmaSquareds[j];
          cepsilon = cepsilons[j];
          cnepsilon = cnepsilons[j];
          cmepsilon = cmepsilons[j];
          if constexpr (applyShift) {
            shift6 = shift6s[j];
          }
        }

        const auto ownedStateJ = ownedStatePtr2[j];

        const SoAFloatPrecision drx = x1ptr[i] - x2ptr[j];
        const SoAFloatPrecision dry = y1ptr[i] - y2ptr[j];
        const SoAFloatPrecision drz = z1ptr[i] - z2ptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

        // Mask away if distance is too large or any particle is a dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

        const SoAFloatPrecision invdr2 = 1. / dr2;
        const SoAFloatPrecision fract = sigmaSquared * invdr2;

        SoAFloatPrecision miem = 1;
        SoAFloatPrecision mien = 1;
        auto m_exp = _mexp;
        auto n_exp = _nexp;

        if constexpr (mode==0) {
          miem = m_exp & 1 ? sqrt(fract) : 1;
          mien = n_exp & 1 ? sqrt(fract) : 1;

          for (size_t k = 1; k < m_exp; k += 2) {
            miem = miem * fract;
          }
          for (size_t k = 1; k < n_exp; k += 2) {
            mien = mien * fract;
          }
        }
        else if constexpr (mode==1){
          miem = m_exp & 1 ? sqrt(fract) : 1;
          mien = n_exp & 1 ? sqrt(fract) : 1;

          m_exp>>=1;
          n_exp>>=1;

          while (m_exp||n_exp) {
            if (m_exp & 1) {
              miem = miem * fract;
            }
            if (n_exp & 1) {
              mien = mien * fract;
            }
            fract = fract * fract;
            m_exp >>= 1;
            n_exp >>= 1;
          }
        }
        else if constexpr (mode==2){
          auto chain = doubleAdditionChain;
          SoAFloatPrecision base =  sqrt(fract);
          SoAFloatPrecision a = 1.0, b = base;
          for(size_t k=0;k<chain_len;k++,chain>>=1){
            if(chain&1){
              SoAFloatPrecision tmp = b*a;
              a = b;
              b = tmp;
            }
            else{
              chain>>=1;
              k++;
              b =  b *b;
              if(chain&1){
                b = b *base;
              }
            }
          }
          miem = a;
          mien = b;
        }

        const SoAFloatPrecision mie = (cnepsilon * mien - cmepsilon * miem);
        const SoAFloatPrecision fac = mask ?  (invdr2 * mie) : 0;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
        if (newton3) {
          fx2ptr[j] -= fx;
          fy2ptr[j] -= fy;
          fz2ptr[j] -= fz;
        }

        if constexpr (calculateGlobals) {
          SoAFloatPrecision virialx = drx * fx;
          SoAFloatPrecision virialy = dry * fy;
          SoAFloatPrecision virialz = drz * fz;
          SoAFloatPrecision potentialEnergy6 = (cepsilon * (miem-mien) + shift6) * mask;

          // Add to the potential energy sum for each particle which is owned.
          // This results in obtaining 12 * the potential energy for the SoA.
          const SoAFloatPrecision energyFactor =
              (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
              (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
          potentialEnergySum += potentialEnergy6 * energyFactor;
          virialSumX += virialx * energyFactor;
          virialSumY += virialy * energyFactor;
          virialSumZ += virialz * energyFactor;
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
    if (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // SoAFunctorPairImpl obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.
      if constexpr (newton3) {
        _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergySum * 0.5;
        _aosThreadData[threadnum].virialSumN3[0] += virialSumX * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += virialSumY * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += virialSumZ * 0.5;
      } else {
        _aosThreadData[threadnum].potentialEnergySumNoN3 += potentialEnergySum;
        _aosThreadData[threadnum].virialSumNoN3[0] += virialSumX;
        _aosThreadData[threadnum].virialSumNoN3[1] += virialSumY;
        _aosThreadData[threadnum].virialSumNoN3[2] += virialSumZ;
      }
    }
  }

 public:
  // clang-format off
  /**
   * @copydoc autopas::Functor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.getNumberOfParticles() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param epsilon24
   * @param sigmaSquared
   */
  void setParticleProperties(SoAFloatPrecision epsilon, SoAFloatPrecision sigmaSquared) {
    double c = static_cast<double>(_nexp) / static_cast<double>(_nexp - _mexp);
    c = c * pow(static_cast<double>(_nexp) / static_cast<double>(_mexp), static_cast<double>(_mexp) / (_nexp - _mexp));

    _cepsilon = c * epsilon;
    _cnepsilon = _nexp * _cepsilon;
    _cmepsilon = _mexp * _cepsilon;
    _sigmaSquared = sigmaSquared;
    if (applyShift) {
      _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShiftMie(epsilon, _sigmaSquared, _cutoffSquared, _nexp, _mexp);
    } else {
      _shift6 = 0.;
    }
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
   * @note molAType and molBType make no difference for MieFunctorFixed, but are kept to have a consistent interface for other
   * functors where they may.
   * @return the number of floating point operations
   */
   //TODO: correct this!
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
   * Get the potential Energy.
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
   * Get the virial.
   * @return
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

 private:
  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa.template begin<Particle::AttributeNames::typeId>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmaSquared = _sigmaSquared;
    SoAFloatPrecision cepsilon = _cepsilon;
    SoAFloatPrecision cnepsilon = _cnepsilon;
    SoAFloatPrecision cmepsilon = _cmepsilon;


    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    SoAFloatPrecision fxacc = 0;
    SoAFloatPrecision fyacc = 0;
    SoAFloatPrecision fzacc = 0;
    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    // checks whether particle i is owned.
    const auto ownedStateI = ownedStatePtr[indexFirst];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }

    // this is a magic number, that should correspond to at least
    // vectorization width*N have testet multiple sizes:
    // 4: does not give a speedup, slower than original AoSFunctor
    // 8: small speedup compared to AoS
    // 12: highest speedup compared to Aos
    // 16: smaller speedup
    // in theory this is a variable, we could auto-tune over...
#ifdef __AVX512F__
    // use a multiple of 8 for avx
    constexpr size_t vecsize = 16;
#else
    // for everything else 12 is faster
    constexpr size_t vecsize = 12;
#endif
    size_t joff = 0;

    // if the size of the verlet list is larger than the given size vecsize,
    // we will use a vectorized version.
    if (neighborListSize >= vecsize) {
      alignas(64) std::array<SoAFloatPrecision, vecsize> xtmp, ytmp, ztmp, xArr, yArr, zArr, fxArr, fyArr, fzArr;
      alignas(64) std::array<autopas::OwnershipState, vecsize> ownedStateArr{};
      // broadcast of the position of particle i
      for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
        xtmp[tmpj] = xptr[indexFirst];
        ytmp[tmpj] = yptr[indexFirst];
        ztmp[tmpj] = zptr[indexFirst];
      }
      // loop over the verlet list from 0 to x*vecsize
      for (; joff < neighborListSize - vecsize + 1; joff += vecsize) {
        // in each iteration we calculate the interactions of particle i with
        // vecsize particles in the neighborlist of particle i starting at
        // particle joff

        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> sigmaSquareds;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> cepsilons;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> cnepsilons;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> cmepsilons;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> shift6s;
        if constexpr (useMixing) {
          for (size_t j = 0; j < vecsize; j++) {
            sigmaSquareds[j] =
                _PPLibrary->getMixingSigmaSquaredMie(typeptr1[indexFirst], typeptr2[neighborListPtr[joff + j]]);
            cepsilons[j] = _PPLibrary->getMixingCEpsilonMie(typeptr1[indexFirst], typeptr2[neighborListPtr[joff + j]]);
            cnepsilons[j] = _PPLibrary->getMixingCNEpsilonMie(typeptr1[indexFirst], typeptr2[neighborListPtr[joff + j]]);
            cmepsilons[j] = _PPLibrary->getMixingCMEpsilonMie(typeptr1[indexFirst], typeptr2[neighborListPtr[joff + j]]);
            if constexpr (applyShift) {
              shift6s[j] = _PPLibrary->getMixingShift6Mie(typeptr1[indexFirst], typeptr2[neighborListPtr[joff + j]]);
            }
          }
        }

        // gather position of particle j
#pragma omp simd safelen(vecsize)
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          xArr[tmpj] = xptr[neighborListPtr[joff + tmpj]];
          yArr[tmpj] = yptr[neighborListPtr[joff + tmpj]];
          zArr[tmpj] = zptr[neighborListPtr[joff + tmpj]];
          ownedStateArr[tmpj] = ownedStatePtr[neighborListPtr[joff + tmpj]];
        }
        // do omp simd with reduction of the interaction
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ) safelen(vecsize)
        for (size_t j = 0; j < vecsize; j++) {
          if constexpr (useMixing) {
            sigmaSquared = sigmaSquareds[j];
            cepsilon = cepsilons[j];
            cnepsilon = cepsilons[j];
            cmepsilon = cepsilons[j];
            if constexpr (applyShift) {
              shift6 = shift6s[j];
            }
          }
          // const size_t j = currentList[jNeighIndex];

          const auto ownedStateJ = ownedStateArr[j];

          const SoAFloatPrecision drx = xtmp[j] - xArr[j];
          const SoAFloatPrecision dry = ytmp[j] - yArr[j];
          const SoAFloatPrecision drz = ztmp[j] - zArr[j];

          const SoAFloatPrecision drx2 = drx * drx;
          const SoAFloatPrecision dry2 = dry * dry;
          const SoAFloatPrecision drz2 = drz * drz;

          const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

          // Mask away if distance is too large or any particle is a dummy.
          // Particle ownedStateI was already checked previously.
          const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

          const SoAFloatPrecision invdr2 = 1. / dr2;
          const SoAFloatPrecision fract = sigmaSquared * invdr2;

          SoAFloatPrecision miem = 1;
          SoAFloatPrecision mien = 1;
          auto m_exp = _mexp;
          auto n_exp = _nexp;

          if constexpr (mode==0) {
            miem = m_exp & 1 ? sqrt(fract) : 1;
            mien = n_exp & 1 ? sqrt(fract) : 1;

            for (size_t k = 1; k < m_exp; k += 2) {
              miem = miem * fract;
            }
            for (size_t k = 1; k < n_exp; k += 2) {
              mien = mien * fract;
            }
          }
          else if constexpr (mode==1){
            miem = m_exp & 1 ? sqrt(fract) : 1;
            mien = n_exp & 1 ? sqrt(fract) : 1;

            m_exp>>=1;
            n_exp>>=1;

            while (m_exp||n_exp) {
              if (m_exp & 1) {
                miem = miem * fract;
              }
              if (n_exp & 1) {
                mien = mien * fract;
              }
              fract = fract * fract;
              m_exp >>= 1;
              n_exp >>= 1;
            }
          }
          else if constexpr (mode==2){
            auto chain = doubleAdditionChain;
            SoAFloatPrecision base =  sqrt(fract);
            SoAFloatPrecision a = 1.0, b = base;
            for(size_t k=0;k<chain_len;k++,chain>>=1){
              if(chain&1){
                SoAFloatPrecision tmp = b*a;
                a = b;
                b = tmp;
              }
              else{
                chain>>=1;
                k++;
                b =  b *b;
                if(chain&1){
                  b = b *base;
                }
              }
            }
            miem = a;
            mien = b;
          }

          const SoAFloatPrecision mie = (cnepsilon * mien - cmepsilon * miem);
          const SoAFloatPrecision fac = mask ?  (invdr2 * mie) : 0;

          const SoAFloatPrecision fx = drx * fac;
          const SoAFloatPrecision fy = dry * fac;
          const SoAFloatPrecision fz = drz * fac;

          fxacc += fx;
          fyacc += fy;
          fzacc += fz;
          if (newton3) {
            fxArr[j] = fx;
            fyArr[j] = fy;
            fzArr[j] = fz;
          }
          if (calculateGlobals) {
            SoAFloatPrecision virialx = drx * fx;
            SoAFloatPrecision virialy = dry * fy;
            SoAFloatPrecision virialz = drz * fz;
            SoAFloatPrecision potentialEnergy6 = mask ? (cepsilon * (miem - mien) + shift6) : 0.;

            const SoAFloatPrecision energyFactor =
                (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
                (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
            potentialEnergySum += potentialEnergy6 * energyFactor;
            virialSumX += virialx * energyFactor;
            virialSumY += virialy * energyFactor;
            virialSumZ += virialz * energyFactor;
          }
        }
        // scatter the forces to where they belong, this is only needed for newton3
        if (newton3) {
#pragma omp simd safelen(vecsize)
          for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
            const size_t j = neighborListPtr[joff + tmpj];
            fxptr[j] -= fxArr[tmpj];
            fyptr[j] -= fyArr[tmpj];
            fzptr[j] -= fzArr[tmpj];
          }
        }
      }
    }
    // this loop goes over the remainder and uses no optimizations
    for (size_t jNeighIndex = joff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;
      if constexpr (useMixing) {
        sigmaSquared = _PPLibrary->getMixingSigmaSquaredMie(typeptr1[indexFirst], typeptr2[j]);
        cepsilon = _PPLibrary->getMixingCEpsilonMie(typeptr1[indexFirst], typeptr2[j]);
        cnepsilon = _PPLibrary->getMixingCEpsilonMie(typeptr1[indexFirst], typeptr2[j]);
        cmepsilon = _PPLibrary->getMixingCEpsilonMie(typeptr1[indexFirst], typeptr2[j]);
        if constexpr (applyShift) {
          shift6 = _PPLibrary->getMixingShift6Mie(typeptr1[indexFirst], typeptr2[j]);
        }
      }

      const auto ownedStateJ = ownedStatePtr[j];
      if (ownedStateJ == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision drx = xptr[indexFirst] - xptr[j];
      const SoAFloatPrecision dry = yptr[indexFirst] - yptr[j];
      const SoAFloatPrecision drz = zptr[indexFirst] - zptr[j];

      const SoAFloatPrecision drx2 = drx * drx;
      const SoAFloatPrecision dry2 = dry * dry;
      const SoAFloatPrecision drz2 = drz * drz;

      const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

      if (dr2 > cutoffSquared) {
        continue;
      }
      const SoAFloatPrecision invdr2 = 1. / dr2;
      const SoAFloatPrecision fract = sigmaSquared * invdr2;

      SoAFloatPrecision miem = 1;
      SoAFloatPrecision mien = 1;
      auto m_exp = _mexp;
      auto n_exp = _nexp;

      if constexpr (mode==0) {
        miem = m_exp & 1 ? sqrt(fract) : 1;
        mien = n_exp & 1 ? sqrt(fract) : 1;

        for (size_t k = 1; k < m_exp; k += 2) {
          miem = miem * fract;
        }
        for (size_t k = 1; k < n_exp; k += 2) {
          mien = mien * fract;
        }
      }
      else if constexpr (mode==1){
        miem = m_exp & 1 ? sqrt(fract) : 1;
        mien = n_exp & 1 ? sqrt(fract) : 1;

        m_exp>>=1;
        n_exp>>=1;

        while (m_exp||n_exp) {
          if (m_exp & 1) {
            miem = miem * fract;
          }
          if (n_exp & 1) {
            mien = mien * fract;
          }
          fract = fract * fract;
          m_exp >>= 1;
          n_exp >>= 1;
        }
      }
      else if constexpr (mode==2){
        auto chain = doubleAdditionChain;
        SoAFloatPrecision base =  sqrt(fract);
        SoAFloatPrecision a = 1.0, b = base;
        for(size_t k=0;k<chain_len;k++,chain>>=1){
          if(chain&1){
            SoAFloatPrecision tmp = b*a;
            a = b;
            b = tmp;
          }
          else{
            chain>>=1;
            k++;
            b =  b *b;
            if(chain&1){
              b = b *base;
            }
          }
        }
        miem = a;
        mien = b;
      }

      const SoAFloatPrecision mie = (cnepsilon * mien - cmepsilon * miem);
      const SoAFloatPrecision fac = (invdr2 * mie);

      const SoAFloatPrecision fx = drx * fac;
      const SoAFloatPrecision fy = dry * fac;
      const SoAFloatPrecision fz = drz * fac;

      fxacc += fx;
      fyacc += fy;
      fzacc += fz;
      if (newton3) {
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;
      }
      if (calculateGlobals) {
        SoAFloatPrecision virialx = drx * fx;
        SoAFloatPrecision virialy = dry * fy;
        SoAFloatPrecision virialz = drz * fz;
        SoAFloatPrecision potentialEnergy6 = (cepsilon * (miem - mien) + shift6);

        // Add to the potential energy sum for each particle which is owned.
        // This results in obtaining 12 * the potential energy for the SoA.
        const SoAFloatPrecision energyFactor =
            (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
            (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
        potentialEnergySum += potentialEnergy6 * energyFactor;
        virialSumX += virialx * energyFactor;
        virialSumY += virialy * energyFactor;
        virialSumZ += virialz * energyFactor;
      }
    }

    if (fxacc != 0 or fyacc != 0 or fzacc != 0) {
      fxptr[indexFirst] += fxacc;
      fyptr[indexFirst] += fyacc;
      fzptr[indexFirst] += fzacc;
    }

    if (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // SoAFunctorVerlet obtains the potential energy * 12. For non-newton3, this sum is divided by 12 in
      // post-processing. For newton3, this sum is only divided by 6 in post-processing, so must be divided by 2 here.

      if (newton3) {
        _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergySum * 0.5;
        _aosThreadData[threadnum].virialSumN3[0] += virialSumX * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += virialSumY * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += virialSumZ * 0.5;
      } else {
        _aosThreadData[threadnum].potentialEnergySumNoN3 += potentialEnergySum;
        _aosThreadData[threadnum].virialSumNoN3[0] += virialSumX;
        _aosThreadData[threadnum].virialSumNoN3[1] += virialSumY;
        _aosThreadData[threadnum].virialSumNoN3[2] += virialSumZ;
      }
    }
  }

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

  const double _cutoffSquared;
  static constexpr uint16_t _nexp = 12;
  static constexpr uint16_t _mexp = 6;
  int doubleAdditionChain;
  uint16_t chain_len;
  // not const because they might be reset through PPL
  double _cepsilon, _cnepsilon, _cmepsilon, _sigmaSquared, _shift6 = 0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

  static constexpr int mode = 0;
  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;
};
}  // namespace mdLib

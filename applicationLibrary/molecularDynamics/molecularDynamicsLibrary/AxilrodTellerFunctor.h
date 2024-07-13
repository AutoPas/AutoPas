/**
 * @file AxilrodTellerFunctor.h
 * @author M. Muehlhaeusser
 * @date 25/07/23
 */

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
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
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool relevantForTuning = true>
class AxilrodTellerFunctor
    : public autopas::TriwiseFunctor<
          Particle, AxilrodTellerFunctor<Particle, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
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
  AxilrodTellerFunctor() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit AxilrodTellerFunctor(double cutoff, void * /*dummy*/)
      : autopas::TriwiseFunctor<
            Particle, AxilrodTellerFunctor<Particle, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
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
  explicit AxilrodTellerFunctor(double cutoff) : AxilrodTellerFunctor(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like nu.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit AxilrodTellerFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : AxilrodTellerFunctor(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  /**
   * Returns name of functor. Intended for use with the iteration logger, to differentiate between calls to
   * iterateTriwise using different functors in the logs.
   * @return name of functor.
   */
  virtual std::string getName() { return "AxilrodTellerFunctorAutoVec"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle &i, Particle &j, Particle &k, bool newton3) {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy() or k.isDummy()) {
      return;
    }
    auto nu = _nu;
    if constexpr (useMixing) {
      nu = _PPLibrary->getMixingNu(i.getTypeId(), j.getTypeId(), k.getTypeId());
    }
    auto drij = j.getR() - i.getR();
    auto drjk = k.getR() - j.getR();
    auto drki = i.getR() - k.getR();

    double dr2ij = autopas::utils::ArrayMath::dot(drij, drij);
    double dr2jk = autopas::utils::ArrayMath::dot(drjk, drjk);
    double dr2ki = autopas::utils::ArrayMath::dot(drki, drki);

    // Check cutoff
    if (dr2ij > _cutoffSquared or dr2jk > _cutoffSquared or dr2ki > _cutoffSquared) {
      return;
    }

    // Dot products of distances belonging to one particle
    double dr2i = autopas::utils::ArrayMath::dot(drij, drki);
    double dr2j = autopas::utils::ArrayMath::dot(drij, drjk);
    double dr2k = autopas::utils::ArrayMath::dot(drjk, drki);

    double dr2ijk = dr2i * dr2j * dr2k;

    double dr2 = dr2ij * dr2jk * dr2ki;
    double dr5 = dr2 * dr2 * std::sqrt(dr2);
    double invdr5 = nu / dr5;

    auto fi = drjk * dr2i * (dr2j - dr2k) + drij * (dr2j * dr2k - dr2jk * dr2ki + 5.0 * dr2ijk / dr2ij) +
              drki * (-dr2j * dr2k + dr2ij * dr2jk - 5.0 * dr2ijk / dr2ki);
    fi *= 3.0 * invdr5;
    i.addF(fi);

    auto fj = fi;
    auto fk = fi;
    if (newton3) {
      fj = drki * dr2j * (dr2k - dr2i) + drij * (-dr2i * dr2k + dr2jk * dr2ki - 5.0 * dr2ijk / dr2ij) +
           drjk * (dr2i * dr2k - dr2ij * dr2ki + 5.0 * dr2ijk / dr2jk);
      fj *= 3.0 * invdr5;
      j.addF(fj);

      /* auto fk = drij * dr2k * (dr2i - dr2j)
                + drjk * (- dr2i * dr2j + dr2ij * dr2ki - 5.0 * dr2ijk / dr2jk)
                + drki * (dr2i * dr2j - dr2ij * dr2jk + 5.0 * dr2ijk / dr2ki);
      fk *= 3.0 * invdr5; */
      fk = (fi + fj) * (-1.0);
      k.addF(fk);
    }

    if (calculateGlobals) {
      // Virial is calculated as f_i * r_i
      auto virialI = fi * i.getR();
      // Calculate third of total potential energy from 3-body interaction
      double potentialEnergy = invdr5 * (dr2 - 3.0 * dr2ijk) / 3.0;

      const int threadnum = autopas::autopas_get_thread_num();
      if (i.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialI;
      }
      // for non-newton3 particles j and/or k will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        auto virialJ = fj * j.getR();
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialJ;
      }
      if (newton3 and k.isOwned()) {
        auto virialK = fk * k.getR();
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialK;
      }
    }
  }

  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImpl<true>(soa);
    } else {
      SoAFunctorSingleImpl<false>(soa);
    }
  }

  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

  /**
   * Functor for structure of arrays (SoA)
   *
   * This functor calculates the forces
   * between all particles of soa1 and soa2 and soa3.
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
   * @param soa3 Third structure of arrays.
   * @param newton3 defines whether or whether not to use newton 3
   */
  void SoAFunctorTriple(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                        autopas::SoAView<SoAArraysType> soa3, const bool newton3) {
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorTriple() is not implemented.");
    if (newton3) {
      SoAFunctorTripleImpl<true>(soa1, soa2, soa3);
    } else {
      SoAFunctorTripleImpl<false>(soa1, soa2, soa3);
    }
  }

 private:
  /**
   * Implementation function of SoAFunctorSingle(soa, newton3)
   *
   * @tparam newton3
   * @param soa
   */
  template <bool newton3>
  void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorSingle() is not implemented.");
    //  TODO
    if (soa.size() == 0) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    for (size_t i = soa.size() - 1; static_cast<long>(i) >= 2; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        continue;
      }
      SoAFloatPrecision fxiacc = 0.;
      SoAFloatPrecision fyiacc = 0.;
      SoAFloatPrecision fziacc = 0.;

      const SoAFloatPrecision xi = xptr[i];
      const SoAFloatPrecision yi = yptr[i];
      const SoAFloatPrecision zi = zptr[i];

      for (size_t j = i - 1; static_cast<long>(j) >= 1; --j) {
        if (ownedStatePtr[j] == autopas::OwnershipState::dummy) {
          continue;
        }
        SoAFloatPrecision fxjacc = 0.;
        SoAFloatPrecision fyjacc = 0.;
        SoAFloatPrecision fzjacc = 0.;

        const SoAFloatPrecision xj = xptr[j];
        const SoAFloatPrecision yj = yptr[j];
        const SoAFloatPrecision zj = zptr[j];

        const SoAFloatPrecision drxij = xj - xi;
        const SoAFloatPrecision dryij = yj - yi;
        const SoAFloatPrecision drzij = zj - zi;

        const SoAFloatPrecision drxij2 = drxij * drxij;
        const SoAFloatPrecision dryij2 = dryij * dryij;
        const SoAFloatPrecision drzij2 = drzij * drzij;

        const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
        if (drij2 > _cutoffSquared) {
          continue;
        }

        for (size_t k = 0; k < j; ++k) {
          if (ownedStatePtr[k] == autopas::OwnershipState::dummy) {
            continue;
          }
          const SoAFloatPrecision xk = xptr[k];
          const SoAFloatPrecision yk = yptr[k];
          const SoAFloatPrecision zk = zptr[k];

          const SoAFloatPrecision drxjk = xk - xj;
          const SoAFloatPrecision dryjk = yk - yj;
          const SoAFloatPrecision drzjk = zk - zj;

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          const SoAFloatPrecision drxki = xi - xk;
          const SoAFloatPrecision dryki = yi - yk;
          const SoAFloatPrecision drzki = zi - zk;

          const SoAFloatPrecision drxki2 = drxki * drxki;
          const SoAFloatPrecision dryki2 = dryki * dryki;
          const SoAFloatPrecision drzki2 = drzki * drzki;

          const SoAFloatPrecision drki2 = drxki2 + dryki2 + drzki2;

          if (drjk2 > _cutoffSquared or drki2 > _cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision drxi2 = drxij * drxki;
          const SoAFloatPrecision dryi2 = dryij * dryki;
          const SoAFloatPrecision drzi2 = drzij * drzki;

          const SoAFloatPrecision dri2 = drxi2 + dryi2 + drzi2;

          const SoAFloatPrecision drxj2 = drxij * drxjk;
          const SoAFloatPrecision dryj2 = dryij * dryjk;
          const SoAFloatPrecision drzj2 = drzij * drzjk;

          const SoAFloatPrecision drj2 = drxj2 + dryj2 + drzj2;

          const SoAFloatPrecision drxk2 = drxjk * drxki;
          const SoAFloatPrecision dryk2 = dryjk * dryki;
          const SoAFloatPrecision drzk2 = drzjk * drzki;

          const SoAFloatPrecision drk2 = drxk2 + dryk2 + drzk2;

          const SoAFloatPrecision drijk2 = dri2 * drj2 * drk2;

          const SoAFloatPrecision dr2 = drij2 * drjk2 * drki2;
          const SoAFloatPrecision dr5 = dr2 * dr2 * std::sqrt(dr2);

          auto nu = _nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(typeptr[i], typeptr[j], typeptr[k]);
          }
          const SoAFloatPrecision invdr5 = nu / dr5;
          const SoAFloatPrecision invdr53 = invdr5 * 3.0;

          const SoAFloatPrecision fxi =
              (drxjk * dri2 * (drj2 - drk2) + drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               drxki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          const SoAFloatPrecision fyi =
              (dryjk * dri2 * (drj2 - drk2) + dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               dryki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          const SoAFloatPrecision fzi =
              (drzjk * dri2 * (drj2 - drk2) + drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               drzki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          fxiacc += fxi;
          fyiacc += fyi;
          fziacc += fzi;

          const SoAFloatPrecision fxj =
              (drxki * drj2 * (drk2 - dri2) + drxij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
               drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
              invdr53;
          const SoAFloatPrecision fyj =
              (dryki * drj2 * (drk2 - dri2) + dryij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
               dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
              invdr53;
          const SoAFloatPrecision fzj =
              (drzki * drj2 * (drk2 - dri2) + drzij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
               drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
              invdr53;
          fxjacc += fxj;
          fyjacc += fyj;
          fzjacc += fzj;

          const SoAFloatPrecision nfxk = fxi + fxj;
          const SoAFloatPrecision nfyk = fyi + fyj;
          const SoAFloatPrecision nfzk = fzi + fzj;
          fxptr[k] -= nfxk;
          fyptr[k] -= nfyk;
          fzptr[k] -= nfzk;

          if constexpr (calculateGlobals) {
            const double potentialEnergy = invdr5 * (dr2 - 3.0 * drijk2);
            potentialEnergySum += potentialEnergy;
            virialSumX += fxi * xi + fxj * xj - nfxk * xk;
            virialSumY += fyi * yi + fyj * yj - nfyk * yk;
            virialSumZ += fzi * zi + fzj * zj - nfzk * zk;
          }
        }
        fxptr[j] += fxjacc;
        fyptr[j] += fyjacc;
        fzptr[j] += fzjacc;
      }
      fxptr[i] += fxiacc;
      fyptr[i] += fyiacc;
      fzptr[i] += fziacc;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();
      _aosThreadData[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadData[threadnum].virialSum[0] += virialSumX;
      _aosThreadData[threadnum].virialSum[1] += virialSumY;
      _aosThreadData[threadnum].virialSum[2] += virialSumZ;
    }
  }

  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    // TODO: should always calculate forces for all particles in soa1, even when newton3 == false
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorPair() is not implemented.");
    // TODO
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      SoAFloatPrecision fxiacc = 0;
      SoAFloatPrecision fyiacc = 0;
      SoAFloatPrecision fziacc = 0;

      const SoAFloatPrecision xi = x1ptr[i];
      const SoAFloatPrecision yi = y1ptr[i];
      const SoAFloatPrecision zi = z1ptr[i];

      // particle 2 from soa1 and 3 from soa2

      for (size_t j = i + 1; j < soa1.size(); ++j) {
        if (ownedState1ptr[j] == autopas::OwnershipState::dummy) {
          continue;
        }

        SoAFloatPrecision fxjacc = 0.;
        SoAFloatPrecision fyjacc = 0.;
        SoAFloatPrecision fzjacc = 0.;

        const SoAFloatPrecision xj = x1ptr[j];
        const SoAFloatPrecision yj = y1ptr[j];
        const SoAFloatPrecision zj = z1ptr[j];

        const SoAFloatPrecision drxij = xj - xi;
        const SoAFloatPrecision dryij = yj - yi;
        const SoAFloatPrecision drzij = zj - zi;

        const SoAFloatPrecision drxij2 = drxij * drxij;
        const SoAFloatPrecision dryij2 = dryij * dryij;
        const SoAFloatPrecision drzij2 = drzij * drzij;

        const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
        if (drij2 > _cutoffSquared) {
          continue;
        }

        for (size_t k = 0; k < soa2.size(); ++k) {
          if (ownedState2ptr[k] == autopas::OwnershipState::dummy) {
            continue;
          }
          const SoAFloatPrecision xk = x2ptr[k];
          const SoAFloatPrecision yk = y2ptr[k];
          const SoAFloatPrecision zk = z2ptr[k];

          const SoAFloatPrecision drxjk = xk - xj;
          const SoAFloatPrecision dryjk = yk - yj;
          const SoAFloatPrecision drzjk = zk - zj;

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          const SoAFloatPrecision drxki = xi - xk;
          const SoAFloatPrecision dryki = yi - yk;
          const SoAFloatPrecision drzki = zi - zk;

          const SoAFloatPrecision drxki2 = drxki * drxki;
          const SoAFloatPrecision dryki2 = dryki * dryki;
          const SoAFloatPrecision drzki2 = drzki * drzki;

          const SoAFloatPrecision drki2 = drxki2 + dryki2 + drzki2;

          if (drjk2 > _cutoffSquared or drki2 > _cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision drxi2 = drxij * drxki;
          const SoAFloatPrecision dryi2 = dryij * dryki;
          const SoAFloatPrecision drzi2 = drzij * drzki;

          const SoAFloatPrecision dri2 = drxi2 + dryi2 + drzi2;

          const SoAFloatPrecision drxj2 = drxij * drxjk;
          const SoAFloatPrecision dryj2 = dryij * dryjk;
          const SoAFloatPrecision drzj2 = drzij * drzjk;

          const SoAFloatPrecision drj2 = drxj2 + dryj2 + drzj2;

          const SoAFloatPrecision drxk2 = drxjk * drxki;
          const SoAFloatPrecision dryk2 = dryjk * dryki;
          const SoAFloatPrecision drzk2 = drzjk * drzki;

          const SoAFloatPrecision drk2 = drxk2 + dryk2 + drzk2;

          const SoAFloatPrecision drijk2 = dri2 * drj2 * drk2;

          const SoAFloatPrecision dr2 = drij2 * drjk2 * drki2;
          const SoAFloatPrecision dr5 = dr2 * dr2 * std::sqrt(dr2);

          auto nu = _nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(type1ptr[i], type1ptr[j], type2ptr[k]);
          }
          const SoAFloatPrecision invdr5 = nu / dr5;
          const SoAFloatPrecision invdr53 = invdr5 * 3.0;

          const SoAFloatPrecision fxi =
              (drxjk * dri2 * (drj2 - drk2) + drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               drxki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          const SoAFloatPrecision fyi =
              (dryjk * dri2 * (drj2 - drk2) + dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               dryki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          const SoAFloatPrecision fzi =
              (drzjk * dri2 * (drj2 - drk2) + drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               drzki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          fxiacc += fxi;
          fyiacc += fyi;
          fziacc += fzi;

          const SoAFloatPrecision fxj =
              (drxki * drj2 * (drk2 - dri2) + drxij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
               drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
              invdr53;
          const SoAFloatPrecision fyj =
              (dryki * drj2 * (drk2 - dri2) + dryij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
               dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
              invdr53;
          const SoAFloatPrecision fzj =
              (drzki * drj2 * (drk2 - dri2) + drzij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
               drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
              invdr53;

          fxjacc += fxj;
          fyjacc += fyj;
          fzjacc += fzj;

          if (newton3) {
            const SoAFloatPrecision nfxk = fxi + fxj;
            const SoAFloatPrecision nfyk = fyi + fyj;
            const SoAFloatPrecision nfzk = fzi + fzj;

            fx2ptr[k] -= nfxk;
            fy2ptr[k] -= nfyk;
            fz2ptr[k] -= nfzk;

            if constexpr (calculateGlobals) {
              virialSumX -= nfxk * xk;
              virialSumY -= nfyk * yk;
              virialSumZ -= nfzk * zk;
            }
          }

          if constexpr (calculateGlobals) {
            double potentialEnergy = (newton3 ? invdr53 : 2.0 * invdr5) * ((1.0 / 3.0) * dr2 - drijk2);
            potentialEnergySum += potentialEnergy;
            virialSumX += fxi * xi + fxj * xj;
            virialSumY += fyi * yi + fyj * yj;
            virialSumZ += fzi * zi + fzj * zj;
          }
        }
        fx1ptr[j] += fxjacc;
        fy1ptr[j] += fyjacc;
        fz1ptr[j] += fzjacc;
      }

      // both particles 2 and 3 from soa2

      for (size_t j = soa2.size() - 1; static_cast<long>(j) >= 1; --j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy) {
          continue;
        }

        SoAFloatPrecision fxjacc = 0.;
        SoAFloatPrecision fyjacc = 0.;
        SoAFloatPrecision fzjacc = 0.;

        const SoAFloatPrecision xj = x2ptr[j];
        const SoAFloatPrecision yj = y2ptr[j];
        const SoAFloatPrecision zj = z2ptr[j];

        const SoAFloatPrecision drxij = xj - xi;
        const SoAFloatPrecision dryij = yj - yi;
        const SoAFloatPrecision drzij = zj - zi;

        const SoAFloatPrecision drxij2 = drxij * drxij;
        const SoAFloatPrecision dryij2 = dryij * dryij;
        const SoAFloatPrecision drzij2 = drzij * drzij;

        const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
        if (drij2 > _cutoffSquared) {
          continue;
        }
        // particle 3 from soa 2
        for (size_t k = 0; k < j; ++k) {
          if (ownedState2ptr[k] == autopas::OwnershipState::dummy) {
            continue;
          }
          const SoAFloatPrecision xk = x2ptr[k];
          const SoAFloatPrecision yk = y2ptr[k];
          const SoAFloatPrecision zk = z2ptr[k];

          const SoAFloatPrecision drxjk = xk - xj;
          const SoAFloatPrecision dryjk = yk - yj;
          const SoAFloatPrecision drzjk = zk - zj;

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          const SoAFloatPrecision drxki = xi - xk;
          const SoAFloatPrecision dryki = yi - yk;
          const SoAFloatPrecision drzki = zi - zk;

          const SoAFloatPrecision drxki2 = drxki * drxki;
          const SoAFloatPrecision dryki2 = dryki * dryki;
          const SoAFloatPrecision drzki2 = drzki * drzki;

          const SoAFloatPrecision drki2 = drxki2 + dryki2 + drzki2;

          if (drjk2 > _cutoffSquared or drki2 > _cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision drxi2 = drxij * drxki;
          const SoAFloatPrecision dryi2 = dryij * dryki;
          const SoAFloatPrecision drzi2 = drzij * drzki;

          const SoAFloatPrecision dri2 = drxi2 + dryi2 + drzi2;

          const SoAFloatPrecision drxj2 = drxij * drxjk;
          const SoAFloatPrecision dryj2 = dryij * dryjk;
          const SoAFloatPrecision drzj2 = drzij * drzjk;

          const SoAFloatPrecision drj2 = drxj2 + dryj2 + drzj2;

          const SoAFloatPrecision drxk2 = drxjk * drxki;
          const SoAFloatPrecision dryk2 = dryjk * dryki;
          const SoAFloatPrecision drzk2 = drzjk * drzki;

          const SoAFloatPrecision drk2 = drxk2 + dryk2 + drzk2;

          const SoAFloatPrecision drijk2 = dri2 * drj2 * drk2;

          const SoAFloatPrecision dr2 = drij2 * drjk2 * drki2;
          const SoAFloatPrecision dr5 = dr2 * dr2 * std::sqrt(dr2);

          auto nu = _nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(type1ptr[i], type2ptr[j], type2ptr[k]);
          }
          const SoAFloatPrecision invdr5 = nu / dr5;
          const SoAFloatPrecision invdr53 = invdr5 * 3.0;

          const SoAFloatPrecision fxi =
              (drxjk * dri2 * (drj2 - drk2) + drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               drxki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          const SoAFloatPrecision fyi =
              (dryjk * dri2 * (drj2 - drk2) + dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               dryki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          const SoAFloatPrecision fzi =
              (drzjk * dri2 * (drj2 - drk2) + drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               drzki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          fxiacc += fxi;
          fyiacc += fyi;
          fziacc += fzi;

          if (newton3) {
            const SoAFloatPrecision fxj =
                (drxki * drj2 * (drk2 - dri2) + drxij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                 drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
                invdr53;
            const SoAFloatPrecision fyj =
                (dryki * drj2 * (drk2 - dri2) + dryij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                 dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
                invdr53;
            const SoAFloatPrecision fzj =
                (drzki * drj2 * (drk2 - dri2) + drzij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                 drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
                invdr53;
            fxjacc += fxj;
            fyjacc += fyj;
            fzjacc += fzj;

            const SoAFloatPrecision nfxk = fxi + fxj;
            const SoAFloatPrecision nfyk = fyi + fyj;
            const SoAFloatPrecision nfzk = fzi + fzj;
            fx2ptr[k] -= nfxk;
            fy2ptr[k] -= nfyk;
            fz2ptr[k] -= nfzk;

            if constexpr (calculateGlobals) {
              virialSumX += fxj * xj - nfxk * xk;
              virialSumY += fyj * yj - nfyk * yk;
              virialSumZ += fzj * zj - nfzk * zk;
            }
          }

          if constexpr (calculateGlobals) {
            double potentialEnergy = (newton3 ? invdr53 : invdr5) * ((1.0 / 3.0) * dr2 - drijk2);
            potentialEnergySum += potentialEnergy;
            virialSumX += fxi * xi;
            virialSumY += fyi * yi;
            virialSumZ += fzi * zi;
          }
        }
        if constexpr (newton3) {
          fx2ptr[j] += fxjacc;
          fy2ptr[j] += fyjacc;
          fz2ptr[j] += fzjacc;
        }
      }

      fx1ptr[i] += fxiacc;
      fy1ptr[i] += fyiacc;
      fz1ptr[i] += fziacc;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();
      _aosThreadData[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadData[threadnum].virialSum[0] += virialSumX;
      _aosThreadData[threadnum].virialSum[1] += virialSumY;
      _aosThreadData[threadnum].virialSum[2] += virialSumZ;
    }
  }

  /**
   * Implementation function of SoAFunctorTriple(soa1, soa2, soa3, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   * @param soa3
   */
  template <bool newton3>
  void SoAFunctorTripleImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                   autopas::SoAView<SoAArraysType> soa3) {
    // autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorTriple() is not implemented.");
    //   TODO
    if (soa1.size() == 0 or soa2.size() == 0 or soa3.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x3ptr = soa3.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y3ptr = soa3.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z3ptr = soa3.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState3ptr = soa3.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx3ptr = soa3.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy3ptr = soa3.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz3ptr = soa3.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type3ptr = soa3.template begin<Particle::AttributeNames::typeId>();

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        continue;
      }
      SoAFloatPrecision fxiacc = 0.;
      SoAFloatPrecision fyiacc = 0.;
      SoAFloatPrecision fziacc = 0.;

      const SoAFloatPrecision xi = x1ptr[i];
      const SoAFloatPrecision yi = y1ptr[i];
      const SoAFloatPrecision zi = z1ptr[i];

      for (size_t j = 0; j < soa2.size(); ++j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy) {
          continue;
        }
        SoAFloatPrecision fxjacc = 0.;
        SoAFloatPrecision fyjacc = 0.;
        SoAFloatPrecision fzjacc = 0.;

        const SoAFloatPrecision xj = x2ptr[j];
        const SoAFloatPrecision yj = y2ptr[j];
        const SoAFloatPrecision zj = z2ptr[j];

        const SoAFloatPrecision drxij = xj - xi;
        const SoAFloatPrecision dryij = yj - yi;
        const SoAFloatPrecision drzij = zj - zi;

        const SoAFloatPrecision drxij2 = drxij * drxij;
        const SoAFloatPrecision dryij2 = dryij * dryij;
        const SoAFloatPrecision drzij2 = drzij * drzij;

        const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
        if (drij2 > _cutoffSquared) {
          continue;
        }

        for (size_t k = 0; k < soa3.size(); ++k) {
          if (ownedState3ptr[k] == autopas::OwnershipState::dummy) {
            continue;
          }
          const SoAFloatPrecision xk = x3ptr[k];
          const SoAFloatPrecision yk = y3ptr[k];
          const SoAFloatPrecision zk = z3ptr[k];

          const SoAFloatPrecision drxjk = xk - xj;
          const SoAFloatPrecision dryjk = yk - yj;
          const SoAFloatPrecision drzjk = zk - zj;

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          const SoAFloatPrecision drxki = xi - xk;
          const SoAFloatPrecision dryki = yi - yk;
          const SoAFloatPrecision drzki = zi - zk;

          const SoAFloatPrecision drxki2 = drxki * drxki;
          const SoAFloatPrecision dryki2 = dryki * dryki;
          const SoAFloatPrecision drzki2 = drzki * drzki;

          const SoAFloatPrecision drki2 = drxki2 + dryki2 + drzki2;

          if (drjk2 > _cutoffSquared or drki2 > _cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision drxi2 = drxij * drxki;
          const SoAFloatPrecision dryi2 = dryij * dryki;
          const SoAFloatPrecision drzi2 = drzij * drzki;

          const SoAFloatPrecision dri2 = drxi2 + dryi2 + drzi2;

          const SoAFloatPrecision drxj2 = drxij * drxjk;
          const SoAFloatPrecision dryj2 = dryij * dryjk;
          const SoAFloatPrecision drzj2 = drzij * drzjk;

          const SoAFloatPrecision drj2 = drxj2 + dryj2 + drzj2;

          const SoAFloatPrecision drxk2 = drxjk * drxki;
          const SoAFloatPrecision dryk2 = dryjk * dryki;
          const SoAFloatPrecision drzk2 = drzjk * drzki;

          const SoAFloatPrecision drk2 = drxk2 + dryk2 + drzk2;

          const SoAFloatPrecision drijk2 = dri2 * drj2 * drk2;

          const SoAFloatPrecision dr2 = drij2 * drjk2 * drki2;
          const SoAFloatPrecision dr5 = dr2 * dr2 * std::sqrt(dr2);

          auto nu = _nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(type1ptr[i], type2ptr[j], type3ptr[k]);
          }
          const SoAFloatPrecision invdr5 = nu / dr5;
          const SoAFloatPrecision invdr53 = invdr5 * 3.0;

          const SoAFloatPrecision fxi =
              (drxjk * dri2 * (drj2 - drk2) + drxij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               drxki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          const SoAFloatPrecision fyi =
              (dryjk * dri2 * (drj2 - drk2) + dryij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               dryki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          const SoAFloatPrecision fzi =
              (drzjk * dri2 * (drj2 - drk2) + drzij * (drj2 * drk2 - drjk2 * drki2 + 5.0 * drijk2 / drij2) +
               drzki * (-drj2 * drk2 + drij2 * drjk2 - 5.0 * drijk2 / drki2)) *
              invdr53;
          fxiacc += fxi;
          fyiacc += fyi;
          fziacc += fzi;

          if (newton3) {
            const SoAFloatPrecision fxj =
                (drxki * drj2 * (drk2 - dri2) + drxij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                 drxjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
                invdr53;
            const SoAFloatPrecision fyj =
                (dryki * drj2 * (drk2 - dri2) + dryij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                 dryjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
                invdr53;
            const SoAFloatPrecision fzj =
                (drzki * drj2 * (drk2 - dri2) + drzij * (-dri2 * drk2 + drjk2 * drki2 - 5.0 * drijk2 / drij2) +
                 drzjk * (dri2 * drk2 - drij2 * drki2 + 5.0 * drijk2 / drjk2)) *
                invdr53;
            fxjacc += fxj;
            fyjacc += fyj;
            fzjacc += fzj;

            const auto nfxk = fxi + fxj;
            const auto nfyk = fyi + fyj;
            const auto nfzk = fzi + fzj;

            fx3ptr[k] -= nfxk;
            fy3ptr[k] -= nfyk;
            fz3ptr[k] -= nfzk;

            if constexpr (calculateGlobals) {
              virialSumX += fxj * xj - nfxk * xk;
              virialSumY += fyj * yj - nfyk * yk;
              virialSumZ += fzj * zj - nfzk * zk;
            }
          }

          if constexpr (calculateGlobals) {
            double potentialEnergy = (newton3 ? invdr53 : invdr5) * ((1.0 / 3.0) * dr2 - drijk2);
            potentialEnergySum += potentialEnergy;
            virialSumX += fxi * xi;
            virialSumY += fyi * yi;
            virialSumZ += fzi * zi;
          }
        }
        if constexpr (newton3) {
          fx2ptr[j] += fxjacc;
          fy2ptr[j] += fyjacc;
          fz2ptr[j] += fzjacc;
        }
      }
      fx1ptr[i] += fxiacc;
      fy1ptr[i] += fyiacc;
      fz1ptr[i] += fziacc;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();
      _aosThreadData[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadData[threadnum].virialSum[0] += virialSumX;
      _aosThreadData[threadnum].virialSum[1] += virialSumY;
      _aosThreadData[threadnum].virialSum[2] += virialSumZ;
    }
  }

 public:
  // clang-format off
  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorVerlet() is not implemented.");
  }

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param nu
   */
  void setParticleProperties(SoAFloatPrecision nu) { _nu = nu; }

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
   * @param molCType molecule C's type id
   * @param newton3 is newton3 applied.
   * @note The molecule types make no difference for AxilrodTellerFunctor, but are kept to have a consistent interface
   * for other functors where they may.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, size_t molCType, bool newton3) {
    //
    // Kernel: 56 = 18 (three dot products) + 9 (coefficients) + 29 (force calculation) sum
    // Adding to particle forces: 3
    // For Newton3: 29 (second force calculation) + 3 (adding force) + 6 (adding force to third p)
    // Total = 56 + 3 + ( 29 + 3 + 6 ) = 59 or 97
    return newton3 ? 97ul : 59ul;
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
      // Accumulate potential energy and virial values.
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
        _virialSum += _aosThreadData[i].virialSum;
      }

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
    autopas::utils::ExceptionHandler::exception("AxilrodTellerFunctor::SoAFunctorVerletImpl() is not implemented.");
  }

  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    // variables
    std::array<double, 3> virialSum;
    double potentialEnergySum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  const double _cutoffSquared;

  // not const because they might be reset through PPL
  double _nu = 0.0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

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

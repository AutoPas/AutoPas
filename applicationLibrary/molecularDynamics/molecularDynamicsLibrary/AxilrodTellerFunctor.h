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
 * The Axilrod-Teller potential
 * ---
 * The reference paper of Axilrod and Teller can be found here: https://doi.org/10.1063/1.1723844
 * \image html 3_body_sketch.png "Sketch of three particles that are used in the Axilrod-Teller Functor" width=400px
 *
 * The Axilrod-Teller potential is a model for the interactions of three molecules which appear when the van
 * der Waals forces are approximated to the third order. It is usually combined with a model for pairwise interaction as
 * e.g. the Lennard-Jones potential.
 *
 * \f[
 * U_{AT} = \nu \frac{3 \cos\gamma_1 \cos\gamma_2 \cos\gamma_3 + 1}{r_{12}^3 r_{23}^3 r_{31}^3}
 * \f]
 *
 * , where \f$r_{ij}\f$ is the distance between particles \f$i\f$ and \f$j\f$ and \f$\gamma_i\f$ is the angle between
 * the sides \f$r_{ij}\f$ and \f$r_{ik}\f$. \f$\nu\f$ is a material dependent parameter of the order \f$V\alpha^3\f$,
 * where \f$V\f$ is the ionization energy and \f$\alpha\f$ the polarizability.
 *
 * The cosines can also be expressed as:
 *
 * \f[
 *  \cos\gamma_1 = \frac{ \vec{r}_{12} \cdot \vec{r}_{13}}{|\vec{r}_{12}||\vec{r}_{13}|}
 * \f]
 *
 * , where \f$\vec{r}_{ij}\f$ is the vector from particle \f$i\f$ to particle \f$j\f$ (\f$i \longrightarrow j\f$ ).
 * It is calculated as \f$\vec{x}_j - \vec{x}_i\f$, where \f$\vec{x}_i\f$ is the position of particle \f$i\f$.
 *
 * Therefore, the potential can also be expressed as:
 *
 * \f[
 * U_{AT} = \nu\frac{-3 (\vec{r}_{12} \cdot \vec{r}_{31}) (\vec{r}_{12} \cdot \vec{r}_{23}) (\vec{r}_{31} \cdot
 * \vec{r}_{23}) + r_{12}^2 r_{23}^2 r_{31}^2}{r_{12}^5 r_{23}^5 r_{31}^5} \f]
 *
 * Note that we have \f$-3\f$ because we use the circular vectors \f$\vec{r}_ {12}, \vec{r}_ {23}, \vec{r}_ {31}\f$.
 *
 * The derivative can be calculated by applying the chain rule and leads to a resulting Force exerted on particle
 * \f$1\f$:
 *
 * \f[
 * \vec{F}_ {1} = - \frac{\partial U_ {AT}}{\partial \vec{x}_ 1}
 * \f]
 *
 * \f[
 * \vec{F}_ {1} = \frac{3}{r_ {12}^5 r_ {23}^5 r_ {31}^5}\cdot
 * \left[ \left( -5\frac{<>_ 1<>_ 2<>_ 3}{r_ {12}^2} - <>_ 1<>_ 3 + r_ {23}^2r_ {31}^2\right)\cdot \vec{r}_ {12}
 *          +\left( 5\frac{<>_ 1<>_ 2<>_ 3}{r_ {23}^2} + <>_ 1<>_ 3 - r_ {12}^2r_ {31}^2\right)\cdot\vec{r}_ {23}
 *          +\left( <>_ 2<>_ 3 - <>_ 2<>_ 1 \right)\cdot \vec{r}_ {31} \right]
 * \f]
 *
 * , where \f$<>_ 1=\vec{r}_ {12}\cdot\vec{r}_ {31}\f$ and so on. The terms are already ordered to show the contribution
 * from all three distance vectors.
 *
 * **Newton's third law**
 *
 * To apply Newton's third law, the force on particle \f$2\f$ needs to be calculated in a similar fashion as for
 * particle \f$1\f$. The force on particle \f$3\f$ can then be written as the negative sum of the other two forces:
 *
 * \f[
 * \vec{F}_3 = -(\vec{F}_1 + \vec{F}_2)
 * \f]
 *
 */

/**
 * A functor to handle Axilrod-Teller(-Muto) interactions between three particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 */
template <class Particle, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool countFLOPs = false>
class AxilrodTellerFunctor
    : public autopas::TriwiseFunctor<
          Particle, AxilrodTellerFunctor<Particle, useMixing, useNewton3, calculateGlobals, countFLOPs>> {
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
      : autopas::TriwiseFunctor<Particle,
                                AxilrodTellerFunctor<Particle, useMixing, useNewton3, calculateGlobals, countFLOPs>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadDataGlobals(),
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
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

  std::string getName() final { return "AxilrodTellerFunctorAutoVec"; }

  bool isRelevantForTuning() final { return true; }

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

    const auto threadnum = autopas::autopas_get_thread_num();

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
    }

    auto nu = _nu;
    if constexpr (useMixing) {
      nu = _PPLibrary->getMixingNu(i.getTypeId(), j.getTypeId(), k.getTypeId());
    }

    const auto displacementIJ = j.getR() - i.getR();
    const auto displacementJK = k.getR() - j.getR();
    const auto displacementKI = i.getR() - k.getR();

    const double distSquaredIJ = autopas::utils::ArrayMath::dot(displacementIJ, displacementIJ);
    const double distSquaredJK = autopas::utils::ArrayMath::dot(displacementJK, displacementJK);
    const double distSquaredKI = autopas::utils::ArrayMath::dot(displacementKI, displacementKI);

    // Check cutoff for every distance
    if (distSquaredIJ > _cutoffSquared or distSquaredJK > _cutoffSquared or distSquaredKI > _cutoffSquared) {
      return;
    }

    // Calculate prefactor
    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    const double factor = 3.0 * nu / allDistsTo5;

    // Dot products of both distance vectors going from one particle
    const double IJDotKI = autopas::utils::ArrayMath::dot(displacementIJ, displacementKI);
    const double IJDotJK = autopas::utils::ArrayMath::dot(displacementIJ, displacementJK);
    const double JKDotKI = autopas::utils::ArrayMath::dot(displacementJK, displacementKI);

    const double allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    const auto forceIDirectionJK = displacementJK * IJDotKI * (IJDotJK - JKDotKI);
    const auto forceIDirectionIJ =
        displacementIJ * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI + 5.0 * allDotProducts / distSquaredIJ);
    const auto forceIDirectionKI =
        displacementKI * (-IJDotJK * JKDotKI + distSquaredIJ * distSquaredJK - 5.0 * allDotProducts / distSquaredKI);

    const auto forceI = (forceIDirectionJK + forceIDirectionIJ + forceIDirectionKI) * factor;
    i.addF(forceI);

    auto forceJ = forceI;
    auto forceK = forceI;
    if (newton3) {
      const auto forceJDirectionKI = displacementKI * IJDotJK * (JKDotKI - IJDotKI);
      const auto forceJDirectionIJ =
          displacementIJ * (-IJDotKI * JKDotKI + distSquaredJK * distSquaredKI - 5.0 * allDotProducts / distSquaredIJ);
      const auto forceJDirectionJK =
          displacementJK * (IJDotKI * JKDotKI - distSquaredIJ * distSquaredKI + 5.0 * allDotProducts / distSquaredJK);

      forceJ = (forceJDirectionKI + forceJDirectionIJ + forceJDirectionJK) * factor;
      j.addF(forceJ);

      forceK = (forceI + forceJ) * (-1.0);
      k.addF(forceK);
    }

    if constexpr (countFLOPs) {
      if (newton3) {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsN3;
      } else {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsNoN3;
      }
    }

    if constexpr (calculateGlobals) {
      // Add 3 * potential energy to every owned particle of the interaction.
      // Division to the correct value is handled in endTraversal().
      const double potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);

      // Virial is calculated as f_i * r_i
      // see Thompson et al.: https://doi.org/10.1063/1.3245303
      const auto virialI = forceI * i.getR();
      if (i.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy3;
        _aosThreadDataGlobals[threadnum].virialSum += virialI;
      }
      // for non-newton3 particles j and/or k will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        const auto virialJ = forceJ * j.getR();
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy3;
        _aosThreadDataGlobals[threadnum].virialSum += virialJ;
      }
      if (newton3 and k.isOwned()) {
        const auto virialK = forceK * k.getR();
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy3;
        _aosThreadDataGlobals[threadnum].virialSum += virialK;
      }
      if constexpr (countFLOPs) {
        if (newton3) {
          ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsN3;
        } else {
          ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3;
        }
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
        // skip dummy particles
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
          // skip dummy particles
          continue;
        }

        SoAFloatPrecision fxjacc = 0.;
        SoAFloatPrecision fyjacc = 0.;
        SoAFloatPrecision fzjacc = 0.;

        const SoAFloatPrecision xj = xptr[j];
        const SoAFloatPrecision yj = yptr[j];
        const SoAFloatPrecision zj = zptr[j];

        // calculate distance i-j
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
            // skip dummy particles
            continue;
          }

          const SoAFloatPrecision xk = xptr[k];
          const SoAFloatPrecision yk = yptr[k];
          const SoAFloatPrecision zk = zptr[k];

          // calculate distance j-k
          const SoAFloatPrecision drxjk = xk - xj;
          const SoAFloatPrecision dryjk = yk - yj;
          const SoAFloatPrecision drzjk = zk - zj;

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          // calculate distance k-i
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

          // all particles are within cutoff -> calculate forces
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
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy3 = invdr53 * (dr2 - 3.0 * drijk2);
            potentialEnergySum += potentialEnergy3;
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
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
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

    // particle 1 always from soa1
    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        // skip dummy particles
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
          // skip dummy particles
          continue;
        }

        SoAFloatPrecision fxjacc = 0.;
        SoAFloatPrecision fyjacc = 0.;
        SoAFloatPrecision fzjacc = 0.;

        const SoAFloatPrecision xj = x1ptr[j];
        const SoAFloatPrecision yj = y1ptr[j];
        const SoAFloatPrecision zj = z1ptr[j];

        // calculate distance i-j
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

        // particle 3 always from soa2
        for (size_t k = 0; k < soa2.size(); ++k) {
          if (ownedState2ptr[k] == autopas::OwnershipState::dummy) {
            // skip dummy particles
            continue;
          }

          const SoAFloatPrecision xk = x2ptr[k];
          const SoAFloatPrecision yk = y2ptr[k];
          const SoAFloatPrecision zk = z2ptr[k];

          // calculate distance j-k
          const SoAFloatPrecision drxjk = xk - xj;
          const SoAFloatPrecision dryjk = yk - yj;
          const SoAFloatPrecision drzjk = zk - zj;

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          // calculate distance k-i
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

          // all particles are within cutoff -> calculate forces
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
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy3 = (newton3 ? 3.0 : 2.0) * invdr53 * (dr2 - 3.0 * drijk2);
            potentialEnergySum += potentialEnergy3;
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
          // skip dummy particles
          continue;
        }

        SoAFloatPrecision fxjacc = 0.;
        SoAFloatPrecision fyjacc = 0.;
        SoAFloatPrecision fzjacc = 0.;

        const SoAFloatPrecision xj = x2ptr[j];
        const SoAFloatPrecision yj = y2ptr[j];
        const SoAFloatPrecision zj = z2ptr[j];

        // calculate distance i-j
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

        // particle 3 always from soa 2
        for (size_t k = 0; k < j; ++k) {
          if (ownedState2ptr[k] == autopas::OwnershipState::dummy) {
            // skip dummy particles
            continue;
          }

          const SoAFloatPrecision xk = x2ptr[k];
          const SoAFloatPrecision yk = y2ptr[k];
          const SoAFloatPrecision zk = z2ptr[k];

          // calculate distance j-k
          const SoAFloatPrecision drxjk = xk - xj;
          const SoAFloatPrecision dryjk = yk - yj;
          const SoAFloatPrecision drzjk = zk - zj;

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          // calculate distance k-i
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

          // all particles are within cutoff -> calculate forces
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
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy3 = (newton3 ? 3.0 : 1.0) * invdr53 * (dr2 - 3.0 * drijk2);
            potentialEnergySum += potentialEnergy3;
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
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
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
        // skip dummy particles
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
          // skip dummy particles
          continue;
        }

        SoAFloatPrecision fxjacc = 0.;
        SoAFloatPrecision fyjacc = 0.;
        SoAFloatPrecision fzjacc = 0.;

        const SoAFloatPrecision xj = x2ptr[j];
        const SoAFloatPrecision yj = y2ptr[j];
        const SoAFloatPrecision zj = z2ptr[j];

        // calculate distance i-j
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
            // skip dummy particles
            continue;
          }

          const SoAFloatPrecision xk = x3ptr[k];
          const SoAFloatPrecision yk = y3ptr[k];
          const SoAFloatPrecision zk = z3ptr[k];

          // calculate distance j-k
          const SoAFloatPrecision drxjk = xk - xj;
          const SoAFloatPrecision dryjk = yk - yj;
          const SoAFloatPrecision drzjk = zk - zj;

          const SoAFloatPrecision drxjk2 = drxjk * drxjk;
          const SoAFloatPrecision dryjk2 = dryjk * dryjk;
          const SoAFloatPrecision drzjk2 = drzjk * drzjk;

          const SoAFloatPrecision drjk2 = drxjk2 + dryjk2 + drzjk2;

          // calculate distance k-i
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

          // all particles are within cutoff -> calculate forces
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

            const SoAFloatPrecision nfxk = fxi + fxj;
            const SoAFloatPrecision nfyk = fyi + fyj;
            const SoAFloatPrecision nfzk = fzi + fzj;

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
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy3 = (newton3 ? 3.0 : 1.0) * invdr53 * (dr2 - 3.0 * drijk2);
            potentialEnergySum += potentialEnergy3;
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
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }


 public:
  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param nu The Axilrod-Teller potential parameter
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
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() final {
    _potentialEnergySum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadDataGlobals.size(); ++i) {
      _aosThreadDataGlobals[i].setZero();
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
      for (const auto &data : _aosThreadDataGlobals) {
        _potentialEnergySum += data.potentialEnergySum;
        _virialSum += data.virialSum;
      }

      // For each interaction, we added the full contribution for all three particles. Divide by 3 here, so that each
      // contribution is only counted once per triplet.
      _potentialEnergySum /= 3.;

      // Additionally, we have always calculated 3*potentialEnergy, so we divide by 3 again.
      _potentialEnergySum /= 3.;

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

  /**
   * Gets the number of useful FLOPs.
   *
   * For the three distance squared calculations, this is:
   * - Displacement: 3
   * - DistanceSquared: 5
   * - Total: 8 * 3 = 24
   *
   * For the force kernel, this is:
   * - calculation of prefactor: 7
   * - dot products: 3 * 5 = 15
   * - all dot products: 2
   * - forceIDirectionJK: 5
   * - forceIDirectionIJ: 9
   * - forceIDirectionKI: 9
   * - add force vectors and multiply: 9
   * - add force to mol i: 3
   * - If N3:
   * - forceJDirectionKI: 5
   * - forceJDirectionIJ: 9
   * - forceJDirectionJK: 9
   * - add force vectors and multiply: 9
   * - add force to mol j: 3
   * - sum forceK: 3 (don't count multiplication with -1.0)
   * - add force to mol k: 3
   * - Total: 59 without n3, 100 with n3
   *
   * For the globals calculation, this is:
   * - potential: 3
   * - virial: 3 without n3, 9 with n3
   * - accumulation: 4 without n3, 12 with n3
   * - Total: 10 without n3, 24 with n3
   *
   * @return number of FLOPs since initTraversal() is called.
   */
  [[nodiscard]] size_t getNumFLOPs() const override {
    if constexpr (countFLOPs) {
      const size_t numDistCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });
      const size_t numGlobalCalcsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsN3; });
      const size_t numGlobalCalcsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsNoN3; });

      constexpr size_t numFLOPsPerDistanceCall = 24;
      constexpr size_t numFLOPsPerN3KernelCall = 100;
      constexpr size_t numFLOPsPerNoN3KernelCall = 59;
      constexpr size_t numFLOPsPerN3GlobalCalc = 24;
      constexpr size_t numFLOPsPerNoN3GlobalCalc = 10;

      return numDistCallsAcc * numFLOPsPerDistanceCall + numKernelCallsN3Acc * numFLOPsPerN3KernelCall +
             numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall + numGlobalCalcsN3Acc * numFLOPsPerN3GlobalCalc +
             numGlobalCalcsNoN3Acc * numFLOPsPerNoN3GlobalCalc;
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<size_t>::max();
    }
  }

  [[nodiscard]] double getHitRate() const override {
    if constexpr (countFLOPs) {
      const size_t numDistCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });

      return (static_cast<double>(numKernelCallsNoN3Acc) + static_cast<double>(numKernelCallsN3Acc)) /
             (static_cast<double>(numDistCallsAcc));
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<double>::quiet_NaN();
    }
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
  class AoSThreadDataGlobals {
   public:
    AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
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

  /**
   * This class stores internal data for FLOP counters for each thread. Make sure that this data has proper size, i.e.
   * k*64 Bytes!
   * The FLOP count and HitRate are not counted/calculated directly, but through helper counters (numKernelCallsNoN3,
   * numKernelCallsN3, numDistCalls, numGlobalCalcs) to reduce computational cost in the functors themselves and to
   * improve maintainability (e.g. if the cost of a kernel call changes).
   */
  class AoSThreadDataFLOPs {
   public:
    AoSThreadDataFLOPs() : __remainingTo64{} {}

    /**
     * Set all counters to zero.
     */
    void setZero() {
      numKernelCallsNoN3 = 0;
      numKernelCallsN3 = 0;
      numDistCalls = 0;
      numGlobalCalcsN3 = 0;
      numGlobalCalcsNoN3 = 0;
    }

    /**
     * Number of calls to Lennard-Jones Kernel with newton3 disabled.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numKernelCallsNoN3 = 0;

    /**
     * Number of calls to Lennard-Jones Kernel with newton3 enabled.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numKernelCallsN3 = 0;

    /**
     * Number of distance calculations.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numDistCalls = 0;

    /**
     * Counter for the number of times the globals have been calculated with newton3 enabled.
     */
    size_t numGlobalCalcsN3 = 0;

    /**
     * Counter for the number of times the globals have been calculated without newton3 enabled.
     */
    size_t numGlobalCalcsNoN3 = 0;

   private:
    /**
     * dummy parameter to get the right size (64 bytes)
     */
    double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
  };

  // make sure of the size of AoSThreadDataGlobals
  static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
  static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

  const double _cutoffSquared;

  // Parameter of the Axilrod-Teller potential
  // not const because they might be reset through PPL
  double _nu = 0.0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals;
  std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;
};
}  // namespace mdLib

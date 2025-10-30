/**
 * @file AxilrodTellerMutoFunctor.h
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
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib {

/**
 * The Axilrod-Teller-Muto potential
 * ---
 * The reference paper of Axilrod and Teller can be found here: https://doi.org/10.1063/1.1723844
 * \image html 3_body_sketch.png "Sketch of three particles that are used in the Axilrod-Teller-Muto Functor"
 * width=400px
 *
 * The Axilrod-Teller-Muto potential is a model for the interactions of three molecules which appear when the van
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
 * A functor to handle Axilrod-Teller-Muto interactions between three particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle_T The type of particle.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 */
template <class Particle_T, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool countFLOPs = false>
class AxilrodTellerMutoFunctor
    : public autopas::TriwiseFunctor<
          Particle_T, AxilrodTellerMutoFunctor<Particle_T, useMixing, useNewton3, calculateGlobals, countFLOPs>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  AxilrodTellerMutoFunctor() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit AxilrodTellerMutoFunctor(double cutoff, void * /*dummy*/)
      : autopas::TriwiseFunctor<
            Particle_T, AxilrodTellerMutoFunctor<Particle_T, useMixing, useNewton3, calculateGlobals, countFLOPs>>(
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
      AutoPasLog(DEBUG, "Using AxilrodTellerMutoFunctor with countFLOPs is not tested for SoA datalayout.");
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
  explicit AxilrodTellerMutoFunctor(double cutoff) : AxilrodTellerMutoFunctor(cutoff, nullptr) {
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
  explicit AxilrodTellerMutoFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : AxilrodTellerMutoFunctor(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "AxilrodTellerMutoFunctorAutoVec"; }

  bool isRelevantForTuning() final { return true; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle_T &i, Particle_T &j, Particle_T &k, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy() or k.isDummy()) {
      return;
    }

    const auto threadnum = autopas::autopas_get_thread_num();

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numTripletsCount;
      _aosThreadDataFLOPs[threadnum].numDistCalls += 3;
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
      if (i.isOwned()) {
        const auto virialI = forceI * (displacementKI - displacementIJ);
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy3;
        _aosThreadDataGlobals[threadnum].virialSum += virialI;
      }
      // for non-newton3 particles j and/or k will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        const auto virialJ = forceJ * (displacementIJ - displacementJK);
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy3;
        _aosThreadDataGlobals[threadnum].virialSum += virialJ;
      }
      if (newton3 and k.isOwned()) {
        const auto virialK = forceK * (displacementJK - displacementKI);
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

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.size() <= 2) return;

    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;  // Note: This is not the potential energy but some fixed multiple of it.
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numTripletsCountingSum = 0;
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numGlobalCalcsSum = 0;

    const SoAFloatPrecision const_nu = _nu;
    const size_t soaSize = soa.size();
    if constexpr (countFLOPs) {
      numTripletsCountingSum = soaSize * (soaSize - 1) * (soaSize - 2) / 6;
    }

    for (unsigned int i = 0; i < soaSize - 2; ++i) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }
      SoAFloatPrecision fXAccI = 0.;
      SoAFloatPrecision fYAccI = 0.;
      SoAFloatPrecision fZAccI = 0.;

      const SoAFloatPrecision xi = xptr[i];
      const SoAFloatPrecision yi = yptr[i];
      const SoAFloatPrecision zi = zptr[i];

      for (unsigned int j = i + 1; j < soaSize - 1; ++j) {
        const auto ownedStateJ = ownedStatePtr[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        const SoAFloatPrecision xj = xptr[j];
        const SoAFloatPrecision yj = yptr[j];
        const SoAFloatPrecision zj = zptr[j];

        const SoAFloatPrecision distXIJ = xj - xi;
        const SoAFloatPrecision distYIJ = yj - yi;
        const SoAFloatPrecision distZIJ = zj - zi;
        const SoAFloatPrecision distSquaredIJ = distXIJ * distXIJ + distYIJ * distYIJ + distZIJ * distZIJ;

        if constexpr (countFLOPs) {
          ++numDistanceCalculationSum;
        }
        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        SoAFloatPrecision fXAccJ = 0.;
        SoAFloatPrecision fYAccJ = 0.;
        SoAFloatPrecision fZAccJ = 0.;

        for (unsigned int k = j + 1; k < soaSize; ++k) {
          const auto ownedStateK = ownedStatePtr[k];
          if (ownedStateK == autopas::OwnershipState::dummy) {
            continue;
          }

          SoAFloatPrecision nu = const_nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(typeptr[i], typeptr[j], typeptr[k]);
          }

          const SoAFloatPrecision xk = xptr[k];
          const SoAFloatPrecision yk = yptr[k];
          const SoAFloatPrecision zk = zptr[k];

          const SoAFloatPrecision distXJK = xk - xj;
          const SoAFloatPrecision distYJK = yk - yj;
          const SoAFloatPrecision distZJK = zk - zj;
          const SoAFloatPrecision distSquaredJK = distXJK * distXJK + distYJK * distYJK + distZJK * distZJK;

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredJK > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision distXKI = xi - xk;
          const SoAFloatPrecision distYKI = yi - yk;
          const SoAFloatPrecision distZKI = zi - zk;
          const SoAFloatPrecision distSquaredKI = distXKI * distXKI + distYKI * distYKI + distZKI * distZKI;

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredKI > cutoffSquared) {
            continue;
          }

          SoAFloatPrecision forceIX, forceIY, forceIZ;
          SoAFloatPrecision forceJX, forceJY, forceJZ;
          SoAFloatPrecision factor, allDotProducts, allDistsSquared;
          SoAKernelN3(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, distXKI, distYKI, distZKI, distSquaredIJ,
                      distSquaredJK, distSquaredKI, nu, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ, factor,
                      allDotProducts, allDistsSquared);

          fXAccI += forceIX;
          fYAccI += forceIY;
          fZAccI += forceIZ;

          fXAccJ += forceJX;
          fYAccJ += forceJY;
          fZAccJ += forceJZ;

          const SoAFloatPrecision forceKX = -(forceIX + forceJX);
          const SoAFloatPrecision forceKY = -(forceIY + forceJY);
          const SoAFloatPrecision forceKZ = -(forceIZ + forceJZ);

          fxptr[k] += forceKX;
          fyptr[k] += forceKY;
          fzptr[k] += forceKZ;

          if constexpr (countFLOPs) {
            ++numKernelCallsN3Sum;
          }

          if constexpr (calculateGlobals) {
            const SoAFloatPrecision potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);
            if (ownedStateI == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy3;
              virialSumX += forceIX * (distXKI - distXIJ);
              virialSumY += forceIY * (distYKI - distYIJ);
              virialSumZ += forceIZ * (distZKI - distZIJ);
            }
            if (ownedStateJ == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy3;
              virialSumX += forceJX * (distXIJ - distXJK);
              virialSumY += forceJY * (distYIJ - distYJK);
              virialSumZ += forceJZ * (distZIJ - distZJK);
            }
            if (ownedStateK == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy3;
              virialSumX += forceKX * (distXJK - distXKI);
              virialSumY += forceKY * (distYJK - distYKI);
              virialSumZ += forceKZ * (distZJK - distZKI);
            }

            if constexpr (countFLOPs) {
              ++numGlobalCalcsSum;
            }
          }
        }
        fxptr[j] += fXAccJ;
        fyptr[j] += fYAccJ;
        fzptr[j] += fZAccJ;
      }
      fxptr[i] += fXAccI;
      fyptr[i] += fYAccI;
      fzptr[i] += fZAccI;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numTripletsCount += numTripletsCountingSum;
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsSum;  // Always N3 in Single SoAFunctor
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final {
    if (soa1.size() == 0 || soa2.size() == 0) return;

    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

  void SoAFunctorTriple(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                        autopas::SoAView<SoAArraysType> soa3, bool newton3) final {
    if (soa1.size() == 0 || soa2.size() == 0 || soa3.size() == 0) return;

    if (newton3) {
      SoAFunctorTripleImpl<true>(soa1, soa2, soa3);
    } else {
      SoAFunctorTripleImpl<false>(soa1, soa2, soa3);
    }
  }

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param nu The Axilrod-Teller-Muto potential parameter
   */
  void setParticleProperties(SoAFloatPrecision nu) { _nu = nu; }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 9>{Particle_T::AttributeNames::id,
                                                              Particle_T::AttributeNames::posX,
                                                              Particle_T::AttributeNames::posY,
                                                              Particle_T::AttributeNames::posZ,
                                                              Particle_T::AttributeNames::forceX,
                                                              Particle_T::AttributeNames::forceY,
                                                              Particle_T::AttributeNames::forceZ,
                                                              Particle_T::AttributeNames::typeId,
                                                              Particle_T::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle_T::AttributeNames, 6>{
        Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::typeId, Particle_T::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
        Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
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
      _virialSum *= (1. / 3.);

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
   * - virial: 6 without n3, 18 with n3
   * - accumulation: 4 without n3, 12 with n3
   * - Total: 13 without n3, 33 with n3
   * @note The exact number of FLOPs can slightly deviate due to SoAFunctorPair sometimes doing a N3 kernel call for
   * only 2 out of three particles.
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

      constexpr size_t numFLOPsPerDistanceCall = 8;
      constexpr size_t numFLOPsPerN3KernelCall = 100;
      constexpr size_t numFLOPsPerNoN3KernelCall = 59;
      constexpr size_t numFLOPsPerN3GlobalCalc = 33;
      constexpr size_t numFLOPsPerNoN3GlobalCalc = 13;

      return numDistCallsAcc * numFLOPsPerDistanceCall + numKernelCallsN3Acc * numFLOPsPerN3KernelCall +
             numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall + numGlobalCalcsN3Acc * numFLOPsPerN3GlobalCalc +
             numGlobalCalcsNoN3Acc * numFLOPsPerNoN3GlobalCalc;
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<size_t>::max();
    }
  }

  /**
   * @copydoc autopas::Functor::getHitRate()
   * @note Specifically, the hitrate for this functor is defined as: (# kernel calls) / (# possible triplets)
   */
  [[nodiscard]] double getHitRate() const override {
    if constexpr (countFLOPs) {
      const size_t numTripletsCount =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numTripletsCount; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });

      return (static_cast<double>(numKernelCallsNoN3Acc) + static_cast<double>(numKernelCallsN3Acc)) /
             (static_cast<double>(numTripletsCount));
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

 private:
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xptr1 = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr1 = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr1 = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict xptr2 = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr2 = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr2 = soa2.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fxptr1 = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr1 = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzptr1 = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fxptr2 = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr2 = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzptr2 = soa2.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle_T::AttributeNames::typeId>();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;  // Note: This is not the potential energy but some fixed multiple of it.
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numTripletsCountingSum = 0;
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const SoAFloatPrecision const_nu = _nu;

    size_t soa1Size = soa1.size();
    size_t soa2Size = soa2.size();
    if constexpr (countFLOPs) {
      numTripletsCountingSum = (soa1Size * soa2Size * (soa1Size + soa2Size - 2)) / 2.;
    }

    // Precompute distances between soa1 and soa2
    std::vector<std::array<SoAFloatPrecision, 4>, autopas::AlignedAllocator<std::array<SoAFloatPrecision, 4>>>
        soa1Soa2Dists(soa1Size * soa2Size);
    for (unsigned int i = 0; i < soa1Size; ++i) {
      for (unsigned int j = 0; j < soa2Size; ++j) {
        const SoAFloatPrecision distXIJ = xptr2[j] - xptr1[i];
        const SoAFloatPrecision distYIJ = yptr2[j] - yptr1[i];
        const SoAFloatPrecision distZIJ = zptr2[j] - zptr1[i];
        const SoAFloatPrecision distSquaredIJ = distXIJ * distXIJ + distYIJ * distYIJ + distZIJ * distZIJ;
        soa1Soa2Dists[j + i * soa2Size] = std::array<SoAFloatPrecision, 4>{distXIJ, distYIJ, distZIJ, distSquaredIJ};
      }
    }
    if constexpr (countFLOPs) {
      numDistanceCalculationSum += soa1Size * soa2Size;
    }

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }
      SoAFloatPrecision fXAccI = 0.;
      SoAFloatPrecision fYAccI = 0.;
      SoAFloatPrecision fZAccI = 0.;

      const SoAFloatPrecision xi = xptr1[i];
      const SoAFloatPrecision yi = yptr1[i];
      const SoAFloatPrecision zi = zptr1[i];

      // CASE: Particle i is in soa1, j and k are both in soa2
      for (unsigned int j = 0; j < soa2.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr2[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        const auto &[distXIJ, distYIJ, distZIJ, distSquaredIJ] = soa1Soa2Dists[i * soa2Size + j];

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        for (unsigned int k = j + 1; k < soa2.size(); ++k) {
          const auto ownedStateK = ownedStatePtr2[k];
          if (ownedStateK == autopas::OwnershipState::dummy) {
            continue;
          }

          SoAFloatPrecision nu = const_nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(typeptr1[i], typeptr2[j], typeptr2[k]);
          }

          // Actual return values are IK, not KI
          const auto &[distXIK, distYIK, distZIK, distSquaredKI] = soa1Soa2Dists[i * soa2Size + k];

          // Due to low 3-body hitrate, mostly better than masking
          if (distSquaredKI > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision distXJK = xptr2[k] - xptr2[j];
          const SoAFloatPrecision distYJK = yptr2[k] - yptr2[j];
          const SoAFloatPrecision distZJK = zptr2[k] - zptr2[j];
          const SoAFloatPrecision distSquaredJK = distXJK * distXJK + distYJK * distYJK + distZJK * distZJK;

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredJK > cutoffSquared) {
            continue;
          }

          if constexpr (not newton3) {
            // Compute only force I without Newton3
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision factor, allDotProducts, allDistsSquared;
            SoAKernelNoN3(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, -distXIK, -distYIK, -distZIK,
                          distSquaredIJ, distSquaredJK, distSquaredKI, nu, forceIX, forceIY, forceIZ, factor,
                          allDotProducts, allDistsSquared);

            fXAccI += forceIX;
            fYAccI += forceIY;
            fZAccI += forceIZ;

            if constexpr (countFLOPs) {
              ++numKernelCallsNoN3Sum;
            }
            if constexpr (calculateGlobals) {
              const SoAFloatPrecision potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);

              if (ownedStateI == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceIX * (-distXIK - distXIJ);
                virialSumY += forceIY * (-distYIK - distYIJ);
                virialSumZ += forceIZ * (-distZIK - distZIJ);
              }
              if constexpr (countFLOPs) {
                ++numGlobalCalcsNoN3Sum;
              }
            }
          } else {
            // Compute all forces with Newton3
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision forceJX, forceJY, forceJZ;
            SoAFloatPrecision factor, allDotProducts, allDistsSquared;
            SoAKernelN3(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, -distXIK, -distYIK, -distZIK,
                        distSquaredIJ, distSquaredJK, distSquaredKI, nu, forceIX, forceIY, forceIZ, forceJX, forceJY,
                        forceJZ, factor, allDotProducts, allDistsSquared);

            fXAccI += forceIX;
            fYAccI += forceIY;
            fZAccI += forceIZ;

            fxptr2[j] += forceJX;
            fyptr2[j] += forceJY;
            fzptr2[j] += forceJZ;

            const SoAFloatPrecision forceKX = -(forceIX + forceJX);
            const SoAFloatPrecision forceKY = -(forceIY + forceJY);
            const SoAFloatPrecision forceKZ = -(forceIZ + forceJZ);

            fxptr2[k] += forceKX;
            fyptr2[k] += forceKY;
            fzptr2[k] += forceKZ;

            if constexpr (countFLOPs) {
              ++numKernelCallsN3Sum;
            }

            if constexpr (calculateGlobals) {
              const SoAFloatPrecision potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);

              if (ownedStateI == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceIX * (-distXIK - distXIJ);
                virialSumY += forceIY * (-distYIK - distYIJ);
                virialSumZ += forceIZ * (-distZIK - distZIJ);
              }
              if (ownedStateJ == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceJX * (distXIJ - distXJK);
                virialSumY += forceJY * (distYIJ - distYJK);
                virialSumZ += forceJZ * (distZIJ - distZJK);
              }
              if (ownedStateK == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceKX * (distXJK + distXIK);
                virialSumY += forceKY * (distYJK + distYIK);
                virialSumZ += forceKZ * (distZJK + distZIK);
              }
              if constexpr (countFLOPs) {
                ++numGlobalCalcsN3Sum;
              }
            }
          }
        }
      }

      // CASE: Particle i and j are both in soa1, k is in soa2
      for (unsigned int j = i + 1; j < soa1.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr1[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        SoAFloatPrecision fXAccJ = 0.;
        SoAFloatPrecision fYAccJ = 0.;
        SoAFloatPrecision fZAccJ = 0.;

        const SoAFloatPrecision distXIJ = xptr1[j] - xi;
        const SoAFloatPrecision distYIJ = yptr1[j] - yi;
        const SoAFloatPrecision distZIJ = zptr1[j] - zi;
        const SoAFloatPrecision distSquaredIJ = distXIJ * distXIJ + distYIJ * distYIJ + distZIJ * distZIJ;

        if constexpr (countFLOPs) {
          ++numDistanceCalculationSum;
        }
        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        for (unsigned int k = 0; k < soa2.size(); ++k) {
          const auto ownedStateK = ownedStatePtr2[k];
          if (ownedStateK == autopas::OwnershipState::dummy) {
            continue;
          }

          SoAFloatPrecision nu = const_nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(typeptr1[i], typeptr1[j], typeptr2[k]);
          }

          const auto &[distXJK, distYJK, distZJK, distSquaredJK] = soa1Soa2Dists[j * soa2Size + k];

          // Due to low 3-body hitrate, mostly better than masking
          if (distSquaredJK > cutoffSquared) {
            continue;
          }

          const auto &[distXIK, distYIK, distZIK, distSquaredKI] = soa1Soa2Dists[i * soa2Size + k];

          if (distSquaredKI > cutoffSquared) {
            continue;
          }

          SoAFloatPrecision forceIX, forceIY, forceIZ;
          SoAFloatPrecision forceJX, forceJY, forceJZ;
          SoAFloatPrecision factor, allDotProducts, allDistsSquared;
          SoAKernelN3(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, -distXIK, -distYIK, -distZIK, distSquaredIJ,
                      distSquaredJK, distSquaredKI, nu, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ, factor,
                      allDotProducts, allDistsSquared);

          fXAccI += forceIX;
          fYAccI += forceIY;
          fZAccI += forceIZ;

          fxptr1[j] += forceJX;
          fyptr1[j] += forceJY;
          fzptr1[j] += forceJZ;

          if constexpr (countFLOPs) {
            ++numKernelCallsN3Sum;
          }
          if constexpr (newton3) {
            const SoAFloatPrecision forceKX = -(forceIX + forceJX);
            const SoAFloatPrecision forceKY = -(forceIY + forceJY);
            const SoAFloatPrecision forceKZ = -(forceIZ + forceJZ);
            fxptr2[k] += forceKX;
            fyptr2[k] += forceKY;
            fzptr2[k] += forceKZ;

            if constexpr (calculateGlobals) {
              const SoAFloatPrecision potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);
              if (ownedStateK == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceKX * (distXJK + distXIK);
                virialSumY += forceKY * (distYJK + distYIK);
                virialSumZ += forceKZ * (distZJK + distZIK);
              }
            }
          }

          if constexpr (calculateGlobals) {
            const SoAFloatPrecision potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);
            if (ownedStateI == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy3;
              virialSumX += forceIX * (-distXIK - distXIJ);
              virialSumY += forceIY * (-distYIK - distYIJ);
              virialSumZ += forceIZ * (-distZIK - distZIJ);
            }
            if (ownedStateJ == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy3;
              virialSumX += forceJX * (distXIJ - distXJK);
              virialSumY += forceJY * (distYIJ - distYJK);
              virialSumZ += forceJZ * (distZIJ - distZJK);
            }
            if constexpr (countFLOPs) {
              ++numGlobalCalcsN3Sum;
            }
          }
        }
        fxptr1[j] += fXAccJ;
        fyptr1[j] += fYAccJ;
        fzptr1[j] += fZAccJ;
      }
      fxptr1[i] += fXAccI;
      fyptr1[i] += fYAccI;
      fzptr1[i] += fZAccI;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numTripletsCount += numTripletsCountingSum;
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

  template <bool newton3>
  void SoAFunctorTripleImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                            autopas::SoAView<SoAArraysType> soa3) {
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xptr1 = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr1 = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr1 = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict xptr2 = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr2 = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr2 = soa2.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict xptr3 = soa3.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr3 = soa3.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr3 = soa3.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr3 = soa3.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fxptr1 = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr1 = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzptr1 = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fxptr2 = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr2 = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzptr2 = soa2.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fyptr3 = soa3.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fxptr3 = soa3.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fzptr3 = soa3.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr3 = soa3.template begin<Particle_T::AttributeNames::typeId>();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;  // Note: This is not the potential energy but some fixed multiple of it.
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numTripletsCountingSum = 0;
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const SoAFloatPrecision const_nu = _nu;

    const auto soa1Size = soa1.size();
    const auto soa2Size = soa2.size();
    const auto soa3Size = soa3.size();
    if constexpr (countFLOPs) {
      numTripletsCountingSum = soa1Size * soa2Size * soa3Size;
    }

    // Precompute distances between soa2 and soa3
    std::vector<std::array<SoAFloatPrecision, 4>, autopas::AlignedAllocator<std::array<SoAFloatPrecision, 4>>> jkDists(
        soa2Size * soa3Size);
    for (unsigned int j = 0; j < soa2Size; ++j) {
      for (unsigned int k = 0; k < soa3Size; ++k) {
        const SoAFloatPrecision distXJK = xptr3[k] - xptr2[j];
        const SoAFloatPrecision distYJK = yptr3[k] - yptr2[j];
        const SoAFloatPrecision distZJK = zptr3[k] - zptr2[j];
        const SoAFloatPrecision distSquaredJK = distXJK * distXJK + distYJK * distYJK + distZJK * distZJK;
        jkDists[k + j * soa3Size] = std::array<SoAFloatPrecision, 4>{distXJK, distYJK, distZJK, distSquaredJK};
      }
    }
    if constexpr (countFLOPs) {
      numDistanceCalculationSum += soa2Size * soa3Size;
    }

    for (unsigned int i = 0; i < soa1Size; ++i) {
      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision xi = xptr1[i];
      const SoAFloatPrecision yi = yptr1[i];
      const SoAFloatPrecision zi = zptr1[i];

      SoAFloatPrecision fXAccI = 0.;
      SoAFloatPrecision fYAccI = 0.;
      SoAFloatPrecision fZAccI = 0.;

      for (unsigned int j = 0; j < soa2Size; ++j) {
        const auto ownedStateJ = ownedStatePtr2[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        const SoAFloatPrecision xj = xptr2[j];
        const SoAFloatPrecision yj = yptr2[j];
        const SoAFloatPrecision zj = zptr2[j];

        const SoAFloatPrecision distXIJ = xj - xi;
        const SoAFloatPrecision distYIJ = yj - yi;
        const SoAFloatPrecision distZIJ = zj - zi;
        const SoAFloatPrecision distSquaredIJ = distXIJ * distXIJ + distYIJ * distYIJ + distZIJ * distZIJ;

        if constexpr (countFLOPs) {
          ++numDistanceCalculationSum;
        }
        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        for (unsigned int k = 0; k < soa3Size; ++k) {
          const auto ownedStateK = ownedStatePtr3[k];
          if (ownedStateK == autopas::OwnershipState::dummy) {
            continue;
          }

          SoAFloatPrecision nu = const_nu;
          if constexpr (useMixing) {
            nu = _PPLibrary->getMixingNu(typeptr1[i], typeptr2[j], typeptr3[k]);
          }

          const auto &[distXJK, distYJK, distZJK, distSquaredJK] = jkDists[j * soa3Size + k];

          if (distSquaredJK > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision distXKI = xi - xptr3[k];
          const SoAFloatPrecision distYKI = yi - yptr3[k];
          const SoAFloatPrecision distZKI = zi - zptr3[k];
          const SoAFloatPrecision distSquaredKI = distXKI * distXKI + distYKI * distYKI + distZKI * distZKI;

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredKI > cutoffSquared) {
            continue;
          }

          if constexpr (not newton3) {
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision factor, allDotProducts, allDistsSquared;
            SoAKernelNoN3(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, distXKI, distYKI, distZKI,
                          distSquaredIJ, distSquaredJK, distSquaredKI, nu, forceIX, forceIY, forceIZ, factor,
                          allDotProducts, allDistsSquared);

            fXAccI += forceIX;
            fYAccI += forceIY;
            fZAccI += forceIZ;

            if constexpr (countFLOPs) {
              ++numKernelCallsNoN3Sum;
            }
            if constexpr (calculateGlobals) {
              const SoAFloatPrecision potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);

              if (ownedStateI == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceIX * (distXKI - distXIJ);
                virialSumY += forceIY * (distYKI - distYIJ);
                virialSumZ += forceIZ * (distZKI - distZIJ);
              }

              if constexpr (countFLOPs) {
                ++numGlobalCalcsNoN3Sum;
              }
            }

          } else {
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision forceJX, forceJY, forceJZ;
            SoAFloatPrecision factor, allDotProducts, allDistsSquared;
            SoAKernelN3(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, distXKI, distYKI, distZKI, distSquaredIJ,
                        distSquaredJK, distSquaredKI, nu, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ, factor,
                        allDotProducts, allDistsSquared);

            fXAccI += forceIX;
            fYAccI += forceIY;
            fZAccI += forceIZ;

            fxptr2[j] += forceJX;
            fyptr2[j] += forceJY;
            fzptr2[j] += forceJZ;

            const SoAFloatPrecision forceKX = -(forceIX + forceJX);
            const SoAFloatPrecision forceKY = -(forceIY + forceJY);
            const SoAFloatPrecision forceKZ = -(forceIZ + forceJZ);

            fxptr3[k] += forceKX;
            fyptr3[k] += forceKY;
            fzptr3[k] += forceKZ;

            if constexpr (countFLOPs) {
              ++numKernelCallsN3Sum;
            }
            if constexpr (calculateGlobals) {
              const SoAFloatPrecision potentialEnergy3 = factor * (allDistsSquared - 3.0 * allDotProducts);
              if (ownedStateI == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceIX * (distXKI - distXIJ);
                virialSumY += forceIY * (distYKI - distYIJ);
                virialSumZ += forceIZ * (distZKI - distZIJ);
              }
              if (ownedStateJ == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceJX * (distXIJ - distXJK);
                virialSumY += forceJY * (distYIJ - distYJK);
                virialSumZ += forceJZ * (distZIJ - distZJK);
              }
              if (ownedStateK == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceKX * (distXJK - distXKI);
                virialSumY += forceKY * (distYJK - distYKI);
                virialSumZ += forceKZ * (distZJK - distZKI);
              }
              if constexpr (countFLOPs) {
                ++numGlobalCalcsN3Sum;
              }
            }
          }
        }
      }
      fxptr1[i] += fXAccI;
      fyptr1[i] += fYAccI;
      fzptr1[i] += fZAccI;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numTripletsCount += numTripletsCountingSum;
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    autopas::utils::ExceptionHandler::exception("AxilrodTellerMutoFunctor::SoAFunctorVerletImpl() is not implemented.");
  }

  /**
   * Inline helper to compute force components for particle I.
   * Returns tuple of (forceIX, forceIY, forceIZ, factor, allDotProducts, allDistsSquared)
   */
  __attribute__((always_inline)) inline void SoAKernelNoN3(
      const SoAFloatPrecision &distXIJ, const SoAFloatPrecision &distYIJ, const SoAFloatPrecision &distZIJ,
      const SoAFloatPrecision &distXJK, const SoAFloatPrecision &distYJK, const SoAFloatPrecision &distZJK,
      const SoAFloatPrecision &distXKI, const SoAFloatPrecision &distYKI, const SoAFloatPrecision &distZKI,
      const SoAFloatPrecision &distSquaredIJ, const SoAFloatPrecision &distSquaredJK,
      const SoAFloatPrecision &distSquaredKI, const SoAFloatPrecision &nu, SoAFloatPrecision &forceIX,
      SoAFloatPrecision &forceIY, SoAFloatPrecision &forceIZ, SoAFloatPrecision &factor,
      SoAFloatPrecision &allDotProducts, SoAFloatPrecision &allDistsSquared) const {
    // Calculate prefactor
    allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const SoAFloatPrecision allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    factor = 3.0 * nu / allDistsTo5;

    // Dot products
    const SoAFloatPrecision IJDotKI = distXIJ * distXKI + distYIJ * distYKI + distZIJ * distZKI;
    const SoAFloatPrecision IJDotJK = distXIJ * distXJK + distYIJ * distYJK + distZIJ * distZJK;
    const SoAFloatPrecision JKDotKI = distXJK * distXKI + distYJK * distYKI + distZJK * distZKI;
    allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    // Force I components
    const SoAFloatPrecision factorIDirectionJK = factor * IJDotKI * (IJDotJK - JKDotKI);
    const SoAFloatPrecision factorIDirectionIJ =
        factor * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI + 5.0 * allDotProducts / distSquaredIJ);
    const SoAFloatPrecision factorIDirectionKI =
        factor * (-IJDotJK * JKDotKI + distSquaredIJ * distSquaredJK - 5.0 * allDotProducts / distSquaredKI);

    forceIX = distXJK * factorIDirectionJK + distXIJ * factorIDirectionIJ + distXKI * factorIDirectionKI;
    forceIY = distYJK * factorIDirectionJK + distYIJ * factorIDirectionIJ + distYKI * factorIDirectionKI;
    forceIZ = distZJK * factorIDirectionJK + distZIJ * factorIDirectionIJ + distZKI * factorIDirectionKI;
  }

  /**
   * Inline helper to compute force components for particle I.
   * Returns tuple of (forceIX, forceIY, forceIZ, factor, allDotProducts, allDistsSquared)
   */
  __attribute__((always_inline)) inline auto SoAKernelN3(
      const SoAFloatPrecision &distXIJ, const SoAFloatPrecision &distYIJ, const SoAFloatPrecision &distZIJ,
      const SoAFloatPrecision &distXJK, const SoAFloatPrecision &distYJK, const SoAFloatPrecision &distZJK,
      const SoAFloatPrecision &distXKI, const SoAFloatPrecision &distYKI, const SoAFloatPrecision &distZKI,
      const SoAFloatPrecision &distSquaredIJ, const SoAFloatPrecision &distSquaredJK,
      const SoAFloatPrecision &distSquaredKI, const SoAFloatPrecision &nu, SoAFloatPrecision &forceIX,
      SoAFloatPrecision &forceIY, SoAFloatPrecision &forceIZ, SoAFloatPrecision &forceJX, SoAFloatPrecision &forceJY,
      SoAFloatPrecision &forceJZ, SoAFloatPrecision &factor, SoAFloatPrecision &allDotProducts,
      SoAFloatPrecision &allDistsSquared) const {
    // Calculate prefactor
    allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const SoAFloatPrecision allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    factor = 3.0 * nu / allDistsTo5;

    // Dot products
    const SoAFloatPrecision IJDotKI = distXIJ * distXKI + distYIJ * distYKI + distZIJ * distZKI;
    const SoAFloatPrecision IJDotJK = distXIJ * distXJK + distYIJ * distYJK + distZIJ * distZJK;
    const SoAFloatPrecision JKDotKI = distXJK * distXKI + distYJK * distYKI + distZJK * distZKI;
    allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    // Force I components
    const SoAFloatPrecision factorIDirectionJK = factor * IJDotKI * (IJDotJK - JKDotKI);
    const SoAFloatPrecision factorIDirectionIJ =
        factor * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI + 5.0 * allDotProducts / distSquaredIJ);
    const SoAFloatPrecision factorIDirectionKI =
        factor * (-IJDotJK * JKDotKI + distSquaredIJ * distSquaredJK - 5.0 * allDotProducts / distSquaredKI);

    forceIX = distXJK * factorIDirectionJK + distXIJ * factorIDirectionIJ + distXKI * factorIDirectionKI;
    forceIY = distYJK * factorIDirectionJK + distYIJ * factorIDirectionIJ + distYKI * factorIDirectionKI;
    forceIZ = distZJK * factorIDirectionJK + distZIJ * factorIDirectionIJ + distZKI * factorIDirectionKI;

    // Force J components
    const SoAFloatPrecision factorJDirectionKI = factor * IJDotJK * (JKDotKI - IJDotKI);
    const SoAFloatPrecision factorJDirectionIJ =
        factor * (-IJDotKI * JKDotKI + distSquaredJK * distSquaredKI - 5.0 * allDotProducts / distSquaredIJ);
    const SoAFloatPrecision factorJDirectionJK =
        factor * (IJDotKI * JKDotKI - distSquaredIJ * distSquaredKI + 5.0 * allDotProducts / distSquaredJK);

    forceJX = distXKI * factorJDirectionKI + distXIJ * factorJDirectionIJ + distXJK * factorJDirectionJK;
    forceJY = distYKI * factorJDirectionKI + distYIJ * factorJDirectionIJ + distYJK * factorJDirectionJK;
    forceJZ = distZKI * factorJDirectionKI + distZIJ * factorJDirectionIJ + distZJK * factorJDirectionJK;
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
      numTripletsCount = 0;
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
     * Number of totally traversed triplets counted.
     * Used for calculating the hit rate. Differs from the number of distance calculations because
     * not all 3 distances are always calculated.
     */
    size_t numTripletsCount = 0;

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
    double __remainingTo64[(64 - 6 * sizeof(size_t)) / sizeof(size_t)];
  };

  // make sure of the size of AoSThreadDataGlobals
  static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
  static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

  const double _cutoffSquared;

  // Parameter of the Axilrod-Teller-Muto potential
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

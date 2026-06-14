/**
 * @file AxilrodTellerMutoFunctor.h
 * @author M. Muehlhaeusser
 * @date 25/07/23
 */

#pragma once

#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

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
    const auto *const rawOwnedStatePtr = reinterpret_cast<const int64_t *>(ownedStatePtr);
    constexpr auto rawDummyState = static_cast<int64_t>(autopas::OwnershipState::dummy);
    constexpr auto rawOwnedState = static_cast<int64_t>(autopas::OwnershipState::owned);

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

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numGlobalCalcsSum = 0;

    const SoAFloatPrecision const_nu = _nu;
    const size_t soaSize = soa.size();

    // Store pairwise distances in a flattened triangular matrix as a look-up table.
    const size_t distArraySize = soaSize * (soaSize - 1) / 2;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> soaDistsX(distArraySize);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> soaDistsY(distArraySize);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> soaDistsZ(distArraySize);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> soaDistsR2(distArraySize);
    for (unsigned int i = 0; i < soaSize; ++i) {
      const size_t baseIndex = i * soaSize - (i * (i + 1) / 2) - i - 1;
      for (unsigned int j = i + 1; j < soaSize; ++j) {
        const SoAFloatPrecision distXIJ = xptr[j] - xptr[i];
        const SoAFloatPrecision distYIJ = yptr[j] - yptr[i];
        const SoAFloatPrecision distZIJ = zptr[j] - zptr[i];
        const SoAFloatPrecision distSquaredIJ = distXIJ * distXIJ + distYIJ * distYIJ + distZIJ * distZIJ;
        soaDistsX[baseIndex + j] = distXIJ;
        soaDistsY[baseIndex + j] = distYIJ;
        soaDistsZ[baseIndex + j] = distZIJ;
        soaDistsR2[baseIndex + j] = distSquaredIJ;
      }
    }

    for (unsigned int i = 0; i < soaSize - 2; ++i) {
      const auto ownedStateI = rawOwnedStatePtr[i];
      if (ownedStateI == rawDummyState) {
        continue;
      }
      SoAFloatPrecision fXAccI = 0.;
      SoAFloatPrecision fYAccI = 0.;
      SoAFloatPrecision fZAccI = 0.;

      for (unsigned int j = i + 1; j < soaSize - 1; ++j) {
        const auto ownedStateJ = rawOwnedStatePtr[j];
        if (ownedStateJ == rawDummyState) {
          continue;
        }

        // Get distance from the look-up table
        const size_t baseIndexI = i * soaSize - (i * (i + 1) / 2) - i - 1;
        const auto &distXIJ = soaDistsX[baseIndexI + j];
        const auto &distYIJ = soaDistsY[baseIndexI + j];
        const auto &distZIJ = soaDistsZ[baseIndexI + j];
        const auto &distSquaredIJ = soaDistsR2[baseIndexI + j];

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        SoAFloatPrecision fXAccJ = 0.;
        SoAFloatPrecision fYAccJ = 0.;
        SoAFloatPrecision fZAccJ = 0.;

#pragma omp simd reduction(+:fXAccI, fYAccI, fZAccI, fXAccJ, fYAccJ, fZAccJ, potentialEnergySum, virialSumX, virialSumY, virialSumZ) simdlen(8)
        for (unsigned int k = j + 1; k < soaSize; ++k) {
          const auto ownedStateK = rawOwnedStatePtr[k];
          // if (ownedStateK == rawDummyState) {
          //   continue;
          // }
          const auto &distXIK = soaDistsX[baseIndexI + k];
          const auto &distYIK = soaDistsY[baseIndexI + k];
          const auto &distZIK = soaDistsZ[baseIndexI + k];
          const auto &distSquaredIK = soaDistsR2[baseIndexI + k];

          // if (distSquaredIK > cutoffSquared) {
          //   continue;
          // }

          const size_t baseIndexJ = j * soaSize - (j * (j + 1) / 2) - j - 1;
          const auto &distXJK = soaDistsX[baseIndexJ + k];
          const auto &distYJK = soaDistsY[baseIndexJ + k];
          const auto &distZJK = soaDistsZ[baseIndexJ + k];
          const auto &distSquaredJK = soaDistsR2[baseIndexJ + k];

          // if (distSquaredJK > cutoffSquared) {
          //   continue;
          // }
          double mask =
              (ownedStateK != rawDummyState && distSquaredIK <= cutoffSquared && distSquaredJK <= cutoffSquared) ? 1.0
                                                                                                                 : 0.0;

          SoAFloatPrecision nu = const_nu;

          SoAFloatPrecision forceIX, forceIY, forceIZ;
          SoAFloatPrecision forceJX, forceJY, forceJZ;
          SoAFloatPrecision factor, allDotProducts, allDistsSquared;
          SoAKernelN3(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, -distXIK, -distYIK, -distZIK, distSquaredIJ,
                      distSquaredJK, distSquaredIK, nu, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ, factor,
                      allDotProducts, allDistsSquared);

          fXAccI += forceIX * mask;
          fYAccI += forceIY * mask;
          fZAccI += forceIZ * mask;

          fXAccJ += forceJX * mask;
          fYAccJ += forceJY * mask;
          fZAccJ += forceJZ * mask;

          const SoAFloatPrecision forceKX = -(forceIX + forceJX);
          const SoAFloatPrecision forceKY = -(forceIY + forceJY);
          const SoAFloatPrecision forceKZ = -(forceIZ + forceJZ);

          fxptr[k] += forceKX * mask;
          fyptr[k] += forceKY * mask;
          fzptr[k] += forceKZ * mask;
        }
        fxptr[j] += fXAccJ;
        fyptr[j] += fYAccJ;
        fzptr[j] += fZAccJ;
      }
      fxptr[i] += fXAccI;
      fyptr[i] += fYAccI;
      fzptr[i] += fZAccI;
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

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.size() == 0) return;

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
    const auto *const rawOwnedStatePtr1 = reinterpret_cast<const int64_t *>(ownedStatePtr1);
    const auto *const rawOwnedStatePtr2 = reinterpret_cast<const int64_t *>(ownedStatePtr2);
    constexpr auto rawDummyState = static_cast<int64_t>(autopas::OwnershipState::dummy);

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

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const SoAFloatPrecision const_nu = _nu;

    size_t soa1Size = soa1.size();
    size_t soa2Size = soa2.size();

    // thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsXJK;
    // thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsYJK;
    // thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsZJK;
    // thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsSquaredJK;
    // thread_local std::vector<size_t, autopas::AlignedAllocator<size_t>> indexJ;
    // thread_local std::vector<size_t, autopas::AlignedAllocator<size_t>> indexK;

    //// CASE 1: Particle i is the single particle in SoA1, j and k are in SoA2
    for (auto i = 0; i < soa1Size; ++i) {
      const auto ownedStateI = rawOwnedStatePtr1[i];
      if (ownedStateI == rawDummyState) {
        continue;
      }

      const SoAFloatPrecision xi = xptr1[i];
      const SoAFloatPrecision yi = yptr1[i];
      const SoAFloatPrecision zi = zptr1[i];

      SoAFloatPrecision fxAccI = 0.0;
      SoAFloatPrecision fyAccI = 0.0;
      SoAFloatPrecision fzAccI = 0.0;

      std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsXJK;
      std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsYJK;
      std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsZJK;
      std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsSquaredJK;
      std::vector<size_t, autopas::AlignedAllocator<size_t>> indexJ;
      distsXJK.reserve(soa2Size);
      distsYJK.reserve(soa2Size);
      distsZJK.reserve(soa2Size);
      distsSquaredJK.reserve(soa2Size);
      indexJ.reserve(soa2Size);

      for (auto j = 0; j < soa2Size; ++j) {
        const auto ownedStateJ = rawOwnedStatePtr2[j];
        if (ownedStateJ == rawDummyState) {
          continue;
        }

        const SoAFloatPrecision xj = xptr2[j];
        const SoAFloatPrecision yj = yptr2[j];
        const SoAFloatPrecision zj = zptr2[j];
        const auto distXIJ = xj - xi;
        const auto distYIJ = yj - yi;
        const auto distZIJ = zj - zi;
        const auto distSquaredIJ = distXIJ * distXIJ + distYIJ * distYIJ + distZIJ * distZIJ;

        if (distSquaredIJ <= cutoffSquared) {
          distsXJK.push_back(distXIJ);
          distsYJK.push_back(distYIJ);
          distsZJK.push_back(distZIJ);
          distsSquaredJK.push_back(distSquaredIJ);
          indexJ.push_back(j);
        }
      }
      for (auto j = 0; j < distsXJK.size(); ++j) {
        const auto distXIJ = distsXJK[j];
        const auto distYIJ = distsYJK[j];
        const auto distZIJ = distsZJK[j];
        const auto distSquaredIJ = distsSquaredJK[j];

        SoAFloatPrecision fxAccJ = 0.0;
        SoAFloatPrecision fyAccJ = 0.0;
        SoAFloatPrecision fzAccJ = 0.0;

#pragma omp simd reduction (+:fxAccI, fyAccI, fzAccI, fxAccJ, fyAccJ, fzAccJ, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (auto k = j + 1; k < distsXJK.size(); ++k) {
          const auto distXKI = -distsXJK[k];
          const auto distYKI = -distsYJK[k];
          const auto distZKI = -distsZJK[k];
          const auto distSquaredKI = distsSquaredJK[k];
          const auto idK = indexJ[k];

          SoAFloatPrecision nu = const_nu;

          const auto KIDotIJ = distXIJ * distXKI + distYIJ * distYKI + distZIJ * distZKI;
          const auto distSquaredJK = 2.0 * KIDotIJ + distSquaredIJ + distSquaredKI;

          const double mask = (distSquaredJK <= cutoffSquared) ? 1.0 : 0.0;

          if constexpr (not newton3) {
            // Compute all forces without Newton3
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision factor, allDotProducts, allDistsSquared;
            SoAKernelNoN3Refactored(distXIJ, distYIJ, distZIJ, distXKI, distYKI, distZKI, distSquaredIJ, distSquaredJK,
                                    distSquaredKI, KIDotIJ, nu, forceIX, forceIY, forceIZ, factor, allDotProducts,
                                    allDistsSquared);

            fxAccI += forceIX * mask;
            fyAccI += forceIY * mask;
            fzAccI += forceIZ * mask;

          } else {
            // Compute all forces with Newton3
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision forceJX, forceJY, forceJZ;
            SoAFloatPrecision factor, allDotProducts, allDistsSquared;
            SoAKernelN3Refactored(distXIJ, distYIJ, distZIJ, distXKI, distYKI, distZKI, distSquaredIJ, distSquaredJK,
                                  distSquaredKI, KIDotIJ, nu, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ,
                                  factor, allDotProducts, allDistsSquared);

            fxAccI += forceIX * mask;
            fyAccI += forceIY * mask;
            fzAccI += forceIZ * mask;

            fxAccJ += forceJX * mask;
            fyAccJ += forceJY * mask;
            fzAccJ += forceJZ * mask;

            const SoAFloatPrecision forceKX = forceIX + forceJX;
            const SoAFloatPrecision forceKY = forceIY + forceJY;
            const SoAFloatPrecision forceKZ = forceIZ + forceJZ;

            fxptr2[idK] -= forceKX * mask;
            fyptr2[idK] -= forceKY * mask;
            fzptr2[idK] -= forceKZ * mask;
          }
        }
        if constexpr (newton3) {
          fxptr2[j] += fxAccJ;
          fyptr2[j] += fyAccJ;
          fzptr2[j] += fzAccJ;
        }
      }
      fxptr1[i] += fxAccI;
      fyptr1[i] += fyAccI;
      fzptr1[i] += fzAccI;
    }

    //// CASE 2: Particle i is the single particle in SoA2, j and k are in SoA1
    //     for (auto i = 0; i < soa2Size; ++i) {
    //       const auto ownedStateI = ownedStatePtr2[i];
    //       if (ownedStateI == autopas::OwnershipState::dummy) {
    //         continue;
    //       }
    //       const auto isOwnedI = ownedStateI == autopas::OwnershipState::owned ? 1.0 : 0.0;
    //
    //       const SoAFloatPrecision xi = xptr2[i];
    //       const SoAFloatPrecision yi = yptr2[i];
    //       const SoAFloatPrecision zi = zptr2[i];
    //
    //       SoAFloatPrecision fxAccI = 0.0;
    //       SoAFloatPrecision fyAccI = 0.0;
    //       SoAFloatPrecision fzAccI = 0.0;
    //
    //       // Step 2: Iterate over unique pairs in the neighbor list
    //       for (auto k = 0; k < soa1Size; ++k) {
    //         const auto ownedStateK = ownedStatePtr1[k];
    //         if (ownedStateK == autopas::OwnershipState::dummy) {
    //           continue;
    //         }
    //         const auto isOwnedK = ownedStateK == autopas::OwnershipState::owned ? 1.0 : 0.0;
    //
    //         const auto xk = xptr1[k];
    //         const auto yk = yptr1[k];
    //         const auto zk = zptr1[k];
    //
    //         const auto distXKI = xi - xk;
    //         const auto distYKI = yi - yk;
    //         const auto distZKI = zi - zk;
    //         const auto distSquaredKI = distXKI * distXKI + distYKI * distYKI + distZKI * distZKI;
    //
    //         if (distSquaredKI > cutoffSquared) {
    //           continue;
    //         }
    //
    //         SoAFloatPrecision fxAccK = 0.0;
    //         SoAFloatPrecision fyAccK = 0.0;
    //         SoAFloatPrecision fzAccK = 0.0;
    //
    // #pragma omp simd reduction (+:fxAccK, fyAccK, fzAccK, fxAccI, fyAccI, fzAccI, potentialEnergySum, virialSumX,
    // virialSumY, virialSumZ)
    //         for (auto j = k + 1; j < soa1Size; ++j) {
    //           const auto ownedStateJ = ownedStatePtr1[j];
    //           if (ownedStateJ == autopas::OwnershipState::dummy) {
    //             continue;
    //           }
    //           const auto isOwnedJ = ownedStateJ == autopas::OwnershipState::owned ? 1.0 : 0.0;
    //
    //           const auto xj = xptr1[j];
    //           const auto yj = yptr1[j];
    //           const auto zj = zptr1[j];
    //
    //           // Compute j-k distance
    //           const SoAFloatPrecision distXJK = xk - xj;
    //           const SoAFloatPrecision distYJK = yk - yj;
    //           const SoAFloatPrecision distZJK = zk - zj;
    //           const SoAFloatPrecision distSquaredJK = distXJK * distXJK + distYJK * distYJK + distZJK * distZJK;
    //
    //           if (distSquaredJK > cutoffSquared) {
    //             continue;
    //           }
    //
    //           const auto JKDotKI = distXJK * distXKI + distYJK * distYKI + distZJK * distZKI;
    //           const auto distSquaredIJ = 2.0 * JKDotKI + distSquaredJK + distSquaredKI;
    //
    //           if (distSquaredIJ > cutoffSquared) {
    //             continue;
    //           }
    //
    //           SoAFloatPrecision nu = const_nu;
    //
    //           // Compute all forces with Newton3
    //           SoAFloatPrecision forceJX, forceJY, forceJZ;
    //           SoAFloatPrecision forceKX, forceKY, forceKZ;
    //           SoAFloatPrecision factor, allDotProducts, allDistsSquared;
    //           SoAKernelN3Refactored(distXJK, distYJK, distZJK, distXKI, distYKI, distZKI, distSquaredJK,
    //           distSquaredKI,
    //                                 distSquaredIJ, JKDotKI, nu, forceJX, forceJY, forceJZ, forceKX, forceKY, forceKZ,
    //                                 factor, allDotProducts, allDistsSquared);
    //
    //           fxptr1[j] += forceJX;
    //           fyptr1[j] += forceJY;
    //           fzptr1[j] += forceJZ;
    //
    //           fxAccK += forceKX;
    //           fyAccK += forceKY;
    //           fzAccK += forceKZ;
    //
    //           SoAFloatPrecision localPotentialEnergy = 0.0;
    //           SoAFloatPrecision localVirialX = 0.0;
    //           SoAFloatPrecision localVirialY = 0.0;
    //           SoAFloatPrecision localVirialZ = 0.0;
    //
    //           if constexpr (newton3) {
    //             const SoAFloatPrecision forceIX = forceKX + forceJX;
    //             const SoAFloatPrecision forceIY = forceKY + forceJY;
    //             const SoAFloatPrecision forceIZ = forceKZ + forceJZ;
    //             fxAccI -= forceIX;
    //             fyAccI -= forceIY;
    //             fzAccI -= forceIZ;
    //           }
    //         }
    //         fxptr1[k] += fxAccK;
    //         fyptr1[k] += fyAccK;
    //         fzptr1[k] += fzAccK;
    //       }
    //       if constexpr (newton3) {
    //         fxptr2[i] += fxAccI;
    //         fyptr2[i] += fyAccI;
    //         fzptr2[i] += fzAccI;
    //       }
    //     }
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
    const auto *const rawOwnedStatePtr1 = reinterpret_cast<const int64_t *>(ownedStatePtr1);
    const auto *const rawOwnedStatePtr2 = reinterpret_cast<const int64_t *>(ownedStatePtr2);
    const auto *const rawOwnedStatePtr3 = reinterpret_cast<const int64_t *>(ownedStatePtr3);
    constexpr auto rawDummyState = static_cast<int64_t>(autopas::OwnershipState::dummy);

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

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const SoAFloatPrecision const_nu = _nu;

    const auto soa1Size = soa1.size();
    const auto soa2Size = soa2.size();
    const auto soa3Size = soa3.size();

    // tread_local to keep memory for the whole traversal
    thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsXJK(8);
    thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsYJK(8);
    thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsZJK(8);
    thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distsSquaredJK(8);
    thread_local std::vector<size_t, autopas::AlignedAllocator<size_t>> indexJ(8);
    thread_local std::vector<size_t, autopas::AlignedAllocator<size_t>> indexK(8);

    // Use a manual counter instead of clear()
    size_t numPairs = 0;

    const size_t jkSizeEstimate = soa2Size * soa3Size / 4;  // Heuristic estimate

    // Note: Using resize() so we can use operator[] safely
    if (distsXJK.size() < jkSizeEstimate) {
      distsXJK.resize(jkSizeEstimate);
      distsYJK.resize(jkSizeEstimate);
      distsZJK.resize(jkSizeEstimate);
      distsSquaredJK.resize(jkSizeEstimate);
      indexJ.resize(jkSizeEstimate);
      indexK.resize(jkSizeEstimate);
    }

    for (uint32_t j = 0; j < soa2Size; ++j) {
      if (rawOwnedStatePtr2[j] == rawDummyState) {
        continue;
      }

      // Prefetch xj, yj, zj to avoid re-reading them inside the k-loop
      const SoAFloatPrecision xj = xptr2[j];
      const SoAFloatPrecision yj = yptr2[j];
      const SoAFloatPrecision zj = zptr2[j];

      for (uint32_t k = 0; k < soa3Size; ++k) {
        if (rawOwnedStatePtr3[k] == rawDummyState) {
          continue;
        }

        const SoAFloatPrecision distXJK = xptr3[k] - xj;
        const SoAFloatPrecision distYJK = yptr3[k] - yj;
        const SoAFloatPrecision distZJK = zptr3[k] - zj;
        const SoAFloatPrecision distSquaredJK = distXJK * distXJK + distYJK * distYJK + distZJK * distZJK;

        if (distSquaredJK <= cutoffSquared) {
          // Manual capacity check (rarely triggers, highly predictable)
          if (numPairs == distsXJK.size()) {
            const size_t newSize = distsXJK.size() * 2;  // Double the size
            distsXJK.resize(newSize);
            distsYJK.resize(newSize);
            distsZJK.resize(newSize);
            distsSquaredJK.resize(newSize);
            indexJ.resize(newSize);
            indexK.resize(newSize);
          }

          // Direct array indexing (blazing fast)
          distsXJK[numPairs] = distXJK;
          distsYJK[numPairs] = distYJK;
          distsZJK[numPairs] = distZJK;
          distsSquaredJK[numPairs] = distSquaredJK;
          indexJ[numPairs] = j;
          indexK[numPairs] = k;
          numPairs++;
        }
      }
    }

    // Loop over all valid j-k pairs
    for (size_t jk = 0; jk < numPairs; ++jk) {
      const auto idJ = indexJ[jk];
      const auto idK = indexK[jk];
      const SoAFloatPrecision distXJK = distsXJK[jk];
      const SoAFloatPrecision distYJK = distsYJK[jk];
      const SoAFloatPrecision distZJK = distsZJK[jk];
      const SoAFloatPrecision distSquaredJK = distsSquaredJK[jk];
      const auto xj = xptr2[idJ];
      const auto yj = yptr2[idJ];
      const auto zj = zptr2[idJ];

      SoAFloatPrecision fXAccJ = 0.;
      SoAFloatPrecision fYAccJ = 0.;
      SoAFloatPrecision fZAccJ = 0.;
      SoAFloatPrecision fXAccK = 0.;
      SoAFloatPrecision fYAccK = 0.;
      SoAFloatPrecision fZAccK = 0.;

#pragma omp simd reduction(+:fXAccK, fYAccK, fZAccK, fXAccJ, fYAccJ, fZAccJ, potentialEnergySum, virialSumX, virialSumY, virialSumZ) simdlen(8)
      for (size_t i = 0; i < soa1Size; ++i) {
        const SoAFloatPrecision distXIJ = xj - xptr1[i];
        const SoAFloatPrecision distYIJ = yj - yptr1[i];
        const SoAFloatPrecision distZIJ = zj - zptr1[i];
        const SoAFloatPrecision distSquaredIJ = distXIJ * distXIJ + distYIJ * distYIJ + distZIJ * distZIJ;

        const SoAFloatPrecision IJDotJK = distXIJ * distXJK + distYIJ * distYJK + distZIJ * distZJK;
        const SoAFloatPrecision distSquaredKI = 2.0 * IJDotJK + distSquaredIJ + distSquaredJK;

        const double mask =
            (rawOwnedStatePtr1[i] != rawDummyState && distSquaredIJ <= cutoffSquared && distSquaredKI <= cutoffSquared)
                ? 1.0
                : 0.0;

        SoAFloatPrecision forceIX, forceIY, forceIZ;
        SoAFloatPrecision factor, allDotProducts, allDistsSquared;
        if constexpr (not newton3) {
          SoAKernelNoN3Refactored(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, distSquaredIJ, distSquaredJK,
                                  distSquaredKI, IJDotJK, const_nu, forceIX, forceIY, forceIZ, factor, allDotProducts,
                                  allDistsSquared);
          fxptr1[i] += forceIX * mask;
          fyptr1[i] += forceIY * mask;
          fzptr1[i] += forceIZ * mask;
        } else {
          SoAFloatPrecision forceJX, forceJY, forceJZ;
          SoAKernelN3Refactored(distXIJ, distYIJ, distZIJ, distXJK, distYJK, distZJK, distSquaredIJ, distSquaredJK,
                                distSquaredKI, IJDotJK, const_nu, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ,
                                factor, allDotProducts, allDistsSquared);
          fxptr1[i] += forceIX * mask;
          fyptr1[i] += forceIY * mask;
          fzptr1[i] += forceIZ * mask;

          fXAccJ += forceJX * mask;
          fYAccJ += forceJY * mask;
          fZAccJ += forceJZ * mask;

          fXAccK -= (forceIX + forceJX) * mask;
          fYAccK -= (forceIY + forceJY) * mask;
          fZAccK -= (forceIZ + forceJZ) * mask;
        }
      }
      if constexpr (newton3) {
        fxptr2[idJ] += fXAccJ;
        fyptr2[idJ] += fYAccJ;
        fzptr2[idJ] += fZAccJ;
        fxptr3[idK] += fXAccK;
        fyptr3[idK] += fYAccK;
        fzptr3[idK] += fZAccK;
      }
    }
  }

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    const auto threadNum = autopas::autopas_get_thread_num();
    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typePtr = soa.template begin<Particle_T::AttributeNames::typeId>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;  // Note: This is not the potential energy but some fixed multiple of it.
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    SoAFloatPrecision fxAccI = 0;
    SoAFloatPrecision fyAccI = 0;
    SoAFloatPrecision fzAccI = 0;
    auto i = indexFirst;

    const SoAFloatPrecision const_nu = _nu;
    const size_t soaSize = soa.size();

    const auto ownedStateI = ownedStatePtr[i];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }
    const auto isOwnedI = ownedStatePtr[i] == autopas::OwnershipState::owned ? 1.0 : 0.0;

    const auto numNeighbors = neighborList.size();
    std::vector<int, autopas::AlignedAllocator<int>> idNeighbors{};
    idNeighbors.reserve(numNeighbors);
    std::vector<int, autopas::AlignedAllocator<int>> typeNeighbors{};
    typeNeighbors.reserve(numNeighbors);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distXNeighbors{};
    distXNeighbors.reserve(numNeighbors);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distYNeighbors{};
    distYNeighbors.reserve(numNeighbors);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distZNeighbors{};
    distZNeighbors.reserve(numNeighbors);
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> distSquaredNeighbors{};
    distSquaredNeighbors.reserve(numNeighbors);
    std::vector<autopas::OwnershipState, autopas::AlignedAllocator<autopas::OwnershipState>> ownershipNeighbors{};
    ownershipNeighbors.reserve(numNeighbors);

    for (size_t neighborIdx = 0; neighborIdx < neighborList.size(); ++neighborIdx) {
      const auto j = neighborList[neighborIdx];
      const auto distX = xPtr[j] - xPtr[i];
      const auto distY = yPtr[j] - yPtr[i];
      const auto distZ = zPtr[j] - zPtr[i];
      const auto distSquared = distX * distX + distY * distY + distZ * distZ;
      if (distSquared <= cutoffSquared and ownedStatePtr[j] != autopas::OwnershipState::dummy and i != j) {
        idNeighbors.push_back(j);
        typeNeighbors.push_back(typePtr[j]);
        distXNeighbors.push_back(distX);
        distYNeighbors.push_back(distY);
        distZNeighbors.push_back(distZ);
        distSquaredNeighbors.push_back(distSquared);
        ownershipNeighbors.push_back(ownedStatePtr[j]);
      }
    }
    SoAFloatPrecision fxAccJ = 0;
    SoAFloatPrecision fyAccJ = 0;
    SoAFloatPrecision fzAccJ = 0;

    const auto actualNeighborSize = idNeighbors.size();
    for (size_t j = 0; j < actualNeighborSize; ++j) {
      const auto distXIJ = distXNeighbors[j];
      const auto distYIJ = distYNeighbors[j];
      const auto distZIJ = distZNeighbors[j];
      const auto distSquaredIJ = distSquaredNeighbors[j];
      const auto isOwnedJ = ownedStatePtr[j] == autopas::OwnershipState::owned ? 1.0 : 0.0;

#pragma omp simd reduction(+ : fxAccI, fyAccI, fzAccI, fxAccJ, fyAccJ, fzAccJ, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
      for (size_t k = j + 1; k < actualNeighborSize; ++k) {
        const auto isOwnedK = ownedStatePtr[k] == autopas::OwnershipState::owned ? 1.0 : 0.0;

        const auto distXKI = -distXNeighbors[k];
        const auto distYKI = -distYNeighbors[k];
        const auto distZKI = -distZNeighbors[k];
        const auto distSquaredKI = distSquaredNeighbors[k];

        const auto IJDotKI = distXKI * distXIJ + distYKI * distYIJ + distZKI * distZIJ;
        const auto distSquaredJK = 2.0 * IJDotKI + distSquaredIJ + distSquaredKI;

        if (distSquaredJK > cutoffSquared) {
          continue;
        }

        SoAFloatPrecision nu = const_nu;

        if constexpr (not newton3) {
          SoAFloatPrecision forceIX, forceIY, forceIZ;
          SoAFloatPrecision factor, allDotProducts, allDistsSquared;
          SoAKernelNoN3Verlet(distXKI, distYKI, distZKI, distXIJ, distYIJ, distZIJ, distSquaredKI, distSquaredIJ,
                              distSquaredJK, IJDotKI, nu, forceIX, forceIY, forceIZ, factor, allDotProducts,
                              allDistsSquared);
          fxAccI += forceIX;
          fyAccI += forceIY;
          fzAccI += forceIZ;
        } else {
          SoAFloatPrecision forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ;
          SoAFloatPrecision factor, allDotProducts, allDistsSquared;
          SoAKernelN3Verlet(distXKI, distYKI, distZKI, distXIJ, distYIJ, distZIJ, distSquaredKI, distSquaredIJ,
                            distSquaredJK, IJDotKI, nu, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ, factor,
                            allDotProducts, allDistsSquared);
          fxAccI += forceIX;
          fyAccI += forceIY;
          fzAccI += forceIZ;

          fxAccJ += forceJX;
          fyAccJ += forceJY;
          fzAccJ += forceJZ;

          const auto forceKX = -(forceIX + forceJX);
          const auto forceKY = -(forceIY + forceJY);
          const auto forceKZ = -(forceIZ + forceJZ);

          fxPtr[k] += forceKX;
          fyPtr[k] += forceKY;
          fzPtr[k] += forceKZ;
        }
      }
      if constexpr (newton3) {
        fxPtr[j] += fxAccJ;
        fyPtr[j] += fyAccJ;
        fzPtr[j] += fzAccJ;
        fxAccJ = 0.0;
        fyAccJ = 0.0;
        fzAccJ = 0.0;
      }
    }

    fxPtr[i] += fxAccI;
    fyPtr[i] += fyAccI;
    fzPtr[i] += fzAccI;
  }

  __attribute__((always_inline)) inline void SoAKernelNoN3Refactored(
      const SoAFloatPrecision &distXIJ, const SoAFloatPrecision &distYIJ, const SoAFloatPrecision &distZIJ,
      const SoAFloatPrecision &distXJK, const SoAFloatPrecision &distYJK, const SoAFloatPrecision &distZJK,
      const SoAFloatPrecision &distSquaredIJ, const SoAFloatPrecision &distSquaredJK,
      const SoAFloatPrecision &distSquaredKI, const SoAFloatPrecision &IJDotJK, const SoAFloatPrecision &nu,
      SoAFloatPrecision &forceIX, SoAFloatPrecision &forceIY, SoAFloatPrecision &forceIZ, SoAFloatPrecision &factor,
      SoAFloatPrecision &allDotProducts, SoAFloatPrecision &allDistsSquared) const {
    // Calculate prefactor
    allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;

    const SoAFloatPrecision allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    factor = 3.0 * nu / allDistsTo5;

    // Dot products
    const SoAFloatPrecision IJDotKI = -(distSquaredIJ + IJDotJK);
    const SoAFloatPrecision JKDotKI = -(IJDotJK + distSquaredJK);
    const SoAFloatPrecision IJdotJKdotKI = IJDotJK * JKDotKI;
    allDotProducts = IJDotKI * IJdotJKdotKI;

    // Force I components
    const SoAFloatPrecision factorIDirectionJK = IJDotKI * (IJDotJK - JKDotKI);
    const SoAFloatPrecision factorIDirectionIJ =
        (-IJdotJKdotKI * distSquaredIJ + allDistsSquared - 5.0 * allDotProducts) / distSquaredIJ;
    const SoAFloatPrecision factorIDirectionKI =
        (IJdotJKdotKI * distSquaredKI - allDistsSquared + 5.0 * allDotProducts) / distSquaredKI;

    const SoAFloatPrecision finalFactorIJ = factor * (factorIDirectionIJ - factorIDirectionKI);
    const SoAFloatPrecision finalFactorJK = factor * (factorIDirectionJK + factorIDirectionKI);

    forceIX = distXJK * finalFactorJK - distXIJ * finalFactorIJ;
    forceIY = distYJK * finalFactorJK - distYIJ * finalFactorIJ;
    forceIZ = distZJK * finalFactorJK - distZIJ * finalFactorIJ;
  }

  __attribute__((always_inline)) inline void SoAKernelN3Refactored(
      const SoAFloatPrecision &distXIJ, const SoAFloatPrecision &distYIJ, const SoAFloatPrecision &distZIJ,
      const SoAFloatPrecision &distXJK, const SoAFloatPrecision &distYJK, const SoAFloatPrecision &distZJK,
      const SoAFloatPrecision &distSquaredIJ, const SoAFloatPrecision &distSquaredJK,
      const SoAFloatPrecision &distSquaredKI, const SoAFloatPrecision &IJDotJK, const SoAFloatPrecision &nu,
      SoAFloatPrecision &forceIX, SoAFloatPrecision &forceIY, SoAFloatPrecision &forceIZ, SoAFloatPrecision &forceJX,
      SoAFloatPrecision &forceJY, SoAFloatPrecision &forceJZ, SoAFloatPrecision &factor,
      SoAFloatPrecision &allDotProducts, SoAFloatPrecision &allDistsSquared) const {
    // Calculate prefactor
    allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;

    const SoAFloatPrecision allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    factor = 3.0 * nu / allDistsTo5;

    // Dot products
    const SoAFloatPrecision IJDotKI = -(distSquaredIJ + IJDotJK);
    const SoAFloatPrecision JKDotKI = -(IJDotJK + distSquaredJK);
    const SoAFloatPrecision IJdotJKdotKI = IJDotJK * JKDotKI;
    allDotProducts = IJDotKI * IJdotJKdotKI;
    const SoAFloatPrecision allDots5 = allDotProducts * 5.0;

    // Force I components
    const SoAFloatPrecision factorIDirectionJK = IJDotKI * (IJDotJK - JKDotKI);
    const SoAFloatPrecision factorIDirectionIJ =
        (-IJdotJKdotKI * distSquaredIJ + allDistsSquared - allDots5) / distSquaredIJ;
    const SoAFloatPrecision factorIDirectionKI =
        (IJdotJKdotKI * distSquaredKI - allDistsSquared + allDots5) / distSquaredKI;

    const SoAFloatPrecision finalFactorIJ = factor * (factorIDirectionIJ - factorIDirectionKI);
    const SoAFloatPrecision finalFactorJK = factor * (factorIDirectionJK + factorIDirectionKI);

    forceIX = distXJK * finalFactorJK - distXIJ * finalFactorIJ;
    forceIY = distYJK * finalFactorJK - distYIJ * finalFactorIJ;
    forceIZ = distZJK * finalFactorJK - distZIJ * finalFactorIJ;

    const SoAFloatPrecision factorJDirectionKI = IJDotJK * (JKDotKI - IJDotKI);
    const SoAFloatPrecision factorJDirectionIJ =
        (-IJDotKI * JKDotKI * distSquaredIJ + allDistsSquared - allDots5) / distSquaredIJ;
    const SoAFloatPrecision factorJDirectionJK =
        (IJDotKI * JKDotKI * distSquaredJK - allDistsSquared + allDots5) / distSquaredJK;

    const SoAFloatPrecision finalFactorIJOnJ = factor * (factorJDirectionIJ - factorJDirectionKI);
    const SoAFloatPrecision finalFactorJKOnJ = factor * (factorJDirectionJK - factorJDirectionKI);

    forceJX = distXJK * finalFactorJKOnJ + distXIJ * finalFactorIJOnJ;
    forceJY = distYJK * finalFactorJKOnJ + distYIJ * finalFactorIJOnJ;
    forceJZ = distZJK * finalFactorJKOnJ + distZIJ * finalFactorIJOnJ;
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

  __attribute__((always_inline)) inline void SoAKernelNoN3Verlet(
      const SoAFloatPrecision &distXKI, const SoAFloatPrecision &distYKI, const SoAFloatPrecision &distZKI,
      const SoAFloatPrecision &distXIJ, const SoAFloatPrecision &distYIJ, const SoAFloatPrecision &distZIJ,
      const SoAFloatPrecision &distSquaredKI, const SoAFloatPrecision &distSquaredIJ,
      const SoAFloatPrecision &distSquaredJK, const SoAFloatPrecision &IJDotKI, const SoAFloatPrecision &nu,
      SoAFloatPrecision &forceIX, SoAFloatPrecision &forceIY, SoAFloatPrecision &forceIZ, SoAFloatPrecision &factor,
      SoAFloatPrecision &allDotProducts, SoAFloatPrecision &allDistsSquared) const {
    // Calculate prefactor
    allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const SoAFloatPrecision allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    factor = 3.0 * nu / allDistsTo5;

    // Dot products using scalar relationships
    const SoAFloatPrecision IJDotJK = -distSquaredIJ - IJDotKI;
    const SoAFloatPrecision JKDotKI = -distSquaredKI - IJDotKI;
    allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    // Force I components
    const SoAFloatPrecision factorIDirectionJK = factor * IJDotKI * (IJDotJK - JKDotKI);
    const SoAFloatPrecision factorIDirectionIJ =
        factor * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI + 5.0 * allDotProducts / distSquaredIJ);
    const SoAFloatPrecision factorIDirectionKI =
        factor * (-IJDotJK * JKDotKI + distSquaredIJ * distSquaredJK - 5.0 * allDotProducts / distSquaredKI);

    // Compute coordinate components of JK vector
    const SoAFloatPrecision distXJK = -(distXIJ + distXKI);
    const SoAFloatPrecision distYJK = -(distYIJ + distYKI);
    const SoAFloatPrecision distZJK = -(distZIJ + distZKI);

    forceIX = distXJK * factorIDirectionJK + distXIJ * factorIDirectionIJ + distXKI * factorIDirectionKI;
    forceIY = distYJK * factorIDirectionJK + distYIJ * factorIDirectionIJ + distYKI * factorIDirectionKI;
    forceIZ = distZJK * factorIDirectionJK + distZIJ * factorIDirectionIJ + distZKI * factorIDirectionKI;
  }

  /**
   * Inline helper to compute force components for particles I and J (Newton 3 enabled) in a Verlet list.
   */
  __attribute__((always_inline)) inline void SoAKernelN3Verlet(
      const SoAFloatPrecision &distXKI, const SoAFloatPrecision &distYKI, const SoAFloatPrecision &distZKI,
      const SoAFloatPrecision &distXIJ, const SoAFloatPrecision &distYIJ, const SoAFloatPrecision &distZIJ,
      const SoAFloatPrecision &distSquaredKI, const SoAFloatPrecision &distSquaredIJ,
      const SoAFloatPrecision &distSquaredJK, const SoAFloatPrecision &IJDotKI, const SoAFloatPrecision &nu,
      SoAFloatPrecision &forceIX, SoAFloatPrecision &forceIY, SoAFloatPrecision &forceIZ, SoAFloatPrecision &forceJX,
      SoAFloatPrecision &forceJY, SoAFloatPrecision &forceJZ, SoAFloatPrecision &factor,
      SoAFloatPrecision &allDotProducts, SoAFloatPrecision &allDistsSquared) const {
    // Calculate prefactor
    allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const SoAFloatPrecision allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
    factor = 3.0 * nu / allDistsTo5;

    // Dot products using scalar relationships
    const SoAFloatPrecision IJDotJK = -distSquaredIJ - IJDotKI;
    const SoAFloatPrecision JKDotKI = -distSquaredKI - IJDotKI;
    allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    // Force I components
    const SoAFloatPrecision factorIDirectionJK = factor * IJDotKI * (IJDotJK - JKDotKI);
    const SoAFloatPrecision factorIDirectionIJ =
        factor * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI + 5.0 * allDotProducts / distSquaredIJ);
    const SoAFloatPrecision factorIDirectionKI =
        factor * (-IJDotJK * JKDotKI + distSquaredIJ * distSquaredJK - 5.0 * allDotProducts / distSquaredKI);

    // Compute coordinate components of JK vector
    const SoAFloatPrecision distXJK = -(distXIJ + distXKI);
    const SoAFloatPrecision distYJK = -(distYIJ + distYKI);
    const SoAFloatPrecision distZJK = -(distZIJ + distZKI);

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

/**
 * @file KryptonExtendedATMFunctor.h
 * @author muehlhaeusser
 * @date 08.08.2024
 */

#pragma once

#include <array>

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
 * A functor to handle the interactions between three krypton atoms using the extended AxilrodTeller potential.
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * The functor follows the potential described by Jäger et al. in: https://doi.org/10.1063/1.4943959
 *
 * @note  All calculations assume units to be in Angström and Kelvin.
 *
 * @tparam Particle_T The type of particle.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 */
template <class Particle_T, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool countFLOPs = false>
class KryptonExtendedATMFunctor
    : public autopas::TriwiseFunctor<Particle_T,
                                     KryptonExtendedATMFunctor<Particle_T, useNewton3, calculateGlobals, countFLOPs>> {
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
  KryptonExtendedATMFunctor() = delete;

  /**
   * Constructor
   * @param cutoff
   */
  explicit KryptonExtendedATMFunctor(double cutoff)
      : autopas::TriwiseFunctor<Particle_T,
                                KryptonExtendedATMFunctor<Particle_T, useNewton3, calculateGlobals, countFLOPs>>(
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

  std::string getName() final { return "KryptonExtendedATMFunctor"; }

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
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
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

    // Actual distances
    const double distIJ = std::sqrt(distSquaredIJ);
    const double distJK = std::sqrt(distSquaredJK);
    const double distKI = std::sqrt(distSquaredKI);

    // Numerators of cosine representation (cos_i = (r_ij^2 + r_ik^2 - r_jk^2) / (2 * r_ij * r_ik)
    const double numKI = distSquaredIJ + distSquaredJK - distSquaredKI;
    const double numJK = distSquaredIJ + distSquaredKI - distSquaredJK;
    const double numIJ = distSquaredJK + distSquaredKI - distSquaredIJ;

    const double numerator = numKI * numJK * numIJ;

    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDists = distIJ * distJK * distKI;
    const double allDistsTripled = allDistsSquared * allDists;

    // Gradient factors of 1. / (rrr)^3
    const double allDistsTriplesGradientIJ = 3. / (allDistsTripled * distSquaredIJ);
    const double allDistsTriplesGradientKI = -3. / (allDistsTripled * distSquaredKI);

    // Product of all cosines multiplied with 3: 3 * cos(a)cos(b)cos(c)
    const double cosines = (3. / 8.) * numerator / allDistsSquared;
    const double cosinesGradientIJ =
        (3. / 4.) * ((numerator / distSquaredIJ - numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);
    const double cosinesGradientKI =
        (3. / 4.) * ((-numerator / distSquaredKI + numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);

    // Gradient factors corresponding to the normal ATM term
    const auto fullATMGradientIJ =
        _nu * ((1. + cosines) * allDistsTriplesGradientIJ + cosinesGradientIJ / allDistsTripled);
    const auto fullATMGradientKI =
        _nu * ((1. + cosines) * allDistsTriplesGradientKI + cosinesGradientKI / allDistsTripled);

    const double expTerm = std::exp(-_alpha * (distIJ + distJK + distKI));

    // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
    std::array<double, 6> sumFactors{};
    double sum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      sumFactors[n] = _constantsA[n] * std::pow((distIJ * distJK * distKI), 2. * n / 3.);
      sum += sumFactors[n];
    }

    // Gradient factor of the sum in ij-direction
    double ijSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - _alpha);
    }

    // Gradient factor of the sum in ki-direction
    double kiSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      kiSum += sumFactors[n] * (2. * n / (3. * distKI) - _alpha);
    }

    // Total gradient factors for the exponential term times the cosines term
    const double fullExpGradientIJ = expTerm * (-(1. + cosines) * ijSum / distIJ + cosinesGradientIJ * sum);
    const double fullExpGradientKI = expTerm * ((1. + cosines) * kiSum / distKI + cosinesGradientKI * sum);

    // Assembling the forces
    const auto forceIDirectionIJ = displacementIJ * (fullATMGradientIJ + fullExpGradientIJ);
    const auto forceIDirecationKI = displacementKI * (fullATMGradientKI + fullExpGradientKI);

    const auto forceI = (forceIDirectionIJ + forceIDirecationKI) * (-1.0);

    i.addF(forceI);

    auto forceJ = forceI;
    auto forceK = forceI;
    if (newton3) {
      // Calculate all components for jk-direction
      const double allDistsTriplesGradientJK = 3. / (allDistsTripled * distSquaredJK);
      const double cosinesGradientJK =
          (3. / 4.) * ((numerator / distSquaredJK + numKI * numIJ - numJK * numIJ - numJK * numKI) / allDistsSquared);
      const auto fullATMGradientJK =
          _nu * ((1. + cosines) * allDistsTriplesGradientJK + cosinesGradientJK / allDistsTripled);

      double jkSum = 0.0;
      for (auto n = 0; n < sumFactors.size(); n++) {
        jkSum += sumFactors[n] * (2. * n / (3. * distJK) - _alpha);
      }
      const double fullExpGradientJK = expTerm * (-(1. + cosines) * jkSum / distJK + cosinesGradientJK * sum);

      // Assembling the forces
      const auto forceJDirectionIJ = displacementIJ * (-fullATMGradientIJ - fullExpGradientIJ);
      const auto forceJDirectionJK = displacementJK * (fullATMGradientJK + fullExpGradientJK);

      forceJ = (forceJDirectionIJ + forceJDirectionJK) * (-1.0);
      j.addF(forceJ);

      // Using newton's third law for the force on particle k
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
      const double potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

      // Virial is calculated as f_i * r_i
      // see Thompson et al.: https://doi.org/10.1063/1.3245303
      const auto virialI = forceI * i.getR();
      if (i.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadDataGlobals[threadnum].virialSum += virialI;
      }
      // for non-newton3 particles j and/or k will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        const auto virialJ = forceJ * j.getR();
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadDataGlobals[threadnum].virialSum += virialJ;
      }
      if (newton3 and k.isOwned()) {
        const auto virialK = forceK * k.getR();
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
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
    const SoAFloatPrecision const_alpha = _alpha;

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

        const SoAFloatPrecision displacementXIJ = xj - xi;
        const SoAFloatPrecision displacementYIJ = yj - yi;
        const SoAFloatPrecision displacementZIJ = zj - zi;
        const SoAFloatPrecision distSquaredIJ =
            displacementXIJ * displacementXIJ + displacementYIJ * displacementYIJ + displacementZIJ * displacementZIJ;

        if constexpr (countFLOPs) {
          ++numDistanceCalculationSum;
        }
        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        const SoAFloatPrecision distIJ = std::sqrt(distSquaredIJ);

        SoAFloatPrecision fXAccJ = 0.;
        SoAFloatPrecision fYAccJ = 0.;
        SoAFloatPrecision fZAccJ = 0.;

        for (unsigned int k = j + 1; k < soaSize; ++k) {
          const auto ownedStateK = ownedStatePtr[k];
          if (ownedStateK == autopas::OwnershipState::dummy) {
            continue;
          }

          const SoAFloatPrecision xk = xptr[k];
          const SoAFloatPrecision yk = yptr[k];
          const SoAFloatPrecision zk = zptr[k];

          const SoAFloatPrecision displacementXJK = xk - xj;
          const SoAFloatPrecision displacementYJK = yk - yj;
          const SoAFloatPrecision displacementZJK = zk - zj;
          const SoAFloatPrecision distSquaredJK =
              displacementXJK * displacementXJK + displacementYJK * displacementYJK + displacementZJK * displacementZJK;

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredJK > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision displacementXKI = xi - xk;
          const SoAFloatPrecision displacementYKI = yi - yk;
          const SoAFloatPrecision displacementZKI = zi - zk;
          const SoAFloatPrecision distSquaredKI =
              displacementXKI * displacementXKI + displacementYKI * displacementYKI + displacementZKI * displacementZKI;

          const SoAFloatPrecision distJK = std::sqrt(distSquaredJK);
          const SoAFloatPrecision distKI = std::sqrt(distSquaredKI);

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredKI > cutoffSquared) {
            continue;
          }

          SoAFloatPrecision forceIX, forceIY, forceIZ;
          SoAFloatPrecision forceJX, forceJY, forceJZ;
          SoAFloatPrecision allDistsTripled, cosines, expTerm, sum;
          SoAKernelN3(distIJ, distJK, distKI, displacementXIJ, displacementYIJ, displacementZIJ, displacementXJK,
                      displacementYJK, displacementZJK, displacementXKI, displacementYKI, displacementZKI,
                      distSquaredIJ, distSquaredJK, distSquaredKI, const_nu, const_alpha, forceIX, forceIY, forceIZ,
                      forceJX, forceJY, forceJZ, allDistsTripled, cosines, expTerm, sum);

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
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

            // Virial is calculated as f_i * r_i
            // see Thompson et al.: https://doi.org/10.1063/1.3245303

            if (ownedStateI == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy;
              virialSumX += forceIX * xi;
              virialSumY += forceIY * yi;
              virialSumZ += forceIZ * zi;
            }
            if (ownedStateJ == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy;
              virialSumX += forceJX * xj;
              virialSumY += forceJY * yj;
              virialSumZ += forceJZ * zj;
            }
            if (ownedStateK == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy;
              virialSumX += forceKX * xk;
              virialSumY += forceKY * yk;
              virialSumZ += forceKZ * zk;
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

      constexpr size_t numFLOPsPerDistanceCall = 0;
      constexpr size_t numFLOPsPerN3KernelCall = 0;
      constexpr size_t numFLOPsPerNoN3KernelCall = 0;
      constexpr size_t numFLOPsPerN3GlobalCalc = 0;
      constexpr size_t numFLOPsPerNoN3GlobalCalc = 0;

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
    const SoAFloatPrecision const_alpha = _alpha;

    size_t soa1Size = soa1.size();
    size_t soa2Size = soa2.size();
    if constexpr (countFLOPs) {
      numTripletsCountingSum = (soa1Size * soa2Size * (soa1Size + soa2Size - 2)) / 2.;
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

        const SoAFloatPrecision xj = xptr2[j];
        const SoAFloatPrecision yj = yptr2[j];
        const SoAFloatPrecision zj = zptr2[j];

        const SoAFloatPrecision displacementXIJ = xj - xi;
        const SoAFloatPrecision displacementYIJ = yj - yi;
        const SoAFloatPrecision displacementZIJ = zj - zi;
        const SoAFloatPrecision distSquaredIJ =
            displacementXIJ * displacementXIJ + displacementYIJ * displacementYIJ + displacementZIJ * displacementZIJ;

        if constexpr (countFLOPs) {
          ++numDistanceCalculationSum;
        }

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        for (unsigned int k = j + 1; k < soa2.size(); ++k) {
          const auto ownedStateK = ownedStatePtr2[k];
          if (ownedStateK == autopas::OwnershipState::dummy) {
            continue;
          }

          const SoAFloatPrecision xk = xptr2[k];
          const SoAFloatPrecision yk = yptr2[k];
          const SoAFloatPrecision zk = zptr2[k];

          const SoAFloatPrecision displacementXJK = xk - xj;
          const SoAFloatPrecision displacementYJK = yk - yj;
          const SoAFloatPrecision displacementZJK = zk - zj;
          const SoAFloatPrecision distSquaredJK =
              displacementXJK * displacementXJK + displacementYJK * displacementYJK + displacementZJK * displacementZJK;

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }

          // Due to low 3-body hitrate, mostly better than masking
          if (distSquaredJK > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision displacementXKI = xi - xk;
          const SoAFloatPrecision displacementYKI = yi - yk;
          const SoAFloatPrecision displacementZKI = zi - zk;
          const SoAFloatPrecision distSquaredKI =
              displacementXKI * displacementXKI + displacementYKI * displacementYKI + displacementZKI * displacementZKI;

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredKI > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision distIJ = std::sqrt(distSquaredIJ);
          const SoAFloatPrecision distJK = std::sqrt(distSquaredJK);
          const SoAFloatPrecision distKI = std::sqrt(distSquaredKI);

          if constexpr (not newton3) {
            // Compute only force I without Newton3
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision allDistsTripled, cosines, expTerm, sum;
            SoAKernelNoN3(distIJ, distJK, distKI, displacementXIJ, displacementYIJ, displacementZIJ, displacementXJK,
                          displacementYJK, displacementZJK, displacementXKI, displacementYKI, displacementZKI,
                          distSquaredIJ, distSquaredJK, distSquaredKI, const_nu, const_alpha, forceIX, forceIY, forceIZ,
                          allDistsTripled, cosines, expTerm, sum);

            fXAccI += forceIX;
            fYAccI += forceIY;
            fZAccI += forceIZ;

            if constexpr (countFLOPs) {
              ++numKernelCallsNoN3Sum;
            }
            if constexpr (calculateGlobals) {
              const SoAFloatPrecision potentialEnergy3 = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

              if (ownedStateI == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceIX * xi;
                virialSumY += forceIY * yi;
                virialSumZ += forceIZ * zi;
              }
              if constexpr (countFLOPs) {
                ++numGlobalCalcsNoN3Sum;
              }
            }
          } else {
            // Compute all forces with Newton3
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision forceJX, forceJY, forceJZ;
            SoAFloatPrecision allDistsTripled, cosines, expTerm, sum;
            SoAKernelN3(distIJ, distJK, distKI, displacementXIJ, displacementYIJ, displacementZIJ, displacementXJK,
                        displacementYJK, displacementZJK, displacementXKI, displacementYKI, displacementZKI,
                        distSquaredIJ, distSquaredJK, distSquaredKI, const_nu, const_alpha, forceIX, forceIY, forceIZ,
                        forceJX, forceJY, forceJZ, allDistsTripled, cosines, expTerm, sum);

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
              const SoAFloatPrecision potentialEnergy3 = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

              if (ownedStateI == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceIX * xi;
                virialSumY += forceIY * yi;
                virialSumZ += forceIZ * zi;
              }
              if (ownedStateJ == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceJX * xj;
                virialSumY += forceJY * yj;
                virialSumZ += forceJZ * zj;
              }
              if (ownedStateK == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceKX * xk;
                virialSumY += forceKY * yk;
                virialSumZ += forceKZ * zk;
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

        const SoAFloatPrecision xj = xptr1[j];
        const SoAFloatPrecision yj = yptr1[j];
        const SoAFloatPrecision zj = zptr1[j];

        const SoAFloatPrecision displacementXIJ = xj - xi;
        const SoAFloatPrecision displacementYIJ = yj - yi;
        const SoAFloatPrecision displacementZIJ = zj - zi;
        const SoAFloatPrecision distSquaredIJ =
            displacementXIJ * displacementXIJ + displacementYIJ * displacementYIJ + displacementZIJ * displacementZIJ;

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

          const SoAFloatPrecision xk = xptr2[k];
          const SoAFloatPrecision yk = yptr2[k];
          const SoAFloatPrecision zk = zptr2[k];

          const SoAFloatPrecision displacementXJK = xk - xj;
          const SoAFloatPrecision displacementYJK = yk - yj;
          const SoAFloatPrecision displacementZJK = zk - zj;
          const SoAFloatPrecision distSquaredJK =
              displacementXJK * displacementXJK + displacementYJK * displacementYJK + displacementZJK * displacementZJK;

          // Due to low 3-body hitrate, mostly better than masking
          if (distSquaredJK > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision displacementXKI = xi - xk;
          const SoAFloatPrecision displacementYKI = yi - yk;
          const SoAFloatPrecision displacementZKI = zi - zk;
          const SoAFloatPrecision distSquaredKI =
              displacementXKI * displacementXKI + displacementYKI * displacementYKI + displacementZKI * displacementZKI;

          if (distSquaredKI > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision distIJ = std::sqrt(distSquaredIJ);
          const SoAFloatPrecision distJK = std::sqrt(distSquaredJK);
          const SoAFloatPrecision distKI = std::sqrt(distSquaredKI);

          SoAFloatPrecision forceIX, forceIY, forceIZ;
          SoAFloatPrecision forceJX, forceJY, forceJZ;
          SoAFloatPrecision allDistsTripled, cosines, expTerm, sum;
          SoAKernelN3(distIJ, distJK, distKI, displacementXIJ, displacementYIJ, displacementZIJ, displacementXJK,
                      displacementYJK, displacementZJK, displacementXKI, displacementYKI, displacementZKI,
                      distSquaredIJ, distSquaredJK, distSquaredKI, const_nu, const_alpha, forceIX, forceIY, forceIZ,
                      forceJX, forceJY, forceJZ, allDistsTripled, cosines, expTerm, sum);

          fXAccI += forceIX;
          fYAccI += forceIY;
          fZAccI += forceIZ;

          fXAccJ += forceJX;
          fYAccJ += forceJY;
          fZAccJ += forceJZ;

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
              const SoAFloatPrecision potentialEnergy3 = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);
              if (ownedStateK == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy3;
                virialSumX += forceKX * xk;
                virialSumY += forceKY * yk;
                virialSumZ += forceKZ * zk;
              }
            }
          }

          if constexpr (calculateGlobals) {
            const SoAFloatPrecision potentialEnergy3 = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);
            if (ownedStateI == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy3;
              virialSumX += forceIX * xi;
              virialSumY += forceIY * yi;
              virialSumZ += forceIZ * zi;
            }
            if (ownedStateJ == autopas::OwnershipState::owned) {
              potentialEnergySum += potentialEnergy3;
              virialSumX += forceJX * xj;
              virialSumY += forceJY * yj;
              virialSumZ += forceJZ * zj;
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
    auto *const __restrict fxptr3 = soa3.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr3 = soa3.template begin<Particle_T::AttributeNames::forceY>();
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
    const SoAFloatPrecision const_alpha = _alpha;

    const auto soa1Size = soa1.size();
    const auto soa2Size = soa2.size();
    const auto soa3Size = soa3.size();
    if constexpr (countFLOPs) {
      numTripletsCountingSum = soa1Size * soa2Size * soa3Size;
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

        const SoAFloatPrecision displacementXIJ = xj - xi;
        const SoAFloatPrecision displacementYIJ = yj - yi;
        const SoAFloatPrecision displacementZIJ = zj - zi;
        const SoAFloatPrecision distSquaredIJ =
            displacementXIJ * displacementXIJ + displacementYIJ * displacementYIJ + displacementZIJ * displacementZIJ;

        if constexpr (countFLOPs) {
          ++numDistanceCalculationSum;
        }
        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        SoAFloatPrecision fXAccJ = 0.;
        SoAFloatPrecision fYAccJ = 0.;
        SoAFloatPrecision fZAccJ = 0.;

        for (unsigned int k = 0; k < soa3Size; ++k) {
          const auto ownedStateK = ownedStatePtr3[k];
          if (ownedStateK == autopas::OwnershipState::dummy) {
            continue;
          }

          const SoAFloatPrecision xk = xptr3[k];
          const SoAFloatPrecision yk = yptr3[k];
          const SoAFloatPrecision zk = zptr3[k];

          const SoAFloatPrecision displacementXJK = xk - xj;
          const SoAFloatPrecision displacementYJK = yk - yj;
          const SoAFloatPrecision displacementZJK = zk - zj;
          const SoAFloatPrecision distSquaredJK =
              displacementXJK * displacementXJK + displacementYJK * displacementYJK + displacementZJK * displacementZJK;

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredJK > cutoffSquared) {
            continue;
          }

          const SoAFloatPrecision displacementXKI = xi - xk;
          const SoAFloatPrecision displacementYKI = yi - yk;
          const SoAFloatPrecision displacementZKI = zi - zk;
          const SoAFloatPrecision distSquaredKI =
              displacementXKI * displacementXKI + displacementYKI * displacementYKI + displacementZKI * displacementZKI;

          const SoAFloatPrecision distIJ = std::sqrt(distSquaredIJ);
          const SoAFloatPrecision distJK = std::sqrt(distSquaredJK);
          const SoAFloatPrecision distKI = std::sqrt(distSquaredKI);

          if constexpr (countFLOPs) {
            ++numDistanceCalculationSum;
          }
          if (distSquaredKI > cutoffSquared) {
            continue;
          }

          if constexpr (not newton3) {
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision allDistsTripled, cosines, expTerm, sum;
            SoAKernelNoN3(distIJ, distJK, distKI, displacementXIJ, displacementYIJ, displacementZIJ, displacementXJK,
                          displacementYJK, displacementZJK, displacementXKI, displacementYKI, displacementZKI,
                          distSquaredIJ, distSquaredJK, distSquaredKI, const_nu, const_alpha, forceIX, forceIY, forceIZ,
                          allDistsTripled, cosines, expTerm, sum);

            fXAccI += forceIX;
            fYAccI += forceIY;
            fZAccI += forceIZ;

            if constexpr (countFLOPs) {
              ++numKernelCallsNoN3Sum;
            }
            if constexpr (calculateGlobals) {
              const SoAFloatPrecision potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

              if (ownedStateI == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy;
                virialSumX += forceIX * xi;
                virialSumY += forceIY * yi;
                virialSumZ += forceIZ * zi;
              }

              if constexpr (countFLOPs) {
                ++numGlobalCalcsNoN3Sum;
              }
            }

          } else {
            SoAFloatPrecision forceIX, forceIY, forceIZ;
            SoAFloatPrecision forceJX, forceJY, forceJZ;
            SoAFloatPrecision allDistsTripled, cosines, expTerm, sum;
            SoAKernelN3(distIJ, distJK, distKI, displacementXIJ, displacementYIJ, displacementZIJ, displacementXJK,
                        displacementYJK, displacementZJK, displacementXKI, displacementYKI, displacementZKI,
                        distSquaredIJ, distSquaredJK, distSquaredKI, const_nu, const_alpha, forceIX, forceIY, forceIZ,
                        forceJX, forceJY, forceJZ, allDistsTripled, cosines, expTerm, sum);

            fXAccI += forceIX;
            fYAccI += forceIY;
            fZAccI += forceIZ;

            fXAccJ += forceJX;
            fYAccJ += forceJY;
            fZAccJ += forceJZ;

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
              const SoAFloatPrecision potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);
              if (ownedStateI == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy;
                virialSumX += forceIX * xi;
                virialSumY += forceIY * yi;
                virialSumZ += forceIZ * zi;
              }
              if (ownedStateJ == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy;
                virialSumX += forceJX * xj;
                virialSumY += forceJY * yj;
                virialSumZ += forceJZ * zj;
              }
              if (ownedStateK == autopas::OwnershipState::owned) {
                potentialEnergySum += potentialEnergy;
                virialSumX += forceKX * xk;
                virialSumY += forceKY * yk;
                virialSumZ += forceKZ * zk;
              }
              if constexpr (countFLOPs) {
                ++numGlobalCalcsN3Sum;
              }
            }
          }
        }
        fxptr2[j] += fXAccJ;
        fyptr2[j] += fYAccJ;
        fzptr2[j] += fZAccJ;
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

  /**
   * Inline helper to compute force components for particle I.
   * Returns tuple of (forceIX, forceIY, forceIZ, factor, allDotProducts, allDistsSquared)
   */
  __attribute__((always_inline)) inline auto SoAKernelNoN3(
      const SoAFloatPrecision &distIJ, const SoAFloatPrecision &distJK, const SoAFloatPrecision &distKI,
      const SoAFloatPrecision &displacementXIJ, const SoAFloatPrecision &displacementYIJ,
      const SoAFloatPrecision &displacementZIJ, const SoAFloatPrecision &displacementXJK,
      const SoAFloatPrecision &displacementYJK, const SoAFloatPrecision &displacementZJK,
      const SoAFloatPrecision &displacementXKI, const SoAFloatPrecision &displacementYKI,
      const SoAFloatPrecision &displacementZKI, const SoAFloatPrecision &distSquaredIJ,
      const SoAFloatPrecision &distSquaredJK, const SoAFloatPrecision &distSquaredKI, const SoAFloatPrecision &nu,
      const SoAFloatPrecision &alpha, SoAFloatPrecision &forceIX, SoAFloatPrecision &forceIY,
      SoAFloatPrecision &forceIZ, SoAFloatPrecision &allDistsTripled, SoAFloatPrecision &cosines,
      SoAFloatPrecision &expTerm, SoAFloatPrecision &sum) const {
    // Numerators of cosine representation (cos_i = (r_ij^2 + r_ik^2 - r_jk^2) / (2 * r_ij * r_ik)
    const double numKI = distSquaredIJ + distSquaredJK - distSquaredKI;
    const double numJK = distSquaredIJ + distSquaredKI - distSquaredJK;
    const double numIJ = distSquaredJK + distSquaredKI - distSquaredIJ;

    const double numerator = numKI * numJK * numIJ;

    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDists = distIJ * distJK * distKI;
    allDistsTripled = allDistsSquared * allDists;

    // Gradient factors of 1. / (rrr)^3
    const double allDistsTriplesGradientIJ = 3. / (allDistsTripled * distSquaredIJ);
    const double allDistsTriplesGradientKI = -3. / (allDistsTripled * distSquaredKI);

    // Product of all cosines multiplied with 3: 3 * cos(a)cos(b)cos(c)
    cosines = (3. / 8.) * numerator / allDistsSquared;
    const double cosinesGradientIJ =
        (3. / 4.) * ((numerator / distSquaredIJ - numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);
    const double cosinesGradientKI =
        (3. / 4.) * ((-numerator / distSquaredKI + numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);

    // Gradient factors corresponding to the normal ATM term
    const auto fullATMGradientIJ =
        nu * ((1. + cosines) * allDistsTriplesGradientIJ + cosinesGradientIJ / allDistsTripled);
    const auto fullATMGradientKI =
        nu * ((1. + cosines) * allDistsTriplesGradientKI + cosinesGradientKI / allDistsTripled);

    expTerm = std::exp(-alpha * (distIJ + distJK + distKI));

    // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
    std::array<double, 6> sumFactors{};
    sum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      sumFactors[n] = _constantsA[n] * std::pow((distIJ * distJK * distKI), 2. * n / 3.);
      sum += sumFactors[n];
    }

    // Gradient factor of the sum in ij-direction
    double ijSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - alpha);
    }

    // Gradient factor of the sum in ki-direction
    double kiSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      kiSum += sumFactors[n] * (2. * n / (3. * distKI) - alpha);
    }

    // Total gradient factors for the exponential term times the cosines term
    const double fullExpGradientIJ = expTerm * (-(1. + cosines) * ijSum / distIJ + cosinesGradientIJ * sum);
    const double fullExpGradientKI = expTerm * ((1. + cosines) * kiSum / distKI + cosinesGradientKI * sum);

    // Assembling the forces
    const auto forceIDirectionXIJ = displacementXIJ * (fullATMGradientIJ + fullExpGradientIJ);
    const auto forceIDirectionYIJ = displacementYIJ * (fullATMGradientIJ + fullExpGradientIJ);
    const auto forceIDirectionZIJ = displacementZIJ * (fullATMGradientIJ + fullExpGradientIJ);
    const auto forceIDirectionXKI = displacementXKI * (fullATMGradientKI + fullExpGradientKI);
    const auto forceIDirectionYKI = displacementYKI * (fullATMGradientKI + fullExpGradientKI);
    const auto forceIDirectionZKI = displacementZKI * (fullATMGradientKI + fullExpGradientKI);

    forceIX = (forceIDirectionXIJ + forceIDirectionXKI) * (-1.0);
    forceIY = (forceIDirectionYIJ + forceIDirectionYKI) * (-1.0);
    forceIZ = (forceIDirectionZIJ + forceIDirectionZKI) * (-1.0);
  }

  /**
   * Inline helper to compute force components for particle I.
   * Returns tuple of (forceIX, forceIY, forceIZ, factor, allDotProducts, allDistsSquared)
   */
  __attribute__((always_inline)) inline auto SoAKernelN3(
      const SoAFloatPrecision &distIJ, const SoAFloatPrecision &distJK, const SoAFloatPrecision &distKI,
      const SoAFloatPrecision &displacementXIJ, const SoAFloatPrecision &displacementYIJ,
      const SoAFloatPrecision &displacementZIJ, const SoAFloatPrecision &displacementXJK,
      const SoAFloatPrecision &displacementYJK, const SoAFloatPrecision &displacementZJK,
      const SoAFloatPrecision &displacementXKI, const SoAFloatPrecision &displacementYKI,
      const SoAFloatPrecision &displacementZKI, const SoAFloatPrecision &distSquaredIJ,
      const SoAFloatPrecision &distSquaredJK, const SoAFloatPrecision &distSquaredKI, const SoAFloatPrecision &nu,
      const SoAFloatPrecision &alpha, SoAFloatPrecision &forceIX, SoAFloatPrecision &forceIY,
      SoAFloatPrecision &forceIZ, SoAFloatPrecision &forceJX, SoAFloatPrecision &forceJY, SoAFloatPrecision &forceJZ,

      SoAFloatPrecision &allDistsTripled, SoAFloatPrecision &cosines, SoAFloatPrecision &expTerm,
      SoAFloatPrecision &sum) const {
    // Numerators of cosine representation (cos_i = (r_ij^2 + r_ik^2 - r_jk^2) / (2 * r_ij * r_ik)
    const double numKI = distSquaredIJ + distSquaredJK - distSquaredKI;
    const double numJK = distSquaredIJ + distSquaredKI - distSquaredJK;
    const double numIJ = distSquaredJK + distSquaredKI - distSquaredIJ;

    const double numerator = numKI * numJK * numIJ;

    const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    const double allDists = distIJ * distJK * distKI;
    allDistsTripled = allDistsSquared * allDists;

    // Gradient factors of 1. / (rrr)^3
    const double allDistsTriplesGradientIJ = 3. / (allDistsTripled * distSquaredIJ);
    const double allDistsTriplesGradientKI = -3. / (allDistsTripled * distSquaredKI);

    // Product of all cosines multiplied with 3: 3 * cos(a)cos(b)cos(c)
    cosines = (3. / 8.) * numerator / allDistsSquared;
    const double cosinesGradientIJ =
        (3. / 4.) * ((numerator / distSquaredIJ - numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);
    const double cosinesGradientKI =
        (3. / 4.) * ((-numerator / distSquaredKI + numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);

    // Gradient factors corresponding to the normal ATM term
    const auto fullATMGradientIJ =
        nu * ((1. + cosines) * allDistsTriplesGradientIJ + cosinesGradientIJ / allDistsTripled);
    const auto fullATMGradientKI =
        nu * ((1. + cosines) * allDistsTriplesGradientKI + cosinesGradientKI / allDistsTripled);

    expTerm = std::exp(-alpha * (distIJ + distJK + distKI));

    // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
    std::array<double, 6> sumFactors{};
    sum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      sumFactors[n] = _constantsA[n] * std::pow((distIJ * distJK * distKI), 2. * n / 3.);
      sum += sumFactors[n];
    }

    // Gradient factor of the sum in ij-direction
    double ijSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - alpha);
    }

    // Gradient factor of the sum in ki-direction
    double kiSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      kiSum += sumFactors[n] * (2. * n / (3. * distKI) - alpha);
    }

    // Total gradient factors for the exponential term times the cosines term
    const double fullExpGradientIJ = expTerm * (-(1. + cosines) * ijSum / distIJ + cosinesGradientIJ * sum);
    const double fullExpGradientKI = expTerm * ((1. + cosines) * kiSum / distKI + cosinesGradientKI * sum);

    // Assembling the forces
    const auto forceIDirectionXIJ = displacementXIJ * (fullATMGradientIJ + fullExpGradientIJ);
    const auto forceIDirectionYIJ = displacementYIJ * (fullATMGradientIJ + fullExpGradientIJ);
    const auto forceIDirectionZIJ = displacementZIJ * (fullATMGradientIJ + fullExpGradientIJ);
    const auto forceIDirectionXKI = displacementXKI * (fullATMGradientKI + fullExpGradientKI);
    const auto forceIDirectionYKI = displacementYKI * (fullATMGradientKI + fullExpGradientKI);
    const auto forceIDirectionZKI = displacementZKI * (fullATMGradientKI + fullExpGradientKI);

    forceIX = (forceIDirectionXIJ + forceIDirectionXKI) * (-1.0);
    forceIY = (forceIDirectionYIJ + forceIDirectionYKI) * (-1.0);
    forceIZ = (forceIDirectionZIJ + forceIDirectionZKI) * (-1.0);

    // Calculate all components for jk-direction
    const double allDistsTriplesGradientJK = 3. / (allDistsTripled * distSquaredJK);
    const double cosinesGradientJK =
        (3. / 4.) * ((numerator / distSquaredJK + numKI * numIJ - numJK * numIJ - numJK * numKI) / allDistsSquared);
    const auto fullATMGradientJK =
        _nu * ((1. + cosines) * allDistsTriplesGradientJK + cosinesGradientJK / allDistsTripled);

    double jkSum = 0.0;
    for (auto n = 0; n < sumFactors.size(); n++) {
      jkSum += sumFactors[n] * (2. * n / (3. * distJK) - _alpha);
    }
    const double fullExpGradientJK = expTerm * (-(1. + cosines) * jkSum / distJK + cosinesGradientJK * sum);

    // Assembling the forces
    const auto forceJDirectionXIJ = displacementXIJ * (-fullATMGradientIJ - fullExpGradientIJ);
    const auto forceJDirectionYIJ = displacementYIJ * (-fullATMGradientIJ - fullExpGradientIJ);
    const auto forceJDirectionZIJ = displacementZIJ * (-fullATMGradientIJ - fullExpGradientIJ);
    const auto forceJDirectionXJK = displacementXJK * (fullATMGradientJK + fullExpGradientJK);
    const auto forceJDirectionYJK = displacementYJK * (fullATMGradientJK + fullExpGradientJK);
    const auto forceJDirectionZJK = displacementZJK * (fullATMGradientJK + fullExpGradientJK);

    forceJX = (forceJDirectionXIJ + forceJDirectionXJK) * (-1.0);
    forceJY = (forceJDirectionYIJ + forceJDirectionYJK) * (-1.0);
    forceJZ = (forceJDirectionZIJ + forceJDirectionZJK) * (-1.0);
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

  // Parameters of the extended Axilrod-Teller potential for Krypton in Kelvin (K) and Angström (A)
  const double _nu = 1.61525e6;    // K*A^9
  const double _alpha = 1.378382;  // A^-1
  // Units: {K, K*A^-2, K*A^-4, K*A^-6, K*A^-8, K*A^-10}
  const std::array<double, 6> _constantsA = {-0.3081304e8, -0.3519442e8, 0.4928052e7, -0.2182411e6, 0.343088e4, 0.0};

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

/**
 * @file KryptonExtendedATMFunctor.h
 * @author muehlhaeusser
 * @date 08.08.2024
 */

#pragma once

#include <immintrin.h>

#include <array>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/Timer.h"
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
 * @tparam Particle The type of particle.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 */
template <class Particle, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool countFLOPs = false>
class KryptonExtendedATMFunctor
    : public autopas::TriwiseFunctor<Particle,
                                     KryptonExtendedATMFunctor<Particle, useNewton3, calculateGlobals, countFLOPs>> {
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
  KryptonExtendedATMFunctor() = delete;

  /**
   * Constructor
   * @param cutoff
   */
  explicit KryptonExtendedATMFunctor(double cutoff)
      : autopas::TriwiseFunctor<Particle,
                                KryptonExtendedATMFunctor<Particle, useNewton3, calculateGlobals, countFLOPs>>(cutoff),
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

  void AoSFunctor(Particle &i, Particle &j, Particle &k, bool newton3) override {
    using namespace autopas::utils::ArrayMath::literals;

    // if (i.isDummy() or j.isDummy() or k.isDummy()) {
    //   std::cout << "DUMMIEEEES" << std::endl;
    //   return;
    // }

    const auto threadnum = autopas::autopas_get_thread_num();
    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
    }

    const auto iPos = i.getR();
    const auto jPos = j.getR();
    const auto kPos = k.getR();

    const auto displacementIJ = jPos - iPos;
    const double distSquaredIJ = autopas::utils::ArrayMath::dot(displacementIJ, displacementIJ);

    const auto displacementJK = kPos - jPos;
    const double distSquaredJK = autopas::utils::ArrayMath::dot(displacementJK, displacementJK);

    const auto displacementKI = iPos - kPos;
    const double distSquaredKI = autopas::utils::ArrayMath::dot(displacementKI, displacementKI);

    // Check cutoff for every distance
    if (distSquaredIJ > _cutoffSquared or distSquaredJK > _cutoffSquared or distSquaredKI > _cutoffSquared) {
      return;
    }
    // _aosThreadDataGlobals[threadnum]._timerSingle.start();

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
    const double numKIIJ = numKI * numIJ;
    const double numJKIJ = numJK * numIJ;
    const double numJKKI = numJK * numKI;
    const double cosinesGradientIJ =
        (3. / 4.) * ((numerator / distSquaredIJ - numKIIJ - numJKIJ + numJKKI) / allDistsSquared);
    const double cosinesGradientKI =
        (3. / 4.) * ((-numerator / distSquaredKI + numKIIJ - numJKIJ + numJKKI) / allDistsSquared);

    // Gradient factors corresponding to the normal ATM term
    const auto fullATMGradientIJ =
        _nu * ((1. + cosines) * allDistsTriplesGradientIJ + cosinesGradientIJ / allDistsTripled);
    const auto fullATMGradientKI =
        _nu * ((1. + cosines) * allDistsTriplesGradientKI + cosinesGradientKI / allDistsTripled);

    const double expTerm = std::exp(-_alpha * (distIJ + distJK + distKI));

    // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
    std::array<double, 6> sumFactors{};
    const double cbrtAllDists = std::cbrt(allDists);
    sumFactors[0] = _constantsA[0];
    const double cbrtAllDistsSquared = cbrtAllDists * cbrtAllDists;
    sumFactors[1] = _constantsA[1] * cbrtAllDistsSquared;
    const double cbrtAllDistsQuadrupled = cbrtAllDistsSquared * cbrtAllDistsSquared;
    sumFactors[2] = _constantsA[2] * cbrtAllDistsQuadrupled;
    sumFactors[3] = _constantsA[3] * allDistsSquared;
    const double cbrtAllDistsOctupled = cbrtAllDistsQuadrupled * cbrtAllDistsQuadrupled;
    sumFactors[4] = _constantsA[4] * cbrtAllDistsOctupled;
    const double cbrtAllDistsDecupled = cbrtAllDistsOctupled * cbrtAllDistsSquared;
    sumFactors[5] = _constantsA[5] * cbrtAllDistsDecupled;
    const double sum = sumFactors[0] + sumFactors[1] + sumFactors[2] + sumFactors[3] + sumFactors[4] + sumFactors[5];

    // _aosThreadDataGlobals[threadnum]._timerTriple.start();
    // Gradient factor of the sum in ij-direction
    double ijSum = 0.0;
#pragma omp simd reduction(+ : ijSum)
    for (int8_t n = 0; n < 6; n++) {
      ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - _alpha);
    }

    // Gradient factor of the sum in ki-direction
    double kiSum = 0.0;
#pragma omp simd reduction(+ : kiSum)
    for (int8_t n = 0; n < 6; n++) {
      kiSum += sumFactors[n] * (2. * n / (3. * distKI) - _alpha);
    }
    // _aosThreadDataGlobals[threadnum]._timerTriple.stop();

    // Total gradient factors for the exponential term times the cosines term
    const double fullExpGradientIJ = expTerm * (-(1. + cosines) * ijSum / distIJ + cosinesGradientIJ * sum);
    const double fullExpGradientKI = expTerm * ((1. + cosines) * kiSum / distKI + cosinesGradientKI * sum);

    // Assembling the forces
    const auto forceIDirectionIJ = displacementIJ * (fullATMGradientIJ + fullExpGradientIJ);
    const auto forceIDirecationKI = displacementKI * (fullATMGradientKI + fullExpGradientKI);

    const auto forceI = (forceIDirectionIJ + forceIDirecationKI) * (-1.0);

    i.addF(forceI);
    //
    // auto forceJ = forceI;
    // auto forceK = forceI;
    // if (newton3) {
    //   // Calculate all components for jk-direction
    //   const double allDistsTriplesGradientJK = 3. / (allDistsTripled * distSquaredJK);
    //   const double cosinesGradientJK =
    //       (3. / 4.) * ((numerator / distSquaredJK + numKI * numIJ - numJK * numIJ - numJK * numKI) / allDistsSquared);
    //   const auto fullATMGradientJK =
    //       _nu * ((1. + cosines) * allDistsTriplesGradientJK + cosinesGradientJK / allDistsTripled);
    //
    //   double jkSum = 0.0;
    //   for (auto n = 0; n < sumFactors.size(); n++) {
    //     jkSum += sumFactors[n] * (2. * n / (3. * distJK) - _alpha);
    //   }
    //   const double fullExpGradientJK = expTerm * (-(1. + cosines) * jkSum / distJK + cosinesGradientJK * sum);
    //
    //   // Assembling the forces
    //   const auto forceJDirectionIJ = displacementIJ * (-fullATMGradientIJ - fullExpGradientIJ);
    //   const auto forceJDirectionJK = displacementJK * (fullATMGradientJK + fullExpGradientJK);
    //
    //   forceJ = (forceJDirectionIJ + forceJDirectionJK) * (-1.0);
    //   j.addF(forceJ);
    //
    //   // Using newton's third law for the force on particle k
    //   forceK = (forceI + forceJ) * (-1.0);
    //   k.addF(forceK);
    // }

    if constexpr (countFLOPs) {
      if (newton3) {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsN3;
      } else {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsNoN3;
      }
    }

    if constexpr (calculateGlobals) {
      // _aosThreadDataGlobals[threadnum]._timerPair.start();

      // Add 3 * potential energy to every owned particle of the interaction.
      // Division to the correct value is handled in endTraversal().
      const double potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

      // Virial is calculated as f_i * r_i
      // see Thompson et al.: https://doi.org/10.1063/1.3245303
      const auto virialI = forceI * iPos;
      if (i.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadDataGlobals[threadnum].virialSum += virialI;
      }
      // for non-newton3 particles j and/or k will be considered in a separate calculation
      // if (newton3 and j.isOwned()) {
      //   const auto virialJ = forceJ * j.getR();
      //   _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
      //   _aosThreadDataGlobals[threadnum].virialSum += virialJ;
      // }
      // if (newton3 and k.isOwned()) {
      //   const auto virialK = forceK * k.getR();
      //   _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
      //   _aosThreadDataGlobals[threadnum].virialSum += virialK;
      // }
      if constexpr (countFLOPs) {
        if (newton3) {
          ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsN3;
        } else {
          ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3;
        }
      }
    // _aosThreadDataGlobals[threadnum]._timerPair.stop();
    }
    // _aosThreadDataGlobals[threadnum]._timerSingle.stop();
  }

  double newtonRaphsonCbrt(double x) {
    const double epsilon = 1e-10; // Tolerance for convergence
    double guess = x; // Initial guess

    while (true) {
      double error = (guess * guess * guess - x) / (3 * guess * guess);
      if (std::fabs(error) < epsilon) {
        break;
      }
      guess -= error;
    }

    return guess;
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.size() == 0) return;

    const auto threadnum = autopas::autopas_get_thread_num();
    _aosThreadDataGlobals[threadnum]._timerSingle.start();

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;  // Note: This is not the potential energy but some fixed multiple of it.
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    for (unsigned int i = 0; i < soa.size(); ++i) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (unsigned int k = j + 1; k < soa.size(); ++k) {
          const auto ownedStateK = ownedStatePtr[k];

          const auto displacementIJx = xptr[j] - xptr[i];
          const auto displacementIJy = yptr[j] - yptr[i];
          const auto displacementIJz = zptr[j] - zptr[i];
          const auto displacementIJx2 = displacementIJx * displacementIJx;
          const auto displacementIJy2 = displacementIJy * displacementIJy;
          const auto displacementIJz2 = displacementIJz * displacementIJz;

          const auto displacementJKx = xptr[k] - xptr[j];
          const auto displacementJKy = yptr[k] - yptr[j];
          const auto displacementJKz = zptr[k] - zptr[j];
          const auto displacementJKx2 = displacementJKx * displacementJKx;
          const auto displacementJKy2 = displacementJKy * displacementJKy;
          const auto displacementJKz2 = displacementJKz * displacementJKz;

          const auto displacementKIx = xptr[i] - xptr[k];
          const auto displacementKIy = yptr[i] - yptr[k];
          const auto displacementKIz = zptr[i] - zptr[k];
          const auto displacementKIx2 = displacementKIx * displacementKIx;
          const auto displacementKIy2 = displacementKIy * displacementKIy;
          const auto displacementKIz2 = displacementKIz * displacementKIz;

          const auto distSquaredIJ = displacementIJx2 + displacementIJy2 + displacementIJz2;
          const auto distSquaredJK = displacementJKx2 + displacementJKy2 + displacementJKz2;
          const auto distSquaredKI = displacementKIx2 + displacementKIy2 + displacementKIz2;

          const bool mask = distSquaredIJ <= cutoffSquared and distSquaredJK <= cutoffSquared
            and distSquaredKI <= cutoffSquared  and ownedStateK != autopas::OwnershipState::dummy;

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
          const auto forceIDirectionIJx = displacementIJx * (fullATMGradientIJ + fullExpGradientIJ);
          const auto forceIDirectionIJy = displacementIJy * (fullATMGradientIJ + fullExpGradientIJ);
          const auto forceIDirectionIJz = displacementIJz * (fullATMGradientIJ + fullExpGradientIJ);

          const auto forceIDirecationKIx = displacementKIx * (fullATMGradientKI + fullExpGradientKI);
          const auto forceIDirecationKIy = displacementKIy * (fullATMGradientKI + fullExpGradientKI);
          const auto forceIDirecationKIz = displacementKIz * (fullATMGradientKI + fullExpGradientKI);

          const auto forceIx = (forceIDirectionIJx + forceIDirecationKIx) * mask * (-1.0);
          const auto forceIy = (forceIDirectionIJy + forceIDirecationKIy) * mask * (-1.0);
          const auto forceIz = (forceIDirectionIJz + forceIDirecationKIz) * mask * (-1.0);

          fxacc += forceIx;
          fyacc += forceIy;
          fzacc += forceIz;

          if constexpr (calculateGlobals) {
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy = mask * (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

            // Virial is calculated as f_i * r_i
            // see Thompson et al.: https://doi.org/10.1063/1.3245303
            const auto virialIx = mask * forceIx * xptr[i];
            const auto virialIy = mask * forceIy * yptr[i];
            const auto virialIz = mask * forceIz * zptr[i];

            potentialEnergySum += potentialEnergy;
            virialSumX += virialIx;
            virialSumY += virialIy;
            virialSumZ += virialIz;
          }
        }
      }
      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
    _aosThreadDataGlobals[threadnum]._timerSingle.stop();

  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final {

    const auto threadnum = autopas::autopas_get_thread_num();
    _aosThreadDataGlobals[threadnum]._timerPair.start();

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

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      for (unsigned int j = 0; j < soa2.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr2[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (unsigned int k = j + 1; k < soa2.size(); ++k) {
          const auto ownedStateK = ownedStatePtr2[k];

          const auto displacementIJx = x2ptr[j] - x1ptr[i];
          const auto displacementIJy = y2ptr[j] - y1ptr[i];
          const auto displacementIJz = z2ptr[j] - z1ptr[i];
          const auto displacementIJx2 = displacementIJx * displacementIJx;
          const auto displacementIJy2 = displacementIJy * displacementIJy;
          const auto displacementIJz2 = displacementIJz * displacementIJz;

          const auto displacementJKx = x2ptr[k] - x2ptr[j];
          const auto displacementJKy = y2ptr[k] - y2ptr[j];
          const auto displacementJKz = z2ptr[k] - z2ptr[j];
          const auto displacementJKx2 = displacementJKx * displacementJKx;
          const auto displacementJKy2 = displacementJKy * displacementJKy;
          const auto displacementJKz2 = displacementJKz * displacementJKz;

          const auto displacementKIx = x1ptr[i] - x2ptr[k];
          const auto displacementKIy = y1ptr[i] - y2ptr[k];
          const auto displacementKIz = z1ptr[i] - z2ptr[k];
          const auto displacementKIx2 = displacementKIx * displacementKIx;
          const auto displacementKIy2 = displacementKIy * displacementKIy;
          const auto displacementKIz2 = displacementKIz * displacementKIz;

          const auto distSquaredIJ = displacementIJx2 + displacementIJy2 + displacementIJz2;
          const auto distSquaredJK = displacementJKx2 + displacementJKy2 + displacementJKz2;
          const auto distSquaredKI = displacementKIx2 + displacementKIy2 + displacementKIz2;

          const bool mask = distSquaredIJ <= cutoffSquared and distSquaredJK <= cutoffSquared
            and distSquaredKI <= cutoffSquared  and ownedStateK != autopas::OwnershipState::dummy;

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
          const auto forceIDirectionIJx = displacementIJx * (fullATMGradientIJ + fullExpGradientIJ);
          const auto forceIDirectionIJy = displacementIJy * (fullATMGradientIJ + fullExpGradientIJ);
          const auto forceIDirectionIJz = displacementIJz * (fullATMGradientIJ + fullExpGradientIJ);

          const auto forceIDirecationKIx = displacementKIx * (fullATMGradientKI + fullExpGradientKI);
          const auto forceIDirecationKIy = displacementKIy * (fullATMGradientKI + fullExpGradientKI);
          const auto forceIDirecationKIz = displacementKIz * (fullATMGradientKI + fullExpGradientKI);

          const auto forceIx = (forceIDirectionIJx + forceIDirecationKIx) * mask * (-1.0);
          const auto forceIy = (forceIDirectionIJy + forceIDirecationKIy) * mask * (-1.0);
          const auto forceIz = (forceIDirectionIJz + forceIDirecationKIz) * mask * (-1.0);

          fxacc += forceIx;
          fyacc += forceIy;
          fzacc += forceIz;

          if constexpr (calculateGlobals) {
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy = mask * (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

            // Virial is calculated as f_i * r_i
            // see Thompson et al.: https://doi.org/10.1063/1.3245303
            const auto virialIx = mask * forceIx * x1ptr[i];
            const auto virialIy = mask * forceIy * y1ptr[i];
            const auto virialIz = mask * forceIz * z1ptr[i];

            potentialEnergySum += potentialEnergy;
            virialSumX += virialIx;
            virialSumY += virialIy;
            virialSumZ += virialIz;
            // for non-newton3 particles j and/or k will be considered in a separate calculation
            // if (newton3 and ownedStateJ) {
            // const auto virialJ = forceJ * j.getR();
            // _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
            // _aosThreadDataGlobals[threadnum].virialSum += virialJ;
            // }
            // if (newton3 and ownedStateK) {
            //   const auto virialK = forceK * k.getR();
            //   _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
            //   _aosThreadDataGlobals[threadnum].virialSum += virialK;
            // }
            // if constexpr (countFLOPs) {
            //   if (newton3) {
            //     ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsN3;
            //   } else {
            //     ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3;
            //   }
            // }
          }
        }
      }

      // --------------------- //
      for (unsigned int j = i + 1; j < soa1.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr1[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (unsigned int k = 0; k < soa2.size(); ++k) {
          const auto ownedStateK = ownedStatePtr2[k];

          const auto displacementIJx = x1ptr[j] - x1ptr[i];
          const auto displacementIJy = y1ptr[j] - y1ptr[i];
          const auto displacementIJz = z1ptr[j] - z1ptr[i];
          const auto displacementIJx2 = displacementIJx * displacementIJx;
          const auto displacementIJy2 = displacementIJy * displacementIJy;
          const auto displacementIJz2 = displacementIJz * displacementIJz;

          const auto displacementJKx = x2ptr[k] - x1ptr[j];
          const auto displacementJKy = y2ptr[k] - y1ptr[j];
          const auto displacementJKz = z2ptr[k] - z1ptr[j];
          const auto displacementJKx2 = displacementJKx * displacementJKx;
          const auto displacementJKy2 = displacementJKy * displacementJKy;
          const auto displacementJKz2 = displacementJKz * displacementJKz;

          const auto displacementKIx = x1ptr[i] - x2ptr[k];
          const auto displacementKIy = y1ptr[i] - y2ptr[k];
          const auto displacementKIz = z1ptr[i] - z2ptr[k];
          const auto displacementKIx2 = displacementKIx * displacementKIx;
          const auto displacementKIy2 = displacementKIy * displacementKIy;
          const auto displacementKIz2 = displacementKIz * displacementKIz;

          const auto distSquaredIJ = displacementIJx2 + displacementIJy2 + displacementIJz2;
          const auto distSquaredJK = displacementJKx2 + displacementJKy2 + displacementJKz2;
          const auto distSquaredKI = displacementKIx2 + displacementKIy2 + displacementKIz2;

          const bool mask = distSquaredIJ <= cutoffSquared and distSquaredJK <= cutoffSquared
            and distSquaredKI <= cutoffSquared  and ownedStateK != autopas::OwnershipState::dummy;

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
          const auto forceIDirectionIJx = displacementIJx * (fullATMGradientIJ + fullExpGradientIJ);
          const auto forceIDirectionIJy = displacementIJy * (fullATMGradientIJ + fullExpGradientIJ);
          const auto forceIDirectionIJz = displacementIJz * (fullATMGradientIJ + fullExpGradientIJ);

          const auto forceIDirecationKIx = displacementKIx * (fullATMGradientKI + fullExpGradientKI);
          const auto forceIDirecationKIy = displacementKIy * (fullATMGradientKI + fullExpGradientKI);
          const auto forceIDirecationKIz = displacementKIz * (fullATMGradientKI + fullExpGradientKI);

          const auto forceIx = (forceIDirectionIJx + forceIDirecationKIx) * mask * (-1.0);
          const auto forceIy = (forceIDirectionIJy + forceIDirecationKIy) * mask * (-1.0);
          const auto forceIz = (forceIDirectionIJz + forceIDirecationKIz) * mask * (-1.0);

          fxacc += forceIx;
          fyacc += forceIy;
          fzacc += forceIz;

          if constexpr (calculateGlobals) {
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy = mask * (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

            // Virial is calculated as f_i * r_i
            // see Thompson et al.: https://doi.org/10.1063/1.3245303
            const auto virialIx = mask * forceIx * x1ptr[i];
            const auto virialIy = mask * forceIy * y1ptr[i];
            const auto virialIz = mask * forceIz * z1ptr[i];

            potentialEnergySum += potentialEnergy;
            virialSumX += virialIx;
            virialSumY += virialIy;
            virialSumZ += virialIz;
            // for non-newton3 particles j and/or k will be considered in a separate calculation
            // if (newton3 and ownedStateJ) {
            // const auto virialJ = forceJ * j.getR();
            // _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
            // _aosThreadDataGlobals[threadnum].virialSum += virialJ;
            // }
            // if (newton3 and ownedStateK) {
            //   const auto virialK = forceK * k.getR();
            //   _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
            //   _aosThreadDataGlobals[threadnum].virialSum += virialK;
            // }
            // if constexpr (countFLOPs) {
            //   if (newton3) {
            //     ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsN3;
            //   } else {
            //     ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3;
            //   }
            // }
          }
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
    _aosThreadDataGlobals[threadnum]._timerPair.stop();

  }

  void SoAFunctorTriple(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
    autopas::SoAView<SoAArraysType> soa3, const bool newton3) final {
    if (newton3) {
      SoAFunctorTriple<true>(soa1, soa2, soa3);
    } else {
      SoAFunctorTriple<false>(soa1, soa2, soa3);
    }
  }

  template <bool newton3>
  void SoAFunctorTriple(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
  autopas::SoAView<SoAArraysType> soa3) {
    if (soa1.size() == 0 or soa2.size() == 0 or soa3.size() == 0) return;

    const auto threadnum = autopas::autopas_get_thread_num();
    _aosThreadDataGlobals[threadnum]._timerTriple.start();

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x3ptr = soa3.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y3ptr = soa3.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z3ptr = soa3.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr3 = soa3.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx3ptr = soa3.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy3ptr = soa3.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz3ptr = soa3.template begin<Particle::AttributeNames::forceZ>();

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      for (unsigned int j = 0; j < soa2.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr2[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

// #pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ)
        for (unsigned int k = 0; k < soa3.size(); ++k) {
          const auto ownedStateK = ownedStatePtr3[k];

          const auto displacementIJx = x2ptr[j] - x1ptr[i];
          const auto displacementIJy = y2ptr[j] - y1ptr[i];
          const auto displacementIJz = z2ptr[j] - z1ptr[i];
          const auto displacementIJx2 = displacementIJx * displacementIJx;
          const auto displacementIJy2 = displacementIJy * displacementIJy;
          const auto displacementIJz2 = displacementIJz * displacementIJz;

          const auto displacementJKx = x3ptr[k] - x2ptr[j];
          const auto displacementJKy = y3ptr[k] - y2ptr[j];
          const auto displacementJKz = z3ptr[k] - z2ptr[j];
          const auto displacementJKx2 = displacementJKx * displacementJKx;
          const auto displacementJKy2 = displacementJKy * displacementJKy;
          const auto displacementJKz2 = displacementJKz * displacementJKz;

          const auto displacementKIx = x1ptr[i] - x3ptr[k];
          const auto displacementKIy = y1ptr[i] - y3ptr[k];
          const auto displacementKIz = z1ptr[i] - z3ptr[k];
          const auto displacementKIx2 = displacementKIx * displacementKIx;
          const auto displacementKIy2 = displacementKIy * displacementKIy;
          const auto displacementKIz2 = displacementKIz * displacementKIz;

          const auto distSquaredIJ = displacementIJx2 + displacementIJy2 + displacementIJz2;
          const auto distSquaredJK = displacementJKx2 + displacementJKy2 + displacementJKz2;
          const auto distSquaredKI = displacementKIx2 + displacementKIy2 + displacementKIz2;

          if (distSquaredIJ <= _cutoffSquared and distSquaredJK <= _cutoffSquared
            and distSquaredKI <= _cutoffSquared  and ownedStateK != autopas::OwnershipState::dummy) {
            continue;
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
          const auto forceIDirectionIJx = displacementIJx * (fullATMGradientIJ + fullExpGradientIJ);
          const auto forceIDirectionIJy = displacementIJy * (fullATMGradientIJ + fullExpGradientIJ);
          const auto forceIDirectionIJz = displacementIJz * (fullATMGradientIJ + fullExpGradientIJ);

          const auto forceIDirecationKIx = displacementKIx * (fullATMGradientKI + fullExpGradientKI);
          const auto forceIDirecationKIy = displacementKIy * (fullATMGradientKI + fullExpGradientKI);
          const auto forceIDirecationKIz = displacementKIz * (fullATMGradientKI + fullExpGradientKI);

          const auto forceIx = (forceIDirectionIJx + forceIDirecationKIx) * (-1.0);
          const auto forceIy = (forceIDirectionIJy + forceIDirecationKIy) * (-1.0);
          const auto forceIz = (forceIDirectionIJz + forceIDirecationKIz) * (-1.0);

          fxacc += forceIx;
          fyacc += forceIy;
          fzacc += forceIz;

          if constexpr (calculateGlobals) {
            // Add 3 * potential energy to every owned particle of the interaction.
            // Division to the correct value is handled in endTraversal().
            const double potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

            // Virial is calculated as f_i * r_i
            // see Thompson et al.: https://doi.org/10.1063/1.3245303
            const auto virialIx = forceIx * x1ptr[i];
            const auto virialIy = forceIy * y1ptr[i];
            const auto virialIz = forceIz * z1ptr[i];

            potentialEnergySum += potentialEnergy;
            virialSumX += virialIx;
            virialSumY += virialIy;
            virialSumZ += virialIz;
            // for non-newton3 particles j and/or k will be considered in a separate calculation
            // if (newton3 and ownedStateJ) {
              // const auto virialJ = forceJ * j.getR();
              // _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
              // _aosThreadDataGlobals[threadnum].virialSum += virialJ;
            // }
            // if (newton3 and ownedStateK) {
            //   const auto virialK = forceK * k.getR();
            //   _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
            //   _aosThreadDataGlobals[threadnum].virialSum += virialK;
            // }
            // if constexpr (countFLOPs) {
            //   if (newton3) {
            //     ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsN3;
            //   } else {
            //     ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3;
            //   }
            // }
          }
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
    _aosThreadDataGlobals[threadnum]._timerTriple.stop();

  }

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
  const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
    const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {

    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;  // Note: This is not the potential energy but some fixed multiple of it.
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    SoAFloatPrecision fxacc = 0.;
    SoAFloatPrecision fyacc = 0.;
    SoAFloatPrecision fzacc = 0.;
    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    const auto ownedStateI = ownedStatePtr[indexFirst];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }
    constexpr size_t vecsize = 12;

    size_t joff = 0;
    size_t koff = 0;

    // if (neighborListSize >= vecsize) {
    //   alignas(64) std::array<SoAFloatPrecision, vecsize> xtmp, ytmp, ztmp, xArr, yArr, zArr, fxArr, fyArr, fzArr;
    //   alignas(64) std::array<autopas::OwnershipState, vecsize> ownedStateArr{};
    //
    //   for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
    //     xtmp[tmpj] = xptr[indexFirst];
    //     ytmp[tmpj] = yptr[indexFirst];
    //     ztmp[tmpj] = zptr[indexFirst];
    //   }
    // }

    for (size_t jNeighIndex = joff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;
      const auto ownedStateJ = ownedStatePtr[j];
      if (ownedStateJ == autopas::OwnershipState::dummy) {
        continue;
      }
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ) safelen(vecsize)
      for (size_t kNeighIndex = jNeighIndex + 1; kNeighIndex < neighborListSize; ++kNeighIndex) {
        size_t k = neighborList[kNeighIndex];
        // if (indexFirst == k) continue;
        const auto ownedStateK = ownedStatePtr[j];
        // if (ownedStateK == autopas::OwnershipState::dummy) {
        //   continue;
        // }

        const auto displacementIJx = xptr[j] - xptr[indexFirst];
        const auto displacementIJy = yptr[j] - yptr[indexFirst];
        const auto displacementIJz = zptr[j] - zptr[indexFirst];
        const auto displacementIJx2 = displacementIJx * displacementIJx;
        const auto displacementIJy2 = displacementIJy * displacementIJy;
        const auto displacementIJz2 = displacementIJz * displacementIJz;

        const auto displacementJKx = xptr[k] - xptr[j];
        const auto displacementJKy = yptr[k] - yptr[j];
        const auto displacementJKz = zptr[k] - zptr[j];
        const auto displacementJKx2 = displacementJKx * displacementJKx;
        const auto displacementJKy2 = displacementJKy * displacementJKy;
        const auto displacementJKz2 = displacementJKz * displacementJKz;

        const auto displacementKIx = xptr[indexFirst] - xptr[k];
        const auto displacementKIy = yptr[indexFirst] - yptr[k];
        const auto displacementKIz = zptr[indexFirst] - zptr[k];
        const auto displacementKIx2 = displacementKIx * displacementKIx;
        const auto displacementKIy2 = displacementKIy * displacementKIy;
        const auto displacementKIz2 = displacementKIz * displacementKIz;

        const auto distSquaredIJ = displacementIJx2 + displacementIJy2 + displacementIJz2;
        const auto distSquaredJK = displacementJKx2 + displacementJKy2 + displacementJKz2;
        const auto distSquaredKI = displacementKIx2 + displacementKIy2 + displacementKIz2;

        const bool mask = distSquaredIJ <= cutoffSquared and distSquaredJK <= cutoffSquared
          and distSquaredKI <= cutoffSquared  and indexFirst != k and ownedStateK != autopas::OwnershipState::dummy;

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
        const auto forceIDirectionIJx = displacementIJx * (fullATMGradientIJ + fullExpGradientIJ);
        const auto forceIDirectionIJy = displacementIJy * (fullATMGradientIJ + fullExpGradientIJ);
        const auto forceIDirectionIJz = displacementIJz * (fullATMGradientIJ + fullExpGradientIJ);

        const auto forceIDirecationKIx = displacementKIx * (fullATMGradientKI + fullExpGradientKI);
        const auto forceIDirecationKIy = displacementKIy * (fullATMGradientKI + fullExpGradientKI);
        const auto forceIDirecationKIz = displacementKIz * (fullATMGradientKI + fullExpGradientKI);

        const auto forceIx = (forceIDirectionIJx + forceIDirecationKIx) * mask * (-1.0);
        const auto forceIy = (forceIDirectionIJy + forceIDirecationKIy) * mask * (-1.0);
        const auto forceIz = (forceIDirectionIJz + forceIDirecationKIz) * mask * (-1.0);

        fxacc += forceIx;
        fyacc += forceIy;
        fzacc += forceIz;

        if constexpr (calculateGlobals) {
          // Add 3 * potential energy to every owned particle of the interaction.
          // Division to the correct value is handled in endTraversal().
          const double potentialEnergy = mask * (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

          // Virial is calculated as f_i * r_i
          // see Thompson et al.: https://doi.org/10.1063/1.3245303
          const auto virialIy = mask * forceIy * yptr[indexFirst];
          const auto virialIz = mask * forceIz * zptr[indexFirst];
          const auto virialIx = mask * forceIx * xptr[indexFirst];

          potentialEnergySum += potentialEnergy;
          virialSumX += virialIx;
          virialSumY += virialIy;
          virialSumZ += virialIz;
        }
      }
    }
    fxptr[indexFirst] += fxacc;
    fyptr[indexFirst] += fyacc;
    fzptr[indexFirst] += fzacc;
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
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
    double timerSingleSum = 0.;
    double timerPairSum = 0.;
    double timerTripleSum = 0.;
    if (calculateGlobals) {
      // Accumulate potential energy and virial values.
      for (const auto &data : _aosThreadDataGlobals) {
        _potentialEnergySum += data.potentialEnergySum;
        _virialSum += data.virialSum;
        timerSingleSum += data._timerSingle.getTotalTime();
        timerPairSum += data._timerPair.getTotalTime();
        timerTripleSum += data._timerTriple.getTotalTime();
      }

      // For each interaction, we added the full contribution for all three particles. Divide by 3 here, so that each
      // contribution is only counted once per triplet.
      _potentialEnergySum /= 3.;

      _postProcessed = true;

      // std::cout << "Time in AoSFunctor Total: " << timerSingleSum / 1000000 << "ms" << std::endl;
      // std::cout << "Time in AoSFunctor Globals: " << timerPairSum / 1000000 << "ms" << std::endl;
      // std::cout << "Time in AoSFunctor loops: " << timerTripleSum / 1000000 << "ms => " << (timerTripleSum / timerSingleSum) * 100 << "%" << std::endl;
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

    autopas::utils::Timer _timerSingle{};
    autopas::utils::Timer _timerPair{};
    autopas::utils::Timer _timerTriple{};

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
  // static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
  static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

  const double _cutoffSquared;

  // Parameters of the extended Axilrod-Teller potential for Krypton in Kelvin (K) and Angström (A)
  const double _nu = 1.61525e6;   // K*A^9
  const double _alpha = 1.378382;  // A^-1
  // Units: {K, K*A^-2, K*A^-4, K*A^-6, K*A^-8, K*A^-10}
  const std::array<double, 6> _constantsA = {-0.3081304e8, -0.3519442e8, 0.4928052e7, -0.2182411e6, 0.343088e4, 0.0};

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

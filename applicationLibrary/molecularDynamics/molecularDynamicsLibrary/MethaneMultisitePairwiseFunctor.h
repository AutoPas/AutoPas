/**
 * @file MethaneMultisitePairwiseFunctor.h
 * @date 04/11/2025
 * @author muehlhaeusser
 */

#pragma once

#include "MoleculeLJ.h"
#include "MultisiteMoleculeLJ.h"
#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

/**
 * A functor that implements the ab initio methane potential derived by Hellmann et al.
 * (https://doi.org/10.1063/1.2932103).
 *
 * @tparam Particle_T The type of particle.
 * @tparam useCorrection Adds the isotropic correction term
 * @tparam useNewton3 Switch for the functor to support newton3 on, off, or both. See FunctorN3Nodes for possible
 * values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether the auto-tuner should consider this functor.
 * @tparam countFLOPs counts FLOPs and hitrate.
 */
template <class Particle_T, bool useCorrection = true,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class MethaneMultisitePairwiseFunctor
    : public autopas::PairwiseFunctor<
          Particle_T, MethaneMultisitePairwiseFunctor<Particle_T, useCorrection, useNewton3, calculateGlobals,
                                                      relevantForTuning, countFLOPs>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Precision of SoA entries
   */
  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

  /**
   * cutoff^2
   */
  const double _cutoffSquared;

  /**
   * Sum of potential energy. Only calculated if calculateGlobals is true.
   */
  double _potentialEnergySum;

  /**
   * Sum of the virial. Only calculated if calculateGlobals is true.
   */
  std::array<double, 3> _virialSum;

  /**
   * Defines whether or whether not the global values are already processed
   */
  bool _postProcessed;

  ParticlePropertiesLibrary<double, size_t> *_PPLibrary;

 public:
  /**
   * Delete Default constructor
   */
  MethaneMultisitePairwiseFunctor() = delete;

 private:
  /**
   * Internal (actually used) constructor
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit MethaneMultisitePairwiseFunctor(SoAFloatPrecision cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle_T,
                                 MethaneMultisitePairwiseFunctor<Particle_T, useCorrection, useNewton3,
                                                                 calculateGlobals, relevantForTuning, countFLOPs>>(
            cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "FLOP counting is enabled but is not supported for multi-site functors yet!");
    }
  }

  constexpr std::array<double, 3> makeSitePosition(double x, double y, double z) {
    // Precomputed value of 1.0 / sqrt(3) for tetrahedral positions.
    constexpr double invSqrt3 = 0.5773502691896258;
    constexpr double factor = 1.099 * invSqrt3;
    return std::array{x * factor, y * factor, z * factor};
  }

  /**
   * List of relative unrotated Methane Site Positions.
   */
  const std::vector<std::array<double, 3>> _sitePositions{
      {0., 0., 0.},                        // C atom, site "C"
      makeSitePosition(0.88, 0.88, 0.88),  // Along C-H bond, site "H"
      makeSitePosition(0.88, -0.88, -0.88),
      makeSitePosition(-0.88, 0.88, -0.88),
      makeSitePosition(-0.88, -0.88, 0.88),
      makeSitePosition(-0.66, -0.66, -0.66),  // site "E"
      makeSitePosition(-0.66, 0.66, 0.66),
      makeSitePosition(0.66, -0.66, 0.66),
      makeSitePosition(0.66, 0.66, -0.66)};

  // Site IDs: C, H, H, H, H, E, E, E, E
  const std::vector<size_t> _siteIds{{0, 1, 1, 1, 1, 2, 2, 2, 2}};

  const std::array<double, 6> _paramsA{0.262373610e7,  0.265413949e7,  0.241399203e6,
                                       -0.271732286e6, -0.749715218e5, 0.123654939e6};

  const std::array<double, 6> _paramsAlpha{0.16878421e1, 0.28827219e1, 0.35917561e1,
                                           0.16490747e1, 0.20593086e1, 0.21451641e1};

  static constexpr double bCC = 0.168275675e1;
  static constexpr double bCH = 0.288261054e1;
  static constexpr double bHH = 0.384703188e1;
  static constexpr double bCE = 0.155011960e1;
  static constexpr double bHE = 0.266424603e1;
  static constexpr double bEE = 0.304993944e1;

  const std::array<double, 6> _paramsC6{0.112317356e7, -0.139633537e7, 0.294147230e6,
                                        0.127844394e7, 0.169329268e6,  -0.590727146e6};

  const std::array<double, 6> _paramsC8{-0.120939119e9, 0.385078060e8,  -0.264781786e7,
                                        0.174762764e7,  -0.810401688e7, 0.679543866e7};

  // Site charges {C, H, E}
  const std::array<double, 3> _charges{-0.379012e3, 0.94753e2, 0.};

  static constexpr double _bCorr = 0.177e1;
  const double _deltaC6 = 0.45347e5;
  const double _deltaC8 = 0.432463e6;

  // Compute base^power at compile time
  static constexpr double bToX(double base, size_t power) {
    double result = 1.0;
    for (size_t i = 0; i < power; ++i) {
      result *= base;
    }
    return result;
  }

  // Generate full table b^0..9 for all 6 b values
  static constexpr auto generateParameterBPowers() {
    return std::array<std::array<double, 10>, 6>{
        {{bToX(bCC, 0), bToX(bCC, 1), bToX(bCC, 2), bToX(bCC, 3), bToX(bCC, 4), bToX(bCC, 5), bToX(bCC, 6),
          bToX(bCC, 7), bToX(bCC, 8), bToX(bCC, 9)},
         {bToX(bCH, 0), bToX(bCH, 1), bToX(bCH, 2), bToX(bCH, 3), bToX(bCH, 4), bToX(bCH, 5), bToX(bCH, 6),
          bToX(bCH, 7), bToX(bCH, 8), bToX(bCH, 9)},
         {bToX(bHH, 0), bToX(bHH, 1), bToX(bHH, 2), bToX(bHH, 3), bToX(bHH, 4), bToX(bHH, 5), bToX(bHH, 6),
          bToX(bHH, 7), bToX(bHH, 8), bToX(bHH, 9)},
         {bToX(bCE, 0), bToX(bCE, 1), bToX(bCE, 2), bToX(bCE, 3), bToX(bCE, 4), bToX(bCE, 5), bToX(bCE, 6),
          bToX(bCE, 7), bToX(bCE, 8), bToX(bCE, 9)},
         {bToX(bHE, 0), bToX(bHE, 1), bToX(bHE, 2), bToX(bHE, 3), bToX(bHE, 4), bToX(bHE, 5), bToX(bHE, 6),
          bToX(bHE, 7), bToX(bHE, 8), bToX(bHE, 9)},
         {bToX(bEE, 0), bToX(bEE, 1), bToX(bEE, 2), bToX(bEE, 3), bToX(bEE, 4), bToX(bEE, 5), bToX(bEE, 6),
          bToX(bEE, 7), bToX(bEE, 8), bToX(bEE, 9)}}};
  }

  static constexpr auto generateParameterBCorrPowers() {
    return std::array<double, 10>{{bToX(_bCorr, 0), bToX(_bCorr, 1), bToX(_bCorr, 2), bToX(_bCorr, 3), bToX(_bCorr, 4),
                                   bToX(_bCorr, 5), bToX(_bCorr, 6), bToX(_bCorr, 7), bToX(_bCorr, 8),
                                   bToX(_bCorr, 9)}};
  }

  static constexpr auto _paramsB = generateParameterBPowers();
  static constexpr auto _bCorrPowers = generateParameterBCorrPowers();

 public:
  /**
   * @param cutoff
   */
  explicit MethaneMultisitePairwiseFunctor(double cutoff) : MethaneMultisitePairwiseFunctor(cutoff, nullptr) {}

  explicit MethaneMultisitePairwiseFunctor(double cutoff,
                                           ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : MethaneMultisitePairwiseFunctor(cutoff, nullptr) {
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "MethaneMultisitePairwiseFunctor"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  /**
   * Functor for arrays of structures (AoS).
   *
   * @param particleA Particle A
   * @param particleB Particle B
   * @param newton3 Flag for if newton3 is used.
   */
  void AoSFunctor(Particle_T &particleA, Particle_T &particleB, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (particleA.isDummy() or particleB.isDummy()) {
      return;
    }

    // Don't calculate force if particleB outside cutoff of particleA
    const auto displacementCoM = autopas::utils::ArrayMath::sub(particleA.getR(), particleB.getR());
    const auto distanceSquaredCoM = autopas::utils::ArrayMath::dot(displacementCoM, displacementCoM);

    if (distanceSquaredCoM > _cutoffSquared) {
      return;
    }

    const size_t numSitesA = _PPLibrary ? _PPLibrary->getNumSites(particleA.getTypeId()) : 9;
    const size_t numSitesB = _PPLibrary ? _PPLibrary->getNumSites(particleB.getTypeId()) : 9;

    // get siteIds
    const auto siteIdsA = _PPLibrary ? _PPLibrary->getSiteTypes(particleA.getTypeId()) : _siteIds;
    const auto siteIdsB = _PPLibrary ? _PPLibrary->getSiteTypes(particleB.getTypeId()) : _siteIds;

    // get unrotated relative site positions
    const std::vector<std::array<double, 3>> unrotatedSitePositionsA =
        _PPLibrary ? _PPLibrary->getSitePositions(particleA.getTypeId()) : _sitePositions;
    const std::vector<std::array<double, 3>> unrotatedSitePositionsB =
        _PPLibrary ? _PPLibrary->getSitePositions(particleB.getTypeId()) : _sitePositions;

    // calculate correctly rotated relative site positions
    const auto rotatedSitePositionsA =
        autopas::utils::quaternion::rotateVectorOfPositions(particleA.getQuaternion(), unrotatedSitePositionsA);
    const auto rotatedSitePositionsB =
        autopas::utils::quaternion::rotateVectorOfPositions(particleB.getQuaternion(), unrotatedSitePositionsB);

    for (int i = 0; i < numSitesA; i++) {
      const auto siteTypeI = siteIdsA[i];
      for (int j = 0; j < numSitesB; j++) {
        const auto siteTypeJ = siteIdsB[j];
        // Map site type combination to 0-5 index
        const int paramIndex = (siteTypeI >= siteTypeJ) ? siteTypeI * (siteTypeI + 1) / 2 + siteTypeJ
                                                        : siteTypeJ * (siteTypeJ + 1) / 2 + siteTypeI;
        const auto alpha = _paramsAlpha[paramIndex];
        const auto A = _paramsA[paramIndex];
        const auto bPows = _paramsB[paramIndex];
        const auto C6 = _paramsC6[paramIndex];
        const auto C8 = _paramsC8[paramIndex];

        const auto displacement = autopas::utils::ArrayMath::add(
            autopas::utils::ArrayMath::sub(displacementCoM, rotatedSitePositionsB[j]), rotatedSitePositionsA[i]);

        const auto distSquared = autopas::utils::ArrayMath::dot(displacement, displacement);
        const auto dist = std::sqrt(distSquared);
        const auto distInv = 1. / dist;
        const auto distInv2 = distInv * distInv;
        const auto distInv3 = distInv2 * distInv;
        const auto distInv4 = distInv3 * distInv;
        const auto distInv5 = distInv4 * distInv;
        const auto distInv6 = distInv5 * distInv;
        const auto distInv7 = distInv6 * distInv;
        const auto distInv8 = distInv7 * distInv;
        const auto distInv9 = distInv8 * distInv;

        const auto expTerm = alpha * A * std::exp(-alpha * dist);

        const auto expBR = std::exp(-bPows[1] * dist);

        const auto tangToenniesSum6 = 6. * distInv7 + 6. * bPows[1] * distInv6 + 3. * bPows[2] * distInv5 +
                                      bPows[3] * distInv4 + bPows[4] * distInv3 / 4. + bPows[5] * distInv2 / 20. +
                                      bPows[6] * distInv / 120. + bPows[7] / 720.;

        const auto tangToenniesSum8 = 8. * distInv9 + 8. * bPows[1] * distInv8 + 4. * bPows[2] * distInv7 +
                                      4. * bPows[3] * distInv6 / 3. + bPows[4] * distInv5 / 3. +
                                      bPows[5] * distInv4 / 15. + bPows[6] * distInv3 / 90. +
                                      bPows[7] * distInv2 / 630. + bPows[8] * distInv / 5040. + bPows[9] / 40320.;

        const auto tangToenniesTerm6 = C6 * (-expBR * tangToenniesSum6 + 6. * distInv7);
        const auto tangToenniesTerm8 = C8 * (-expBR * tangToenniesSum8 + 8. * distInv9);

        const auto chargeTerm = _charges[siteTypeI] * _charges[siteTypeJ] * distInv2;

        auto forceFactor = (-expTerm + tangToenniesTerm6 + tangToenniesTerm8 - chargeTerm) * distInv;

        // Add correction term once for the C-C interaction (as they are located at the molecule centers)
        if (useCorrection and (i == 0 and j == 0)) {
          constexpr auto bCorrPows = _bCorrPowers;
          const auto expBRCorr = std::exp(-bCorrPows[1] * dist);

          const auto tangToenniesCorrectionSum6 = 6. * distInv7 + 6. * bCorrPows[1] * distInv6 +
                                                  3. * bCorrPows[2] * distInv5 + bCorrPows[3] * distInv4 +
                                                  bCorrPows[4] * distInv3 / 4. + bCorrPows[5] * distInv2 / 20. +
                                                  bCorrPows[6] * distInv / 120. + bCorrPows[7] / 720.;

          const auto tangToenniesCorrectionSum8 = 8. * distInv9 + 8. * bCorrPows[1] * distInv8 +
                                                  4. * bCorrPows[2] * distInv7 + 4. * bCorrPows[3] * distInv6 / 3. +
                                                  bCorrPows[4] * distInv5 / 3. + bCorrPows[5] * distInv4 / 15. +
                                                  bCorrPows[6] * distInv3 / 90. + bCorrPows[7] * distInv2 / 630. +
                                                  bCorrPows[8] * distInv / 5040. + bCorrPows[9] / 40320.;

          const auto tangToenniesTermCorrection6 = _deltaC6 * (-expBRCorr * tangToenniesCorrectionSum6 + 6. * distInv7);
          const auto tangToenniesTermCorrection8 = _deltaC8 * (-expBRCorr * tangToenniesCorrectionSum8 + 8. * distInv9);

          forceFactor += (tangToenniesTermCorrection6 + tangToenniesTermCorrection8) * distInv;
        }

        const auto force = autopas::utils::ArrayMath::mulScalar(displacement, forceFactor);

        // Add force on site to net force
        particleA.addF(force);
        if (newton3) {
          particleB.subF(force);
        }

        // Add torque applied by force
        particleA.addTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsA[i], force));
        if (newton3) {
          particleB.subTorque(autopas::utils::ArrayMath::cross(rotatedSitePositionsB[j], force));
        }

        if (calculateGlobals) {
          // We always add the full contribution for each owned particle and divide the sums by 2 in endTraversal().
          // Potential energy has an additional factor of 6, which is also handled in endTraversal().
          const auto innerSum6 = distInv6 + bPows[1] * distInv5 + bPows[2] * distInv4 / 2. + bPows[3] * distInv3 / 6. +
                                 bPows[4] * distInv2 / 24. + bPows[5] * distInv / 120. + bPows[6] / 720.;
          const auto innerSum8 = innerSum6 * distInv2 + bPows[7] * distInv / 5040. + bPows[8] / 40320.;
          auto potentialEnergy = A * std::exp(-alpha * dist) - C6 * (distInv6 - expBR * innerSum6) -
                                 C8 * (distInv8 - expBR * innerSum8) +
                                 _charges[siteTypeI] * _charges[siteTypeJ] * distInv;
          if (useCorrection and (i == 0 and j == 0)) {
            const auto bCorrPows = _bCorrPowers;
            const auto expBRCorr = std::exp(-bCorrPows[1] * dist);

            const auto innerCorrectionSum6 = distInv6 + bCorrPows[1] * distInv5 + bCorrPows[2] * distInv4 / 2. +
                                             bCorrPows[3] * distInv3 / 6. + bCorrPows[4] * distInv2 / 24. +
                                             bCorrPows[5] * distInv / 120. + bCorrPows[6] / 720.;

            const auto innerCorrectionSum8 =
                innerCorrectionSum6 * distInv2 + bCorrPows[7] * distInv / 5040. + bCorrPows[8] / 40320.;

            potentialEnergy -= _deltaC6 * (distInv6 - expBRCorr * innerCorrectionSum6) +
                               _deltaC8 * (distInv8 - expBRCorr * innerCorrectionSum8);
          }
          const auto virial = displacement * force;

          const auto threadNum = autopas::autopas_get_thread_num();

          if (particleA.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy;
            _aosThreadData[threadNum].virialSum += virial;
          }
          // for non-newton3 the second particle will be considered in a separate calculation
          if (newton3 and particleB.isOwned()) {
            _aosThreadData[threadNum].potentialEnergySum += potentialEnergy;
            _aosThreadData[threadNum].virialSum += virial;
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversing of the soa, however, it still needs to know about newton3
   * to use it correctly for the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.size() == 0) return;
  }

  /*
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

  // clang-format off
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  /**
   * Sets the molecule properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   * @param sitePositions vector of 3D relative unrotated untranslated site positions
   */
  void setParticleProperties(std::vector<std::array<SoAFloatPrecision, 3>> sitePositions) {
    _sitePositions = sitePositions;
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 16>{
        Particle_T::AttributeNames::id,          Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,        Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::forceX,      Particle_T::AttributeNames::forceY,
        Particle_T::AttributeNames::forceZ,      Particle_T::AttributeNames::quaternion0,
        Particle_T::AttributeNames::quaternion1, Particle_T::AttributeNames::quaternion2,
        Particle_T::AttributeNames::quaternion3, Particle_T::AttributeNames::torqueX,
        Particle_T::AttributeNames::torqueY,     Particle_T::AttributeNames::torqueZ,
        Particle_T::AttributeNames::typeId,      Particle_T::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle_T::AttributeNames, 16>{
        Particle_T::AttributeNames::id,          Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,        Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::forceX,      Particle_T::AttributeNames::forceY,
        Particle_T::AttributeNames::forceZ,      Particle_T::AttributeNames::quaternion0,
        Particle_T::AttributeNames::quaternion1, Particle_T::AttributeNames::quaternion2,
        Particle_T::AttributeNames::quaternion3, Particle_T::AttributeNames::torqueX,
        Particle_T::AttributeNames::torqueY,     Particle_T::AttributeNames::torqueZ,
        Particle_T::AttributeNames::typeId,      Particle_T::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 6>{
        Particle_T::AttributeNames::forceX,  Particle_T::AttributeNames::forceY,  Particle_T::AttributeNames::forceZ,
        Particle_T::AttributeNames::torqueX, Particle_T::AttributeNames::torqueY, Particle_T::AttributeNames::torqueZ};
  }

  /**
   * @return useMixing
   */
  constexpr static bool getMixing() { return false; }

  /**
   * Get the number of flops used per kernel call - i.e. number of flops to calculate kernel *given* the two particles
   * lie within the cutoff (i.e. distance^2 / cutoff has been already been calculated).
   * Note: there is currently a large difference between AoS & SoA number of flops. This function returns the AoS
   * number of flops.
   * @param molAType molecule A's type id
   * @param molBType molecule B's type id
   * @param newton3 true if newton3 optimizations enabled
   * @return Number of FLOPs
   */
  unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, bool newton3) {
    // Site-to-site displacement: 6 (3 in the SoA case, but this requires O(N) precomputing site positions)
    // Site-to-site distance squared: 4
    // Compute scale: 9
    // Apply scale to force: With newton3: 6, Without: 3
    // Apply scale to torque: With newton3 18, Without: 9 (0 in SoA case, with O(N) post computing)
    // Site-to-site total: With newton3: 33, Without: 26
    // (SoA total: With N3L: 22, Without N3L: 19)
    // Above multiplied by number sites of i * number sites of j
    const unsigned long siteToSiteFlops = newton3 ? 33ul : 26ul;
    return 0.0;
  }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() final {
    _potentialEnergySum = 0;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); i++) {
      _aosThreadData[i].setZero();
    }
  }

  /**
   * Postprocesses global values, e.g. potential energy & virial
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
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
        _virialSum += _aosThreadData[i].virialSum;
      }

      // For each interaction, we added the full contribution for both particles. Divide by 2 here, so that each
      // contribution is only counted once per pair.
      _potentialEnergySum *= 0.5;
      _virialSum *= 0.5;

      _postProcessed = true;

      AutoPasLog(DEBUG, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(DEBUG, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential energy.
   *
   * @return the potential energy
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

 private:
  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   * @tparam newton3 flag for if newton's third law is used
   * @param soaA structure of arrays A
   * @param soaB structure of arrays B
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soaA, autopas::SoAView<SoAArraysType> soaB) {
    if (soaA.size() == 0 || soaB.size() == 0) return;
  }

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexPrime,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {}

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

  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size (should be multiple of 64)");

  /**
   * Thread buffer for AoS
   */
  std::vector<AoSThreadData> _aosThreadData;
};
}  // namespace mdLib
/**
 * @file AxilrodTellerMutoFunctorHWY.h
 * @author D. Martin
 * @date 13/11/25
 */

#pragma once

#include <hwy/highway.h>

#include "HighwayDefs.h"
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
 * Note: This functor only supports aligned SoAs. A version that also supports unaligned SoAs can be found here:
 * https://github.com/AutoPas/AutoPas/commit/2ccd4f57333e219f237e64e21b00b7a34ecca5fe
 * @tparam Particle_T The type of particle.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 */
template <class Particle_T, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool countFLOPs = false>
class AxilrodTellerMutoFunctorHWY
    : public autopas::TriwiseFunctor<
          Particle_T, AxilrodTellerMutoFunctorHWY<Particle_T, useMixing, useNewton3, calculateGlobals, countFLOPs>> {
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
  AxilrodTellerMutoFunctorHWY() = delete;

  /**
   * Constructor for Functor with mixing enabled/disabled. When using this functor without particlePropertiesLibrary it
   * is necessary to call setParticleProperties() to set internal constants.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit AxilrodTellerMutoFunctorHWY(double cutoff,
                                       std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>>
                                           particlePropertiesLibrary = std::nullopt)
      : autopas::TriwiseFunctor<Particle_T, AxilrodTellerMutoFunctorHWY>(cutoff),
        _cutoffSquared{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadDataGlobals(),
        _postProcessed{false},
        _PPLibrary{particlePropertiesLibrary} {
    const size_t numMaxThreads = autopas::autopas_get_max_threads();
    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals.resize(numMaxThreads);
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs.resize(numMaxThreads);
      AutoPasLog(DEBUG, "Using AxilrodTellerMutoFunctorHWY with countFLOPs is not tested for SoA datalayout.");
    }

    if constexpr (useMixing) {
      if (not _PPLibrary.has_value()) {
        throw std::runtime_error("Mixing is enabled but no ParticlePropertiesLibrary was provided!");
      }
    } else {
      if (_PPLibrary.has_value()) {
        throw std::runtime_error("Mixing is disabled but a ParticlePropertiesLibrary was provided!");
      }
    }

    precomputeBuffers1.resize(numMaxThreads);
    precomputeBuffers2.resize(numMaxThreads);
  }

  std::string getName() final { return "AxilrodTellerMutoFunctorHWY"; }

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
      nu = _PPLibrary->get().getMixingNu(i.getTypeId(), j.getTypeId(), k.getTypeId());
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

    // check if the position, force and ownership SoAs are aligned to enable aligned vector loads
    const bool alignedVectors = areAllSoAsAligned(soa);
    if (not alignedVectors) {
      throw autopas::utils::ExceptionHandler::AutoPasException("AxilrodTellerMutoFunctorHWY only allows aligned SoAs.");
    }

    SoAFunctorSingleImpl(soa);
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final {
    if (soa1.size() == 0 || soa2.size() == 0) return;

    // check if the position, force and ownership SoAs are aligned to enable aligned vector loads
    const bool alignedVectors = areAllSoAsAligned(soa1, soa2);
    if (not alignedVectors) {
      throw autopas::utils::ExceptionHandler::AutoPasException("AxilrodTellerMutoFunctorHWY only allows aligned SoAs.");
    }

    if (newton3) {
      SoAFunctorPairImpl</*newton3*/ true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl</*newton3*/ false>(soa1, soa2);
    }
  }

  void SoAFunctorTriple(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                        autopas::SoAView<SoAArraysType> soa3, bool newton3) final {
    if (soa1.size() == 0 || soa2.size() == 0 || soa3.size() == 0) return;

    // check if the position, force and ownership SoAs are aligned to enable aligned vector loads
    const bool alignedVectors = areAllSoAsAligned(soa1, soa2, soa3);
    if (not alignedVectors) {
      throw autopas::utils::ExceptionHandler::AutoPasException("AxilrodTellerMutoFunctorHWY only allows aligned SoAs.");
    }

    if (newton3) {
      SoAFunctorTripleImpl</*newton3*/ true>(soa1, soa2, soa3);
    } else {
      SoAFunctorTripleImpl</*newton3*/ false>(soa1, soa2, soa3);
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
  struct PrecomputeBuffer;

  void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typePtr = soa.template begin<Particle_T::AttributeNames::typeId>();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    VectorDouble potentialEnergySum =
        _zeroDoubleVec;  // Note: This is not the potential energy but some fixed multiple of it.
    VectorDouble virialSumX = _zeroDoubleVec;
    VectorDouble virialSumY = _zeroDoubleVec;
    VectorDouble virialSumZ = _zeroDoubleVec;

    const size_t soaSize = soa.size();
    const size_t numTriplets =
        countFLOPs ? soaSize * (soaSize - 1) * (soaSize - 2) / 6 : 0;  // Only needed if counting FLOPs
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;
    const SoAFloatPrecision const_nu = _nu;

    // Precompute distances between particles in the soa. We use a row-major packed lower‐triangle data structure to
    // save memory. Only distances in the form of i->j are stored, but not j->i. If j->i is required, it can be
    // calculated by negating i->j.
    auto &intraSoAPairDists = precomputeBuffers1[threadnum];
    intraSoAPairDists.template init</*LowerTriangle*/ true>(soaSize, soaSize);
    // soa1 <-> soa1
    intraSoAPairDists.template fillHighway</*LowerTriangle*/ true>(xPtr, yPtr, zPtr, xPtr, yPtr, zPtr, cutoffSquared);

    if constexpr (countFLOPs) {
      numDistanceCalculationSum += (soaSize * soaSize - soaSize) / 2;
    }

    // We iterate in reverse order over i and j to enable aligned loads in the k-loop.
    for (size_t i = soaSize - 1; i > 1; --i) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }
      VectorDouble fXAccI = _zeroDoubleVec;
      VectorDouble fYAccI = _zeroDoubleVec;
      VectorDouble fZAccI = _zeroDoubleVec;

      for (size_t j = i - 1; j > 0; --j) {
        const auto ownedStateJ = ownedStatePtr[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        // Load precomputed distances for particle i and j.
        auto [distXIJ, distYIJ, distZIJ, distSquaredIJ, invR5IJ] =
            intraSoAPairDists.template get</*LowerTriangle*/ true>(i, j);

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        // prepare vector registers
        VectorDouble fXAccJ = _zeroDoubleVec;
        VectorDouble fYAccJ = _zeroDoubleVec;
        VectorDouble fZAccJ = _zeroDoubleVec;

        const auto distXIJVec = highway::Set(tag_double, distXIJ);
        const auto distYIJVec = highway::Set(tag_double, distYIJ);
        const auto distZIJVec = highway::Set(tag_double, distZIJ);
        const auto distSquaredIJVec = highway::Set(tag_double, distSquaredIJ);
        const auto invR5IJVec = highway::Set(tag_double, invR5IJ);

        const size_t blockEnd = j;
        for (unsigned int k = 0; k < blockEnd; k += _vecLengthDouble) {
          const auto numLanesToProcess = std::min(_vecLengthDouble, blockEnd - k);

          handleKLoopBody</*LowerPackedTrianglePB1*/ true, /*LowerPackedTrianglePB2*/ true, /*newton3*/ true,
                          /*newton3Kernel*/ true>(
              i, j, k, const_nu, intraSoAPairDists, intraSoAPairDists, ownedStatePtr, typePtr, typePtr, typePtr,
              distXIJVec, distYIJVec, distZIJVec, distSquaredIJVec, invR5IJVec, fXAccI, fYAccI, fZAccI, fXAccJ, fYAccJ,
              fZAccJ, ownedStateI, ownedStateJ, fxPtr, fyPtr, fzPtr, virialSumX, virialSumY, virialSumZ,
              potentialEnergySum, numKernelCallsN3Sum, numGlobalCalcsN3Sum, numKernelCallsNoN3Sum,
              numGlobalCalcsNoN3Sum, numLanesToProcess);
        }

        // Reduce force on particle J.
        fxPtr[j] += highway::ReduceSum(tag_double, fXAccJ);
        fyPtr[j] += highway::ReduceSum(tag_double, fYAccJ);
        fzPtr[j] += highway::ReduceSum(tag_double, fZAccJ);
      }
      // Reduce force on particle I.
      fxPtr[i] += highway::ReduceSum(tag_double, fXAccI);
      fyPtr[i] += highway::ReduceSum(tag_double, fYAccI);
      fzPtr[i] += highway::ReduceSum(tag_double, fZAccI);
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numTripletsCount += numTriplets;
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;  // Always N3 in Single SoAFunctor
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += highway::ReduceSum(tag_double, potentialEnergySum);
      _aosThreadDataGlobals[threadnum].virialSum[0] += highway::ReduceSum(tag_double, virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += highway::ReduceSum(tag_double, virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += highway::ReduceSum(tag_double, virialSumZ);
    }
  }

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

    auto *const __restrict fxPtr1 = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr1 = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr1 = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fxPtr2 = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr2 = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr2 = soa2.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle_T::AttributeNames::typeId>();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    VectorDouble potentialEnergySum =
        _zeroDoubleVec;  // Note: This is not the potential energy but some fixed multiple of it.
    VectorDouble virialSumX = _zeroDoubleVec;
    VectorDouble virialSumY = _zeroDoubleVec;
    VectorDouble virialSumZ = _zeroDoubleVec;

    size_t soa1Size = soa1.size();
    size_t soa2Size = soa2.size();
    const size_t numTriplets = countFLOPs ? (soa1Size * soa2Size * (soa1Size + soa2Size - 2)) / 2 : 0;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;
    const SoAFloatPrecision const_nu = _nu;

    // Precompute distances between particles in the SoAs. soa2 <-> soa2 is a lower packed triangle matrix, soa1 <->
    // soa2 is a full matrix.
    auto &intraSoA2PairDists = precomputeBuffers1[threadnum];
    auto &interSoA1SoA2PairDists = precomputeBuffers2[threadnum];
    intraSoA2PairDists.template init</*LowerTriangle*/ true>(soa2Size, soa2Size);
    interSoA1SoA2PairDists.template init</*LowerTriangle*/ false>(soa1Size, soa2Size);

    // soa2 <-> soa2
    intraSoA2PairDists.template fillHighway</*LowerTriangle*/ true>(xptr2, yptr2, zptr2, xptr2, yptr2, zptr2,
                                                                    cutoffSquared);
    // soa1 <-> soa2
    interSoA1SoA2PairDists.template fillHighway</*LowerTriangle*/ false>(xptr1, yptr1, zptr1, xptr2, yptr2, zptr2,
                                                                         cutoffSquared);

    if constexpr (countFLOPs) {
      const size_t packedSizeSoA2 = (soa2Size * soa2Size - soa2Size) / 2;
      numDistanceCalculationSum += soa1Size * soa2Size + packedSizeSoA2;
    }

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision xi = xptr1[i];
      const SoAFloatPrecision yi = yptr1[i];
      const SoAFloatPrecision zi = zptr1[i];

      VectorDouble fXAccI = _zeroDoubleVec;
      VectorDouble fYAccI = _zeroDoubleVec;
      VectorDouble fZAccI = _zeroDoubleVec;

      // CASE: Particle i is in soa1, j and k are both in soa2
      // We iterate in reverse order over j to enable aligned stores in the k-loop.
      for (unsigned int j = soa2Size - 1; j > 0; --j) {
        const auto ownedStateJ = ownedStatePtr2[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        // Load precomputed distances for particle i and j.
        auto [distXIJ, distYIJ, distZIJ, distSquaredIJ, invR5IJ] =
            interSoA1SoA2PairDists.template get</*LowerTriangle*/ false>(i, j);

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        VectorDouble fXAccJ = _zeroDoubleVec;
        VectorDouble fYAccJ = _zeroDoubleVec;
        VectorDouble fZAccJ = _zeroDoubleVec;

        // prepare vector registers
        const auto distXIJVec = highway::Set(tag_double, distXIJ);
        const auto distYIJVec = highway::Set(tag_double, distYIJ);
        const auto distZIJVec = highway::Set(tag_double, distZIJ);
        const auto distSquaredIJVec = highway::Set(tag_double, distSquaredIJ);
        const auto invR5IJVec = highway::Set(tag_double, invR5IJ);

        const size_t blockEnd = j;
        for (unsigned int k = 0; k < blockEnd; k += _vecLengthDouble) {
          const auto numLanesToProcess = std::min(_vecLengthDouble, static_cast<size_t>(blockEnd - k));

          handleKLoopBody</*LowerPackedTrianglePB1*/ true, /*LowerPackedTrianglePB2*/ false, /*newton3*/ newton3,
                          /*newton3Kernel*/ newton3>(
              i, j, k, const_nu, intraSoA2PairDists, interSoA1SoA2PairDists, ownedStatePtr2, typeptr1, typeptr2,
              typeptr2, distXIJVec, distYIJVec, distZIJVec, distSquaredIJVec, invR5IJVec, fXAccI, fYAccI, fZAccI,
              fXAccJ, fYAccJ, fZAccJ, ownedStateI, ownedStateJ, fxPtr2, fyPtr2, fzPtr2, virialSumX, virialSumY,
              virialSumZ, potentialEnergySum, numKernelCallsN3Sum, numGlobalCalcsN3Sum, numKernelCallsNoN3Sum,
              numGlobalCalcsNoN3Sum, numLanesToProcess);
        }
        // Reduce force on particle j.
        fxPtr2[j] += highway::ReduceSum(tag_double, fXAccJ);
        fyPtr2[j] += highway::ReduceSum(tag_double, fYAccJ);
        fzPtr2[j] += highway::ReduceSum(tag_double, fZAccJ);
      }

      // CASE: Particle i and j are both in soa1, k is in soa2
      for (unsigned int j = i + 1; j < soa1.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr1[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        const SoAFloatPrecision xj = xptr1[j];
        const SoAFloatPrecision yj = yptr1[j];
        const SoAFloatPrecision zj = zptr1[j];

        const SoAFloatPrecision distXIJ = xj - xi;
        const SoAFloatPrecision distYIJ = yj - yi;
        const SoAFloatPrecision distZIJ = zj - zi;
        const SoAFloatPrecision distSquaredIJ = distXIJ * distXIJ + distYIJ * distYIJ + distZIJ * distZIJ;

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        const SoAFloatPrecision distIJ = std::sqrt(distSquaredIJ);
        const SoAFloatPrecision r5 = distSquaredIJ * distSquaredIJ * distIJ;
        const SoAFloatPrecision invR5IJ = 1.0 / r5;

        VectorDouble fXAccJ = _zeroDoubleVec;
        VectorDouble fYAccJ = _zeroDoubleVec;
        VectorDouble fZAccJ = _zeroDoubleVec;

        // prepare vector registers
        const auto distXIJVec = highway::Set(tag_double, distXIJ);
        const auto distYIJVec = highway::Set(tag_double, distYIJ);
        const auto distZIJVec = highway::Set(tag_double, distZIJ);
        const auto distSquaredIJVec = highway::Set(tag_double, distSquaredIJ);
        const auto invR5IJVec = highway::Set(tag_double, invR5IJ);

        unsigned int blockEnd = soa2.size();
        for (unsigned int k = 0; k < blockEnd; k += _vecLengthDouble) {
          const auto numLanesToProcess = std::min(_vecLengthDouble, static_cast<size_t>(blockEnd - k));

          handleKLoopBody</*LowerPackedTrianglePB1*/ false, /*LowerPackedTrianglePB2*/ false, /*newton3*/ newton3,
                          /*newton3Kernel*/ true>(
              i, j, k, const_nu, interSoA1SoA2PairDists, interSoA1SoA2PairDists, ownedStatePtr2, typeptr1, typeptr1,
              typeptr2, distXIJVec, distYIJVec, distZIJVec, distSquaredIJVec, invR5IJVec, fXAccI, fYAccI, fZAccI,
              fXAccJ, fYAccJ, fZAccJ, ownedStateI, ownedStateJ, fxPtr2, fyPtr2, fzPtr2, virialSumX, virialSumY,
              virialSumZ, potentialEnergySum, numKernelCallsN3Sum, numGlobalCalcsN3Sum, numKernelCallsNoN3Sum,
              numGlobalCalcsNoN3Sum, numLanesToProcess);
        }

        // Reduce force on particle j.
        fxPtr1[j] += highway::ReduceSum(tag_double, fXAccJ);
        fyPtr1[j] += highway::ReduceSum(tag_double, fYAccJ);
        fzPtr1[j] += highway::ReduceSum(tag_double, fZAccJ);
      }
      // Reduce force on particle i.
      fxPtr1[i] += highway::ReduceSum(tag_double, fXAccI);
      fyPtr1[i] += highway::ReduceSum(tag_double, fYAccI);
      fzPtr1[i] += highway::ReduceSum(tag_double, fZAccI);
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numTripletsCount += numTriplets;
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += highway::ReduceSum(tag_double, potentialEnergySum);

      _aosThreadDataGlobals[threadnum].virialSum[0] += highway::ReduceSum(tag_double, virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += highway::ReduceSum(tag_double, virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += highway::ReduceSum(tag_double, virialSumZ);
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

    auto *const __restrict fxPtr1 = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr1 = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr1 = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fxPtr2 = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr2 = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr2 = soa2.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fyPtr3 = soa3.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fxPtr3 = soa3.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fzPtr3 = soa3.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr3 = soa3.template begin<Particle_T::AttributeNames::typeId>();

    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    VectorDouble potentialEnergySum =
        _zeroDoubleVec;  // Note: This is not the potential energy but some fixed multiple of it.
    VectorDouble virialSumX = _zeroDoubleVec;
    VectorDouble virialSumY = _zeroDoubleVec;
    VectorDouble virialSumZ = _zeroDoubleVec;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;
    const SoAFloatPrecision const_nu = _nu;
    const auto soa1Size = soa1.size();
    const auto soa2Size = soa2.size();
    const auto soa3Size = soa3.size();

    const size_t numTriplets = countFLOPs ? soa1Size * soa2Size * soa3Size : 0;

    // Precompute distances between particles in the SoAs.
    auto &interSoA1SoA3PairDists = precomputeBuffers1[threadnum];
    auto &interSoA2SoA3PairDists = precomputeBuffers2[threadnum];

    interSoA1SoA3PairDists.template init</*LowerTriangle*/ false>(soa1Size, soa3Size);
    interSoA2SoA3PairDists.template init</*LowerTriangle*/ false>(soa2Size, soa3Size);

    // soa1 <-> soa3
    interSoA1SoA3PairDists.template fillHighway</*LowerTriangle*/ false>(xptr1, yptr1, zptr1, xptr3, yptr3, zptr3,
                                                                         cutoffSquared);
    // soa2 <-> soa3
    interSoA2SoA3PairDists.template fillHighway</*LowerTriangle*/ false>(xptr2, yptr2, zptr2, xptr3, yptr3, zptr3,
                                                                         cutoffSquared);

    if constexpr (countFLOPs) {
      numDistanceCalculationSum += soa1Size * soa3Size + soa2Size * soa3Size;
    }

    for (unsigned int i = 0; i < soa1Size; ++i) {
      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision xi = xptr1[i];
      const SoAFloatPrecision yi = yptr1[i];
      const SoAFloatPrecision zi = zptr1[i];

      VectorDouble fXAccI = _zeroDoubleVec;
      VectorDouble fYAccI = _zeroDoubleVec;
      VectorDouble fZAccI = _zeroDoubleVec;

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

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        const SoAFloatPrecision distIJ = std::sqrt(distSquaredIJ);
        const SoAFloatPrecision r5 = distSquaredIJ * distSquaredIJ * distIJ;
        const SoAFloatPrecision invR5IJ = 1.0 / r5;

        VectorDouble fXAccJ = _zeroDoubleVec;
        VectorDouble fYAccJ = _zeroDoubleVec;
        VectorDouble fZAccJ = _zeroDoubleVec;

        // prepare vector registers
        const auto distXIJVec = highway::Set(tag_double, distXIJ);
        const auto distYIJVec = highway::Set(tag_double, distYIJ);
        const auto distZIJVec = highway::Set(tag_double, distZIJ);
        const auto distSquaredIJVec = highway::Set(tag_double, distSquaredIJ);
        const auto invR5IJVec = highway::Set(tag_double, invR5IJ);

        unsigned int blockEnd = soa3.size();
        for (unsigned int k = 0; k < blockEnd; k += _vecLengthDouble) {
          const auto numLanesToProcess = std::min(_vecLengthDouble, static_cast<size_t>(blockEnd - k));

          handleKLoopBody</*LowerPackedTrianglePB1*/ false, /*LowerPackedTrianglePB2*/ false, /*newton3*/ newton3,
                          /*newton3Kernel*/ newton3>(
              i, j, k, const_nu, interSoA2SoA3PairDists, interSoA1SoA3PairDists, ownedStatePtr3, typeptr1, typeptr2,
              typeptr3, distXIJVec, distYIJVec, distZIJVec, distSquaredIJVec, invR5IJVec, fXAccI, fYAccI, fZAccI,
              fXAccJ, fYAccJ, fZAccJ, ownedStateI, ownedStateJ, fxPtr3, fyPtr3, fzPtr3, virialSumX, virialSumY,
              virialSumZ, potentialEnergySum, numKernelCallsN3Sum, numGlobalCalcsN3Sum, numKernelCallsNoN3Sum,
              numGlobalCalcsNoN3Sum, numLanesToProcess);
        }
        // Reduce force on particle j.
        fxPtr2[j] += highway::ReduceSum(tag_double, fXAccJ);
        fyPtr2[j] += highway::ReduceSum(tag_double, fYAccJ);
        fzPtr2[j] += highway::ReduceSum(tag_double, fZAccJ);
      }
      // Reduce force on particle i.
      fxPtr1[i] += highway::ReduceSum(tag_double, fXAccI);
      fyPtr1[i] += highway::ReduceSum(tag_double, fYAccI);
      fzPtr1[i] += highway::ReduceSum(tag_double, fZAccI);
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numTripletsCount += numTriplets;
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += highway::ReduceSum(tag_double, potentialEnergySum);

      _aosThreadDataGlobals[threadnum].virialSum[0] += highway::ReduceSum(tag_double, virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += highway::ReduceSum(tag_double, virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += highway::ReduceSum(tag_double, virialSumZ);
    }
  }

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    autopas::utils::ExceptionHandler::exception(
        "AxilrodTellerMutoFunctorHWY::SoAFunctorVerletImpl() is not implemented.");
  }

  /**
   * This function vectorizes the loop over k. Precalculated distances between particles are loaded into vector
   * registers, the cutoff criterion is checked, the particle properties are constructed, and the kernel is called. If
   * necessary, the globals are calculated.
   *
   * @tparam LowerPackedTrianglePB1.
   * @tparam LowerPackedTrianglePB2.
   * @tparam newton3.
   * @tparam newton3Kernel.
   * @param i current index of the i-loop.
   * @param j current index of the j-loop.
   * @param k current index of the k-loop.
   * @param const_nu nu value passed to the kernel (if a PPL is used, the values are constructed in the function body).
   * @param soaDists1 precomputed pairwise distances.
   * @param soaDists2 precomputed pairwise distances.
   * @param ownedStatePtr Pointer to the ownershipstate SoA of k.
   * @param typePtrI Pointer to the typeID of i.
   * @param typePtrJ Pointer to the typeID of j.
   * @param typePtrK Pointer to the typeID of k.
   * @param distXIJVec All lanes of a highway vector filled with the x-offset between i and j.
   * @param distYIJVec All lanes of a highway vector filled with the y-offset between i and j.
   * @param distZIJVec All lanes of a highway vector filled with the z-offset between i and j.
   * @param distSquaredIJVec All lanes of a highway vector filled with the sqared distance between i and j.
   * @param invR5IJVec All lanes of a highway vector filled with the inverse square root to 5 of the distance between i
   * and j.
   * @param fXAccI Force x accumulator for i as highway vector.
   * @param fYAccI Force y accumulator for i as highway vector.
   * @param fZAccI Force z accumulator for i as highway vector.
   * @param fXAccJ Force x accumulator for j as highway vector.
   * @param fYAccJ Force y accumulator for j as highway vector.
   * @param fZAccJ Force z accumulator for j as highway vector.
   * @param ownedStateI Ownership state of i.
   * @param ownedStateJ Ownership state of j.
   * @param fxPtr Pointer to the x-force SoA for k.
   * @param fyPtr Pointer to the y-force SoA for k.
   * @param fzPtr Pointer to the z-force SoA for k.
   * @param virialSumX x-Virial accumulator.
   * @param virialSumY y-Virial accumulator.
   * @param virialSumZ z-Virial accumulator.
   * @param potentialEnergySum Potential energy accumulator.
   * @param numKernelCallsN3Sum Kernel call accumulator newton3.
   * @param numGlobalCalcsN3Sum Globals calculation accumulator newton3.
   * @param numKernelCallsNoN3Sum Kernel call accumulator no newton3.
   * @param numGlobalCalcsNoN3Sum Globals calculation accumulator no newton3.
   * @param numLanesToProcess The number of lanes to be processed. This should only be less than the length of the
   * vector registers in the remainder case.
   */
  template <bool LowerPackedTrianglePB1, bool LowerPackedTrianglePB2, bool newton3, bool newton3Kernel>
  HWY_INLINE void handleKLoopBody(
      const size_t i, const size_t j, const size_t k, const SoAFloatPrecision const_nu,
      const PrecomputeBuffer &soaDists1, const PrecomputeBuffer &soaDists2,
      const autopas::OwnershipState *const __restrict ownedStatePtr, const size_t *const __restrict typePtrI,
      const size_t *const __restrict typePtrJ, const size_t *const __restrict typePtrK, const VectorDouble &distXIJVec,
      const VectorDouble &distYIJVec, const VectorDouble &distZIJVec, const VectorDouble &distSquaredIJVec,
      const VectorDouble &invR5IJVec, VectorDouble &fXAccI, VectorDouble &fYAccI, VectorDouble &fZAccI,
      VectorDouble &fXAccJ, VectorDouble &fYAccJ, VectorDouble &fZAccJ, const autopas::OwnershipState ownedStateI,
      const autopas::OwnershipState ownedStateJ, double *const __restrict fxPtr, double *const __restrict fyPtr,
      double *const __restrict fzPtr, VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ,
      VectorDouble &potentialEnergySum, size_t &numKernelCallsN3Sum, size_t &numGlobalCalcsN3Sum,
      size_t &numKernelCallsNoN3Sum, size_t &numGlobalCalcsNoN3Sum, size_t numLanesToProcess = _vecLengthDouble) {
    const auto remainderMask = mdLib::highway::FirstN(mdLib::tag_long, numLanesToProcess);

    // load distJKVecfrom precomputed data
    const auto [distXJKVec, distYJKVec, distZJKVec, distSquaredJKVec, invR5JKVec] =
        soaDists1.template loadRowVec<LowerPackedTrianglePB1>(j, k);
    // load distKIVec from precomputed data
    const auto [distXKIVecNeg, distYKIVecNeg, distZKIVecNeg, distSquaredKIVec, invR5KIVec] =
        soaDists2.template loadRowVec<LowerPackedTrianglePB2>(i, k);

    const auto distXKIVec = highway::Neg(distXKIVecNeg);
    const auto distYKIVec = highway::Neg(distYKIVecNeg);
    const auto distZKIVec = highway::Neg(distZKIVecNeg);

    // Since autopas::OwnershipState::dummy == 0, we can use MaskedLoad in the remainder-case, as it sets the remaining
    // lanes to 0
    const auto ownershipK =
        highway::MaskedLoad(remainderMask, tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr[k]));

    // calculate cutoff-masks between i<->k and j<->k
    const auto maskJK = highway::Le(distSquaredJKVec, highway::Set(tag_double, _cutoffSquared));
    const auto maskKI = highway::Le(distSquaredKIVec, highway::Set(tag_double, _cutoffSquared));
    // We check all lanes to see if they have dummy ownership state. In the kernel, only interactions with owned or halo
    // particles are calculated. We have to cast to a double mask because in the next step we perform a logical
    // operation with maskJK and maskKI (double masks). Google Highway only allows logical operations with masks of the
    // same type.
    const auto maskOwnershipK =
        highway::Ne(highway::ConvertTo(tag_double, ownershipK),
                    highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::dummy)));
    // Mask for triplets between i and j and multiple k
    const auto maskK = highway::And(highway::And(maskJK, maskKI), maskOwnershipK);

    if (HWY_UNLIKELY(highway::AllFalse(tag_double, maskK))) {
      return;
    }

    // load nus
    auto nu = highway::Set(tag_double, const_nu);
    if constexpr (useMixing) {
      const auto typeKIndices =
          highway::MaskedLoad(remainderMask, tag_long, reinterpret_cast<const int64_t *>(&typePtrK[k]));
      nu = _PPLibrary->get().getMixingNuHWY(typePtrI[i], typePtrJ[j], typeKIndices, numLanesToProcess);
    }

    // The call to SoAKernelHWY always requires forceJX, forceJY, forceJZ. In the case of newton3==false these buffers
    // are not used.
    VectorDouble forceIX, forceIY, forceIZ;
    VectorDouble forceJX, forceJY, forceJZ;
    VectorDouble factor, allDotProducts, allDistsSquared;
    SoAKernelHWY<newton3Kernel>(distXIJVec, distYIJVec, distZIJVec, distXJKVec, distYJKVec, distZJKVec, distXKIVec,
                                distYKIVec, distZKIVec, distSquaredIJVec, distSquaredJKVec, distSquaredKIVec, nu,
                                invR5IJVec, invR5JKVec, invR5KIVec, forceIX, forceIY, forceIZ, forceJX, forceJY,
                                forceJZ, factor, allDotProducts, allDistsSquared, maskK);
    if constexpr (not newton3) {
      // Compute only force I without Newton3
      fXAccI += forceIX;
      fYAccI += forceIY;
      fZAccI += forceIZ;

      if constexpr (newton3Kernel) {
        fXAccJ += forceJX;
        fYAccJ += forceJY;
        fZAccJ += forceJZ;
      }

      if constexpr (countFLOPs) {
        numKernelCallsNoN3Sum += numLanesToProcess;
      }

      if constexpr (calculateGlobals) {
        const auto remainderMaskDouble = highway::RebindMask(tag_double, remainderMask);

        const auto potentialEnergy3 = factor * (allDistsSquared - _threeDoubleVec * allDotProducts);

        // sum potential energy only on owned particles
        const auto ownedMaskI =
            highway::And(highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                                     highway::Set(tag_double, static_cast<double>(ownedStateI))),
                         remainderMaskDouble);

        potentialEnergySum += highway::IfThenElseZero(ownedMaskI, potentialEnergy3);
        // Virial for i
        const auto virialIX = forceIX * (distXKIVec - distXIJVec);
        const auto virialIY = forceIY * (distYKIVec - distYIJVec);
        const auto virialIZ = forceIZ * (distZKIVec - distZIJVec);

        const auto maskedVirialIX = highway::IfThenElseZero(ownedMaskI, virialIX);
        const auto maskedVirialIY = highway::IfThenElseZero(ownedMaskI, virialIY);
        const auto maskedVirialIZ = highway::IfThenElseZero(ownedMaskI, virialIZ);

        virialSumX += maskedVirialIX;
        virialSumY += maskedVirialIY;
        virialSumZ += maskedVirialIZ;

        if constexpr (newton3Kernel) {
          const auto ownedMaskJ =
              highway::And(highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                                       highway::Set(tag_double, static_cast<double>(ownedStateJ))),
                           remainderMaskDouble);

          potentialEnergySum += highway::IfThenElseZero(ownedMaskJ, potentialEnergy3);

          // Virial for j
          const auto virialJX = forceJX * (distXIJVec - distXJKVec);
          const auto virialJY = forceJY * (distYIJVec - distYJKVec);
          const auto virialJZ = forceJZ * (distZIJVec - distZJKVec);

          const auto maskedVirialJX = highway::IfThenElseZero(ownedMaskJ, virialJX);
          const auto maskedVirialJY = highway::IfThenElseZero(ownedMaskJ, virialJY);
          const auto maskedVirialJZ = highway::IfThenElseZero(ownedMaskJ, virialJZ);

          virialSumX += maskedVirialJX;
          virialSumY += maskedVirialJY;
          virialSumZ += maskedVirialJZ;
        }

        if constexpr (countFLOPs) {
          numGlobalCalcsNoN3Sum += numLanesToProcess;
        }
      }
    } else {
      // reduce forces
      fXAccI += forceIX;
      fYAccI += forceIY;
      fZAccI += forceIZ;

      fXAccJ += forceJX;
      fYAccJ += forceJY;
      fZAccJ += forceJZ;

      const auto forceKX = highway::Neg(forceIX + forceJX);
      const auto forceKY = highway::Neg(forceIY + forceJY);
      const auto forceKZ = highway::Neg(forceIZ + forceJZ);

      // Store force acting on particle k.
      const VectorDouble fxK = loadPacked(tag_double, &fxPtr[k]);
      const VectorDouble fyK = loadPacked(tag_double, &fyPtr[k]);
      const VectorDouble fzK = loadPacked(tag_double, &fzPtr[k]);

      const VectorDouble fxKNew = fxK + forceKX;
      const VectorDouble fyKNew = fyK + forceKY;
      const VectorDouble fzKNew = fzK + forceKZ;

      storePacked(tag_double, &fxPtr[k], fxKNew);
      storePacked(tag_double, &fyPtr[k], fyKNew);
      storePacked(tag_double, &fzPtr[k], fzKNew);

      if constexpr (countFLOPs) {
        numKernelCallsN3Sum += numLanesToProcess;
      }

      if constexpr (calculateGlobals) {
        const auto remainderMaskDouble = highway::RebindMask(tag_double, remainderMask);

        const auto potentialEnergy3 = factor * (allDistsSquared - _threeDoubleVec * allDotProducts);

        // sum potential energy only on owned particles
        const auto ownedMaskI =
            highway::And(highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                                     highway::Set(tag_double, static_cast<double>(ownedStateI))),
                         remainderMaskDouble);
        const auto ownedMaskJ =
            highway::And(highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                                     highway::Set(tag_double, static_cast<double>(ownedStateJ))),
                         remainderMaskDouble);
        const auto ownedMaskK =
            highway::Eq(highway::ConvertTo(tag_double, ownershipK),
                        highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)));

        // potential energy for i,j,k (masked)
        const auto potentialEnergyI = highway::IfThenElseZero(ownedMaskI, potentialEnergy3);
        const auto potentialEnergyJ = highway::IfThenElseZero(ownedMaskJ, potentialEnergy3);
        const auto potentialEnergyK = highway::IfThenElseZero(ownedMaskK, potentialEnergy3);

        potentialEnergySum += potentialEnergyI + potentialEnergyJ + potentialEnergyK;

        // Virial for i
        const auto virialIX = forceIX * (distXKIVec - distXIJVec);
        const auto virialIY = forceIY * (distYKIVec - distYIJVec);
        const auto virialIZ = forceIZ * (distZKIVec - distZIJVec);

        const auto maskedVirialIX = highway::IfThenElseZero(ownedMaskI, virialIX);
        const auto maskedVirialIY = highway::IfThenElseZero(ownedMaskI, virialIY);
        const auto maskedVirialIZ = highway::IfThenElseZero(ownedMaskI, virialIZ);

        // Virial for j
        const auto virialJX = forceJX * (distXIJVec - distXJKVec);
        const auto virialJY = forceJY * (distYIJVec - distYJKVec);
        const auto virialJZ = forceJZ * (distZIJVec - distZJKVec);

        const auto maskedVirialJX = highway::IfThenElseZero(ownedMaskJ, virialJX);
        const auto maskedVirialJY = highway::IfThenElseZero(ownedMaskJ, virialJY);
        const auto maskedVirialJZ = highway::IfThenElseZero(ownedMaskJ, virialJZ);

        // Virial for k
        const auto virialKX = forceKX * (distXJKVec - distXKIVec);
        const auto virialKY = forceKY * (distYJKVec - distYKIVec);
        const auto virialKZ = forceKZ * (distZJKVec - distZKIVec);

        const auto maskedVirialKX = highway::IfThenElseZero(ownedMaskK, virialKX);
        const auto maskedVirialKY = highway::IfThenElseZero(ownedMaskK, virialKY);
        const auto maskedVirialKZ = highway::IfThenElseZero(ownedMaskK, virialKZ);

        // Reduce the virial
        virialSumX += maskedVirialIX + maskedVirialJX + maskedVirialKX;
        virialSumY += maskedVirialIY + maskedVirialJY + maskedVirialKY;
        virialSumZ += maskedVirialIZ + maskedVirialJZ + maskedVirialKZ;

        if constexpr (countFLOPs) {
          numGlobalCalcsN3Sum += numLanesToProcess;
        }
      }
    }
  }

  /**
   * Inline helper to compute force components for particle I.
   * Returns tuple of (forceIX, forceIY, forceIZ, factor, allDotProducts, allDistsSquared)
   */
  template <bool newton3>
  HWY_INLINE auto SoAKernelHWY(const VectorDouble &distXIJ, const VectorDouble &distYIJ, const VectorDouble &distZIJ,
                               const VectorDouble &distXJK, const VectorDouble &distYJK, const VectorDouble &distZJK,
                               const VectorDouble &distXKI, const VectorDouble &distYKI, const VectorDouble &distZKI,
                               const VectorDouble &distSquaredIJ, const VectorDouble &distSquaredJK,
                               const VectorDouble &distSquaredKI, const VectorDouble &nu, const VectorDouble &invR5IJ,
                               const VectorDouble &invR5JK, const VectorDouble &invR5KI, VectorDouble &forceIX,
                               VectorDouble &forceIY, VectorDouble &forceIZ, VectorDouble &forceJX,
                               VectorDouble &forceJY, VectorDouble &forceJZ, VectorDouble &factor,
                               VectorDouble &allDotProducts, VectorDouble &allDistsSquared,
                               const MaskDouble &maskK) const {
    // computes a * b + c * d + e * f
    auto fmaHelper = [](auto a, auto b, auto c, auto d, auto e, auto f) {
      const auto x = highway::MulAdd(c, d, e * f);
      return highway::MulAdd(a, b, x);
    };

    // If globals are needed, still compute allDistsSquared (r_ij^2 * r_jk^2 * r_ki^2)
    if constexpr (calculateGlobals) {
      allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
    }

    // Use precomputed 1/r^5 for each pair: allInvDistsTo5 = (1/r_ij^5)*(1/r_jk^5)*(1/r_ki^5)
    const VectorDouble allInvDistsTo5 = invR5IJ * invR5JK * invR5KI;

    // factor is three * nu * product(inv r^5)
    factor = _threeDoubleVec * nu * allInvDistsTo5;

    // Dot products
    const VectorDouble IJDotKI = fmaHelper(distXIJ, distXKI, distYIJ, distYKI, distZIJ, distZKI);
    const VectorDouble IJDotJK = fmaHelper(distXIJ, distXJK, distYIJ, distYJK, distZIJ, distZJK);
    const VectorDouble JKDotKI = fmaHelper(distXJK, distXKI, distYJK, distYKI, distZJK, distZKI);
    allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    // Force I components
    const VectorDouble factorIDirectionJK = factor * IJDotKI * (IJDotJK - JKDotKI);

    const VectorDouble factorIDirectionIJ = factor * (highway::MulSub(IJDotJK, JKDotKI, distSquaredJK * distSquaredKI) +
                                                      _fiveDoubleVec * allDotProducts / distSquaredIJ);
    const VectorDouble factorIDirectionKI =
        factor * (highway::MulAdd(highway::Neg(IJDotJK), JKDotKI, distSquaredIJ * distSquaredJK) -
                  _fiveDoubleVec * allDotProducts / distSquaredKI);

    forceIX = fmaHelper(distXJK, factorIDirectionJK, distXIJ, factorIDirectionIJ, distXKI, factorIDirectionKI);
    forceIY = fmaHelper(distYJK, factorIDirectionJK, distYIJ, factorIDirectionIJ, distYKI, factorIDirectionKI);
    forceIZ = fmaHelper(distZJK, factorIDirectionJK, distZIJ, factorIDirectionIJ, distZKI, factorIDirectionKI);

    if constexpr (newton3) {
      // Force J components
      const VectorDouble factorJDirectionKI = factor * IJDotJK * (JKDotKI - IJDotKI);
      const VectorDouble factorJDirectionIJ =
          factor * (highway::MulAdd(highway::Neg(IJDotKI), JKDotKI, distSquaredJK * distSquaredKI) -
                    _fiveDoubleVec * allDotProducts / distSquaredIJ);

      const VectorDouble factorJDirectionJK =
          factor * (highway::MulSub(IJDotKI, JKDotKI, distSquaredIJ * distSquaredKI) +
                    _fiveDoubleVec * allDotProducts / distSquaredJK);

      forceJX = fmaHelper(distXKI, factorJDirectionKI, distXIJ, factorJDirectionIJ, distXJK, factorJDirectionJK);
      forceJY = fmaHelper(distYKI, factorJDirectionKI, distYIJ, factorJDirectionIJ, distYJK, factorJDirectionJK);
      forceJZ = fmaHelper(distZKI, factorJDirectionKI, distZIJ, factorJDirectionIJ, distZJK, factorJDirectionJK);
    }

    // apply masks
    forceIX = highway::IfThenElseZero(maskK, forceIX);
    forceIY = highway::IfThenElseZero(maskK, forceIY);
    forceIZ = highway::IfThenElseZero(maskK, forceIZ);

    if constexpr (newton3) {
      forceJX = highway::IfThenElseZero(maskK, forceJX);
      forceJY = highway::IfThenElseZero(maskK, forceJY);
      forceJZ = highway::IfThenElseZero(maskK, forceJZ);
    }

    if constexpr (calculateGlobals) {
      factor = highway::IfThenElseZero(maskK, factor);
      allDotProducts = highway::IfThenElseZero(maskK, allDotProducts);
      allDistsSquared = highway::IfThenElseZero(maskK, allDistsSquared);
    }
  }

  template <class Tag, typename T>
  HWY_INLINE static auto loadPacked(const Tag tag, const T *HWY_RESTRICT ptr) {
    return highway::Load(tag, ptr);
  }

  template <class Tag, typename T, typename Vec>
  HWY_INLINE static void storePacked(const Tag tag, T *HWY_RESTRICT ptr, const Vec &v) {
    highway::Store(v, tag, ptr);
  }

  template <class... SoAs>
  bool areAllSoAsAligned(const SoAs &...soas) {
    return (isSoAAligned(soas) and ...);
  }

  template <class SoA>
  bool isSoAAligned(const SoA &soa) {
    return isAttributeAligned<SoA, Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY,
                              Particle_T::AttributeNames::forceZ>(soa, alignmentSoAFloatHwyVector) and

           isAttributeAligned<SoA, Particle_T::AttributeNames::posX, Particle_T::AttributeNames::posY,
                              Particle_T::AttributeNames::posZ>(soa, alignmentSoAFloatHwyVector) and

           isAttributeAligned<SoA, Particle_T::AttributeNames::ownershipState>(soa, alignmentSoAFloatHwyVector) and

           isAttributeAligned<SoA, Particle_T::AttributeNames::typeId>(soa, alignmentSoAFloatHwyVector);
  }

  template <class SoA, auto... Attributes>
  HWY_INLINE bool isAttributeAligned(const SoA &soa, size_t alignment) {
    return (soa.template isAligned<Attributes>(alignment) and ...);
  }

  /**
   * Data structure that stores precomputed data. In the case of intra-SoA data, a lower packed triangle matrix can be
   * used to save memory space. In the case of inter-SoA data, a full matrix is occupied. Both data structures allow
   * vectorized loads and stores. However, in the case of a lower packed triangle matrix, no aligned vector instructions
   * are used.
   *
   * @tparam LowerTriangle If the matrix should use a lower packed triangle layout
   */
  struct PrecomputeBuffer {
   public:
    using Alloc = autopas::AlignedAllocator<SoAFloatPrecision>;
    using Vec = std::vector<SoAFloatPrecision, Alloc>;

    /**
     * Struct for storing vector registers with precalculated x, y, z displacements, the squared distance, and the
     * inverse r^5 between two particles.
     */
    struct DistRowVec {
      /**
       * x-displacement between two particles.
       */
      VectorDouble dx;
      /**
       * y-displacement between two particles.
       */
      VectorDouble dy;
      /**
       * z-displacement between two particles.
       */
      VectorDouble dz;
      /**
       * The squared distance between two particles.
       */
      VectorDouble r2;
      /**
       * The inverese distance to 5 between two particles (1.0/r^5).
       */
      VectorDouble invR5;
    };

    /**
     * Initializes the PrecomputeBuffer with the specified size. If the internal memory size is smaller than the
     * required size, it is increased.
     *
     * @tparam LowerTriangle
     * @param rows
     * @param cols
     */
    template <bool LowerTriangle>
    void init(const size_t rows, const size_t cols) {
      _nRows = rows;
      _nCols = cols;
      if constexpr (not LowerTriangle) {
        // round up to multiple of vector length
        _nColsPadded = ((_nCols + _vecLengthDouble - 1) / _vecLengthDouble) * _vecLengthDouble;
      }

      const size_t size = [&]() {
        if constexpr (LowerTriangle) {
          assert(rows == cols and "LowerTriangle requires a square matrix");

          const size_t logicalSize = (rows * rows - rows) / 2;

          // SIMD padding so that full vector loads are always safe
          return ((logicalSize + _vecLengthDouble - 1) / _vecLengthDouble) * _vecLengthDouble;
        } else {
          return rows * _nColsPadded;
        }
      }();

      resizeStorage(size);
    }

    // Lower-triangle index
    template <bool LowerTriangle>
    HWY_INLINE const size_t triIndex(const size_t i, const size_t j) const {
      static_assert(LowerTriangle, "triIndex only valid for LowerTriangle=true");
      assert(i > j);
      return i * (i - 1) / 2 + j;
    }

    // Full-matrix index
    template <bool LowerTriangle>
    HWY_INLINE const size_t fullIndex(const size_t i, const size_t j) const {
      static_assert(not LowerTriangle, "fullIndex only valid for LowerTriangle=false");
      return j + i * _nColsPadded;
    }

    template <bool LowerTriangle>
    HWY_INLINE const size_t index(const size_t i, const size_t j) const {
      if constexpr (LowerTriangle) {
        return triIndex<true>(i, j);
      } else {
        return fullIndex<false>(i, j);
      }
    }

    template <bool LowerTriangle>
    HWY_INLINE void set(const size_t i, const size_t j, const SoAFloatPrecision dx, const SoAFloatPrecision dy,
                        const SoAFloatPrecision dz, const SoAFloatPrecision dist2, const SoAFloatPrecision invr5) {
      const size_t idx = index(i, j);
      _dx[idx] = dx;
      _dy[idx] = dy;
      _dz[idx] = dz;
      _squared[idx] = dist2;
      _invR5[idx] = invr5;
    }

    template <bool LowerTriangle>
    HWY_INLINE auto get(const size_t i, const size_t j) const {
      const size_t idx = index<LowerTriangle>(i, j);
      return std::tuple{_dx[idx], _dy[idx], _dz[idx], _squared[idx], _invR5[idx]};
    }

    // Highway-Vector load
    template <bool LowerTriangle>
    HWY_INLINE DistRowVec loadRowVec(const size_t i, const size_t jStart) const {
      const size_t idx = [&]() {
        if constexpr (LowerTriangle) {
          return triIndex<true>(i, jStart);
        } else {
          return fullIndex<false>(i, jStart);
        }
      }();

      DistRowVec result;
      if constexpr (LowerTriangle) {
        result.dx = highway::LoadU(tag_double, &_dx[idx]);
        result.dy = highway::LoadU(tag_double, &_dy[idx]);
        result.dz = highway::LoadU(tag_double, &_dz[idx]);
        result.r2 = highway::LoadU(tag_double, &_squared[idx]);
        result.invR5 = highway::LoadU(tag_double, &_invR5[idx]);
      } else {
        result.dx = highway::Load(tag_double, &_dx[idx]);
        result.dy = highway::Load(tag_double, &_dy[idx]);
        result.dz = highway::Load(tag_double, &_dz[idx]);
        result.r2 = highway::Load(tag_double, &_squared[idx]);
        result.invR5 = highway::Load(tag_double, &_invR5[idx]);
      }
      return result;
    }

    template <bool LowerTriangle>
    __attribute__((no_sanitize("address")))
    HWY_INLINE void fillHighway(const SoAFloatPrecision *xptr1, const SoAFloatPrecision *yptr1,
                                const SoAFloatPrecision *zptr1, const SoAFloatPrecision *xptr2,
                                const SoAFloatPrecision *yptr2, const SoAFloatPrecision *zptr2,
                                const SoAFloatPrecision cutoffSquared) {
      for (size_t i = 0; i < _nRows; ++i) {
        const auto x1v = highway::Set(tag_double, xptr1[i]);
        const auto y1v = highway::Set(tag_double, yptr1[i]);
        const auto z1v = highway::Set(tag_double, zptr1[i]);

        const auto cutoffSquaredV = highway::Set(tag_double, cutoffSquared);

        size_t j = 0;
        size_t rowLength = 0;

        if constexpr (LowerTriangle) {
          rowLength = i;
          if (i == 0) continue;
        } else {
          rowLength = _nColsPadded;
        }

        auto checkLoopCondition = [&]() {
          if constexpr (LowerTriangle) {
            return j + _vecLengthDouble <= rowLength;
          } else {
            return j < rowLength;
          }
        };

        for (; checkLoopCondition(); j += _vecLengthDouble) {
          const auto x2v = loadPacked(tag_double, &xptr2[j]);
          const auto y2v = loadPacked(tag_double, &yptr2[j]);
          const auto z2v = loadPacked(tag_double, &zptr2[j]);

          const auto dx = x2v - x1v;
          const auto dy = y2v - y1v;
          const auto dz = z2v - z1v;

          const auto dist2 = highway::MulAdd(dx, dx, highway::MulAdd(dy, dy, dz * dz));
          const auto cutoffMask = highway::Le(dist2, cutoffSquaredV);
          const auto sqrtR2 = highway::Sqrt(dist2);
          const auto r5 = highway::Mul(highway::Mul(dist2, dist2), sqrtR2);
          const auto invr5 = highway::MaskedDiv(cutoffMask, highway::Set(tag_double, 1.0), r5);

          const size_t idx = index<LowerTriangle>(i, j);

          if constexpr (LowerTriangle) {
            highway::StoreU(dx, tag_double, &_dx[idx]);
            highway::StoreU(dy, tag_double, &_dy[idx]);
            highway::StoreU(dz, tag_double, &_dz[idx]);
            highway::StoreU(dist2, tag_double, &_squared[idx]);
            highway::StoreU(invr5, tag_double, &_invR5[idx]);
          } else {
            highway::Store(dx, tag_double, &_dx[idx]);
            highway::Store(dy, tag_double, &_dy[idx]);
            highway::Store(dz, tag_double, &_dz[idx]);
            highway::Store(dist2, tag_double, &_squared[idx]);
            highway::Store(invr5, tag_double, &_invR5[idx]);
          }
        }

        // Remainder (only used in the LowerTriangle case)
        if (j < rowLength) {
          const size_t width = rowLength - j;

          const auto x2v = highway::LoadN(tag_double, &xptr2[j], width);
          const auto y2v = highway::LoadN(tag_double, &yptr2[j], width);
          const auto z2v = highway::LoadN(tag_double, &zptr2[j], width);

          const auto dx = x2v - x1v;
          const auto dy = y2v - y1v;
          const auto dz = z2v - z1v;

          const auto dist2 = highway::MulAdd(dx, dx, highway::MulAdd(dy, dy, dz * dz));
          const auto cutoffMask = highway::Le(dist2, highway::Set(tag_double, cutoffSquared));
          const auto sqrtR2 = highway::Sqrt(dist2);
          const auto r5 = highway::Mul(highway::Mul(dist2, dist2), sqrtR2);
          const auto invr5 = highway::MaskedDiv(cutoffMask, highway::Set(tag_double, 1.0), r5);

          const size_t idx = index<LowerTriangle>(i, j);

          highway::StoreN(dx, tag_double, &_dx[idx], width);
          highway::StoreN(dy, tag_double, &_dy[idx], width);
          highway::StoreN(dz, tag_double, &_dz[idx], width);
          highway::StoreN(dist2, tag_double, &_squared[idx], width);
          highway::StoreN(invr5, tag_double, &_invR5[idx], width);
        }
      }
    }

   private:
    void resizeStorage(size_t size) {
      if (size > _currentStorageSize) {
        _dx.resize(size);
        _dy.resize(size);
        _dz.resize(size);
        _squared.resize(size);
        _invR5.resize(size);
        _currentStorageSize = size;
      }
    }

    Vec _dx;
    Vec _dy;
    Vec _dz;
    Vec _squared;
    Vec _invR5;
    size_t _currentStorageSize{0};

    size_t _nRows{};
    size_t _nCols{};
    size_t _nColsPadded{};
  };

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

  // optional to hold a reference to the ParticlePropertiesLibrary. If a ParticlePropertiesLibrary is not used the
  // optional is empty.
  std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>> _PPLibrary;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals;
  std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // Helper vectors that are often used
  const VectorDouble _threeDoubleVec = highway::Set(tag_double, 3.0);
  const VectorDouble _fiveDoubleVec = highway::Set(tag_double, 5.0);
  const VectorDouble _zeroDoubleVec = highway::Zero(tag_double);

  // Alignment for SoAFloatPrecision vector register
  const std::size_t alignmentSoAFloatHwyVector = _vecLengthDouble * sizeof(SoAFloatPrecision);

  // Precomute buffers for distances and invR5 between particles. They are class members to avoid excessive memory
  // allocation and deallocation.
  std::vector<PrecomputeBuffer> precomputeBuffers1;
  std::vector<PrecomputeBuffer> precomputeBuffers2;
};
}  // namespace mdLib
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

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit AxilrodTellerMutoFunctorHWY(double cutoff, void * /*dummy*/)
      : autopas::TriwiseFunctor<
            Particle_T, AxilrodTellerMutoFunctorHWY<Particle_T, useMixing, useNewton3, calculateGlobals, countFLOPs>>(
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
      AutoPasLog(DEBUG, "Using AxilrodTellerMutoFunctorHWY with countFLOPs is not tested for SoA datalayout.");
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
  explicit AxilrodTellerMutoFunctorHWY(double cutoff) : AxilrodTellerMutoFunctorHWY(cutoff, nullptr) {
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
  explicit AxilrodTellerMutoFunctorHWY(double cutoff,
                                       ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : AxilrodTellerMutoFunctorHWY(cutoff, nullptr) {
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

    const bool alignedForce =
        isAligned3D<decltype(soa), Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY,
                    Particle_T::AttributeNames::forceZ>(soa, alignmentSoAFloatHwyVector);

    if (newton3) {
      if (alignedForce) {
        SoAFunctorSingleImpl<true, true>(soa);
      } else {
        SoAFunctorSingleImpl<true, false>(soa);
      }
    } else {
      if (alignedForce) {
        SoAFunctorSingleImpl<false, true>(soa);
      } else {
        SoAFunctorSingleImpl<false, false>(soa);
      }
    }
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final {
    if (soa1.size() == 0 || soa2.size() == 0) return;

    const bool alignedForce =
        isAligned3D<decltype(soa2), Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY,
                    Particle_T::AttributeNames::forceZ>(soa2, alignmentSoAFloatHwyVector);

    if (newton3) {
      if (alignedForce) {
        SoAFunctorPairImpl<true, true>(soa1, soa2);
      } else {
        SoAFunctorPairImpl<true, false>(soa1, soa2);
      }
    } else {
      if (alignedForce) {
        SoAFunctorPairImpl<false, true>(soa1, soa2);
      } else {
        SoAFunctorPairImpl<false, false>(soa1, soa2);
      }
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
  struct DistanceMatrix;

  template <bool newton3, bool alignedSoAView>
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
        highway::Zero(tag_double);  // Note: This is not the potential energy but some fixed multiple of it.
    VectorDouble virialSumX = highway::Zero(tag_double);
    VectorDouble virialSumY = highway::Zero(tag_double);
    VectorDouble virialSumZ = highway::Zero(tag_double);

    size_t numTripletsCountingSum = 0;
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numGlobalCalcsSum = 0;

    const SoAFloatPrecision const_nu = _nu;
    const size_t soaSize = soa.size();
    const size_t soaSizeSquared = soaSize * soaSize;

    if constexpr (countFLOPs) {
      numTripletsCountingSum = soaSize * (soaSize - 1) * (soaSize - 2) / 6;
    }

    // Precompute distances between particles in the soa. We use a row-major packed lower‐triangle data structure to
    // save memory. If a distance i->j is required as j->i it can be negated.
    DistanceMatrix intraSoA(soaSize, soaSize, true);
    // soa1 <-> soa1
    for (size_t i = 1; i < soaSize; ++i) {
      for (size_t j = 0; j < i; ++j) {
        const auto dx = xPtr[j] - xPtr[i];
        const auto dy = yPtr[j] - yPtr[i];
        const auto dz = zPtr[j] - zPtr[i];
        const auto dist2 = dx * dx + dy * dy + dz * dz;
        const double invr5 = (dist2 <= cutoffSquared) ? 1.0 / (dist2 * dist2 * std::sqrt(dist2)) : 0.0;
        intraSoA.set(i, j, dx, dy, dz, dist2, invr5);
      }
    }

    const size_t packedSize = (soaSizeSquared - soaSize) / 2;
    if constexpr (countFLOPs) {
      numDistanceCalculationSum += packedSize;
    }

    // We iterate in reverse order over i and j to enable aligned loads in the k-loop.
    for (size_t i = soaSize - 1; i > 1; --i) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }
      VectorDouble fXAccI = highway::Set(tag_double, 0.);
      VectorDouble fYAccI = highway::Set(tag_double, 0.);
      VectorDouble fZAccI = highway::Set(tag_double, 0.);

      for (size_t j = i - 1; j > 0; --j) {
        const auto ownedStateJ = ownedStatePtr[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        auto [distXIJ, distYIJ, distZIJ, distSquaredIJ, invR5IJ] = intraSoA.get(i, j);

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        VectorDouble fXAccJ = highway::Set(tag_double, 0.);
        VectorDouble fYAccJ = highway::Set(tag_double, 0.);
        VectorDouble fZAccJ = highway::Set(tag_double, 0.);

        // prepare vector registers

        const auto distXIJVec = highway::Set(tag_double, distXIJ);
        const auto distYIJVec = highway::Set(tag_double, distYIJ);
        const auto distZIJVec = highway::Set(tag_double, distZIJ);
        const auto distSquaredIJVec = highway::Set(tag_double, distSquaredIJ);
        const auto invR5IJVec = highway::Set(tag_double, invR5IJ);

        const size_t blockEnd = j & ~(_vecLengthDouble - 1);

        unsigned int k = 0;
        for (; k < blockEnd; k += _vecLengthDouble) {
          // const bool remainder = k >= blockEnd;

          handleKLoopBody<newton3, /*remainder*/ false, alignedSoAView>(
              i, j, k, const_nu, intraSoA, intraSoA, ownedStatePtr, typePtr, distXIJVec, distYIJVec, distZIJVec,
              distSquaredIJVec, invR5IJVec, fXAccI, fYAccI, fZAccI, fXAccJ, fYAccJ, fZAccJ, ownedStateI, ownedStateJ,
              fxPtr, fyPtr, fzPtr, virialSumX, virialSumY, virialSumZ, potentialEnergySum, numKernelCallsN3Sum,
              numGlobalCalcsSum);
        }

        const auto restK = j - k;

        handleKLoopBody<newton3, /*remainder*/ true, alignedSoAView>(
            i, j, k, const_nu, intraSoA, intraSoA, ownedStatePtr, typePtr, distXIJVec, distYIJVec, distZIJVec,
            distSquaredIJVec, invR5IJVec, fXAccI, fYAccI, fZAccI, fXAccJ, fYAccJ, fZAccJ, ownedStateI, ownedStateJ,
            fxPtr, fyPtr, fzPtr, virialSumX, virialSumY, virialSumZ, potentialEnergySum, numKernelCallsN3Sum,
            numGlobalCalcsSum, restK);

        fxPtr[j] += highway::ReduceSum(tag_double, fXAccJ);
        fyPtr[j] += highway::ReduceSum(tag_double, fYAccJ);
        fzPtr[j] += highway::ReduceSum(tag_double, fZAccJ);
      }
      fxPtr[i] += highway::ReduceSum(tag_double, fXAccI);
      fyPtr[i] += highway::ReduceSum(tag_double, fYAccI);
      fzPtr[i] += highway::ReduceSum(tag_double, fZAccI);
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numTripletsCount += numTripletsCountingSum;
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsSum;  // Always N3 in Single SoAFunctor
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += highway::ReduceSum(tag_double, potentialEnergySum);
      _aosThreadDataGlobals[threadnum].virialSum[0] += highway::ReduceSum(tag_double, virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += highway::ReduceSum(tag_double, virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += highway::ReduceSum(tag_double, virialSumZ);
    }
  }

  template <bool newton3, bool alignedSoAView>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    const auto threadnum = autopas::autopas_get_thread_num();

    const auto alignment = highway::Lanes(tag_double) * sizeof(SoAFloatPrecision);
    bool alignedPosSoA2 = soa2.template isAligned<Particle_T::AttributeNames::posX>(alignment) and
                          soa2.template isAligned<Particle_T::AttributeNames::posY>(alignment) and
                          soa2.template isAligned<Particle_T::AttributeNames::posZ>(alignment);
    bool alignedForceSoA2 = soa2.template isAligned<Particle_T::AttributeNames::forceX>(alignment) and
                            soa2.template isAligned<Particle_T::AttributeNames::forceY>(alignment) and
                            soa2.template isAligned<Particle_T::AttributeNames::forceZ>(alignment);

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

    VectorDouble potentialEnergySum =
        highway::Zero(tag_double);  // Note: This is not the potential energy but some fixed multiple of it.
    VectorDouble virialSumX = highway::Zero(tag_double);
    VectorDouble virialSumY = highway::Zero(tag_double);
    VectorDouble virialSumZ = highway::Zero(tag_double);

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

    DistanceMatrix intraSoA1Dists(soa1Size, soa1Size, true);
    DistanceMatrix intraSoA2Dists(soa2Size, soa2Size, true);
    DistanceMatrix interSoADists(soa1Size, soa2Size, false);

    // soa1 <-> soa1
    for (size_t i = 1; i < soa1Size; ++i) {
      for (size_t j = 0; j < i; ++j) {
        const auto dx = xptr1[j] - xptr1[i];
        const auto dy = yptr1[j] - yptr1[i];
        const auto dz = zptr1[j] - zptr1[i];
        const auto dist2 = dx * dx + dy * dy + dz * dz;
        const double invr5 = (dist2 <= cutoffSquared) ? 1.0 / (dist2 * dist2 * std::sqrt(dist2)) : 0.0;
        intraSoA1Dists.set(i, j, dx, dy, dz, dist2, invr5);
      }
    }

    // soa2 <-> soa2
    for (size_t i = 1; i < soa2Size; ++i) {
      for (size_t j = 0; j < i; ++j) {
        const auto dx = xptr2[j] - xptr2[i];
        const auto dy = yptr2[j] - yptr2[i];
        const auto dz = zptr2[j] - zptr2[i];
        const auto dist2 = dx * dx + dy * dy + dz * dz;
        const double invr5 = (dist2 <= cutoffSquared) ? 1.0 / (dist2 * dist2 * std::sqrt(dist2)) : 0.0;
        intraSoA2Dists.set(i, j, dx, dy, dz, dist2, invr5);
      }
    }

    // soa1 <-> soa2
    for (size_t i = 0; i < soa1Size; ++i) {
      for (size_t j = 0; j < soa2Size; ++j) {
        const auto dx = xptr2[j] - xptr1[i];
        const auto dy = yptr2[j] - yptr1[i];
        const auto dz = zptr2[j] - zptr1[i];
        const auto dist2 = dx * dx + dy * dy + dz * dz;
        const double invr5 = (dist2 <= cutoffSquared) ? 1.0 / (dist2 * dist2 * std::sqrt(dist2)) : 0.0;
        interSoADists.set(i, j, dx, dy, dz, dist2, invr5);
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
      VectorDouble fXAccI = highway::Set(tag_double, 0.);
      VectorDouble fYAccI = highway::Set(tag_double, 0.);
      VectorDouble fZAccI = highway::Set(tag_double, 0.);

      // CASE: Particle i is in soa1, j and k are both in soa2
      // We iterate in reverse order over j to enable aligned stores in the k-loop.
      for (unsigned int j = soa2Size - 1; j > 0; --j) {
        const auto ownedStateJ = ownedStatePtr2[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        auto [distXIJ, distYIJ, distZIJ, distSquaredIJ, invR5IJ] = interSoADists.get(i, j);

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        VectorDouble fXAccJ = highway::Set(tag_double, 0.);
        VectorDouble fYAccJ = highway::Set(tag_double, 0.);
        VectorDouble fZAccJ = highway::Set(tag_double, 0.);

        // prepare vector registers
        const auto distXIJVec = highway::Set(tag_double, distXIJ);
        const auto distYIJVec = highway::Set(tag_double, distYIJ);
        const auto distZIJVec = highway::Set(tag_double, distZIJ);
        const auto distSquaredIJVec = highway::Set(tag_double, distSquaredIJ);
        const auto invR5IJVec = highway::Set(tag_double, invR5IJ);

        const size_t blockEnd = j & ~(_vecLengthDouble - 1);

        for (unsigned int k = 0; k < j; k += _vecLengthDouble) {
          const bool remainder = k >= blockEnd;
          const auto restK = j - k;

          const auto [distXJKVec, distYJKVec, distZJKVec, distSquaredJKVec, invR5JKVec] =
              intraSoA2Dists.loadRowVec(j, k, restK, remainder);

          const auto [distXKIVecNeg, distYKIVecNeg, distZKIVecNeg, distSquaredKIVec, invR5KIVec] =
              interSoADists.loadRowVec(i, k, restK, remainder);

          const auto distXKIVec = highway::Neg(distXKIVecNeg);
          const auto distYKIVec = highway::Neg(distYKIVecNeg);
          const auto distZKIVec = highway::Neg(distZKIVecNeg);

          const auto mask = highway::FirstN(tag_long, remainder ? restK : _vecLengthDouble);
          const auto ownershipK =
              highway::MaskedLoadOr(highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::dummy)), mask,
                                    tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr2[k]));

          // calculate masks for cutoff between i<->k and j<->k
          const auto maskJK = highway::Le(distSquaredJKVec, highway::Set(tag_double, _cutoffSquared));
          const auto maskKI = highway::Le(distSquaredKIVec, highway::Set(tag_double, _cutoffSquared));
          // ownership mask
          const auto maskOwnershipK =
              highway::Ne(highway::ConvertTo(tag_double, ownershipK),
                          highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::dummy)));
          // mask for triplets between i and j and multiple k
          const auto maskK = highway::And(highway::And(maskJK, maskKI), maskOwnershipK);

          if (highway::AllFalse(tag_double, maskK)) {
            continue;
          }

          // load nus
          auto nu = highway::Set(tag_double, const_nu);
          if constexpr (useMixing) {
            HWY_ALIGN double nus[_vecLengthDouble] = {0.};

            for (size_t n = 0; n < (remainder ? restK : _vecLengthDouble); ++n) {
              nus[n] = _PPLibrary->getMixingNu(typeptr1[i], typeptr2[j], typeptr2[k + n]);
            }
            nu = highway::Load(tag_double, nus);
          }

          if constexpr (not newton3) {
            // Compute only force I without Newton3
            VectorDouble forceIX, forceIY, forceIZ;
            VectorDouble forceJX, forceJY, forceJZ;
            VectorDouble factor, allDotProducts, allDistsSquared;
            SoAKernel_HWY<false>(distXIJVec, distYIJVec, distZIJVec, distXJKVec, distYJKVec, distZJKVec, distXKIVec,
                                 distYKIVec, distZKIVec, distSquaredIJVec, distSquaredJKVec, distSquaredKIVec, nu,
                                 invR5IJVec, invR5JKVec, invR5KIVec, forceIX, forceIY, forceIZ, forceJX, forceJY,
                                 forceJZ, factor, allDotProducts, allDistsSquared, maskK);

            fXAccI += forceIX;
            fYAccI += forceIY;
            fZAccI += forceIZ;

            if constexpr (countFLOPs) {
              ++numKernelCallsNoN3Sum;
            }
            if constexpr (calculateGlobals) {
              const auto potentialEnergy3 = factor * (allDistsSquared - _threeDoubleVec * allDotProducts);

              // sum potential energy only on owned particles
              const auto ownedMaskI =
                  highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                              highway::Set(tag_double, static_cast<double>(ownedStateI)));

              potentialEnergySum += highway::IfThenElse(ownedMaskI, potentialEnergy3, highway::Zero(tag_double));
              // Virial for i
              const auto virialIX = forceIX * (distXKIVec - distXIJVec);
              const auto virialIY = forceIY * (distYKIVec - distYIJVec);
              const auto virialIZ = forceIZ * (distZKIVec - distZIJVec);

              const auto maskedVirialIX = highway::IfThenElse(ownedMaskI, virialIX, highway::Zero(tag_double));
              const auto maskedVirialIY = highway::IfThenElse(ownedMaskI, virialIY, highway::Zero(tag_double));
              const auto maskedVirialIZ = highway::IfThenElse(ownedMaskI, virialIZ, highway::Zero(tag_double));
              virialSumX += maskedVirialIX;
              virialSumY += maskedVirialIY;
              virialSumZ += maskedVirialIZ;

              if constexpr (countFLOPs) {
                numGlobalCalcsNoN3Sum += remainder ? restK : _vecLengthDouble;
              }
            }
          } else {
            // Compute all forces with Newton3
            VectorDouble forceIX, forceIY, forceIZ;
            VectorDouble forceJX, forceJY, forceJZ;
            VectorDouble factor, allDotProducts, allDistsSquared;
            SoAKernel_HWY<true>(distXIJVec, distYIJVec, distZIJVec, distXJKVec, distYJKVec, distZJKVec, distXKIVec,
                                distYKIVec, distZKIVec, distSquaredIJVec, distSquaredJKVec, distSquaredKIVec, nu,
                                invR5IJVec, invR5JKVec, invR5KIVec, forceIX, forceIY, forceIZ, forceJX, forceJY,
                                forceJZ, factor, allDotProducts, allDistsSquared, maskK);

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
            const VectorDouble fx2 =
                remainder ? highway::LoadN(tag_double, &fxptr2[k], restK) : highway::LoadU(tag_double, &fxptr2[k]);
            const VectorDouble fy2 =
                remainder ? highway::LoadN(tag_double, &fyptr2[k], restK) : highway::LoadU(tag_double, &fyptr2[k]);
            const VectorDouble fz2 =
                remainder ? highway::LoadN(tag_double, &fzptr2[k], restK) : highway::LoadU(tag_double, &fzptr2[k]);

            const VectorDouble fx2New = fx2 + forceKX;
            const VectorDouble fy2New = fy2 + forceKY;
            const VectorDouble fz2New = fz2 + forceKZ;

            storePacked<alignedSoAView>(tag_double, &fxptr2[k], fx2New, remainder ? restK : 0);
            storePacked<alignedSoAView>(tag_double, &fyptr2[k], fy2New, remainder ? restK : 0);
            storePacked<alignedSoAView>(tag_double, &fzptr2[k], fz2New, remainder ? restK : 0);

            if constexpr (countFLOPs) {
              ++numKernelCallsN3Sum;
            }

            if constexpr (calculateGlobals) {
              const auto potentialEnergy3 = factor * (allDistsSquared - _threeDoubleVec * allDotProducts);

              // sum potential energy only on owned particles
              const auto ownedMaskI =
                  highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                              highway::Set(tag_double, static_cast<double>(ownedStateI)));
              const auto ownedMaskJ =
                  highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                              highway::Set(tag_double, static_cast<double>(ownedStateJ)));
              const auto ownedMaskK =
                  highway::Eq(highway::ConvertTo(tag_double, ownershipK),
                              highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)));

              // potential energy for i,j,k (masked)
              const auto potentialEnergyI =
                  highway::IfThenElse(ownedMaskI, potentialEnergy3, highway::Zero(tag_double));
              const auto potentialEnergyJ =
                  highway::IfThenElse(ownedMaskJ, potentialEnergy3, highway::Zero(tag_double));
              const auto potentialEnergyK =
                  highway::IfThenElse(ownedMaskK, potentialEnergy3, highway::Zero(tag_double));

              potentialEnergySum += potentialEnergyI + potentialEnergyJ + potentialEnergyK;

              // Virial for i
              const auto virialIX = forceIX * (distXKIVec - distXIJVec);
              const auto virialIY = forceIY * (distYKIVec - distYIJVec);
              const auto virialIZ = forceIZ * (distZKIVec - distZIJVec);

              const auto maskedVirialIX = highway::IfThenElse(ownedMaskI, virialIX, highway::Zero(tag_double));
              const auto maskedVirialIY = highway::IfThenElse(ownedMaskI, virialIY, highway::Zero(tag_double));
              const auto maskedVirialIZ = highway::IfThenElse(ownedMaskI, virialIZ, highway::Zero(tag_double));

              // Virial for j
              const auto virialJX = forceJX * (distXIJVec - distXJKVec);
              const auto virialJY = forceJY * (distYIJVec - distYJKVec);
              const auto virialJZ = forceJZ * (distZIJVec - distZJKVec);

              const auto maskedVirialJX = highway::IfThenElse(ownedMaskJ, virialJX, highway::Zero(tag_double));
              const auto maskedVirialJY = highway::IfThenElse(ownedMaskJ, virialJY, highway::Zero(tag_double));
              const auto maskedVirialJZ = highway::IfThenElse(ownedMaskJ, virialJZ, highway::Zero(tag_double));

              // Virial for k
              const auto virialKX = forceKX * (distXJKVec - distXKIVec);
              const auto virialKY = forceKY * (distYJKVec - distYKIVec);
              const auto virialKZ = forceKZ * (distZJKVec - distZKIVec);

              const auto maskedVirialKX = highway::IfThenElse(ownedMaskK, virialKX, highway::Zero(tag_double));
              const auto maskedVirialKY = highway::IfThenElse(ownedMaskK, virialKY, highway::Zero(tag_double));
              const auto maskedVirialKZ = highway::IfThenElse(ownedMaskK, virialKZ, highway::Zero(tag_double));

              // Reduce the virial
              virialSumX += maskedVirialIX + maskedVirialJX + maskedVirialKX;
              virialSumY += maskedVirialIY + maskedVirialJY + maskedVirialKY;
              virialSumZ += maskedVirialIZ + maskedVirialJZ + maskedVirialKZ;

              if constexpr (countFLOPs) {
                numGlobalCalcsN3Sum += remainder ? restK : _vecLengthDouble;
              }
            }
          }
        }
        fxptr2[j] += highway::ReduceSum(tag_double, fXAccJ);
        fyptr2[j] += highway::ReduceSum(tag_double, fYAccJ);
        fzptr2[j] += highway::ReduceSum(tag_double, fZAccJ);
      }

      // CASE: Particle i and j are both in soa1, k is in soa2
      for (unsigned int j = i + 1; j < soa1.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr1[j];
        if (ownedStateJ == autopas::OwnershipState::dummy) {
          continue;
        }

        // we have to mirror i and j because we iterate with i=0... and j=i+1...
        auto [distXIJ, distYIJ, distZIJ, distSquaredIJ, invR5IJ] = intraSoA1Dists.get(j, i);
        distXIJ = -distXIJ;
        distYIJ = -distYIJ;
        distZIJ = -distZIJ;

        if (distSquaredIJ > cutoffSquared) {
          continue;
        }

        VectorDouble fXAccJ = highway::Set(tag_double, 0.);
        VectorDouble fYAccJ = highway::Set(tag_double, 0.);
        VectorDouble fZAccJ = highway::Set(tag_double, 0.);

        // prepare vector registers
        const auto distXIJVec = highway::Set(tag_double, distXIJ);
        const auto distYIJVec = highway::Set(tag_double, distYIJ);
        const auto distZIJVec = highway::Set(tag_double, distZIJ);
        const auto distSquaredIJVec = highway::Set(tag_double, distSquaredIJ);
        const auto invR5IJVec = highway::Set(tag_double, invR5IJ);

        for (unsigned int k = 0; k < soa2.size(); k += _vecLengthDouble) {
          const bool remainder = (k + _vecLengthDouble > soa2.size());
          const size_t restK = soa2.size() - k;

          const auto [distXJKVec, distYJKVec, distZJKVec, distSquaredJKVec, invR5JKVec] =
              interSoADists.loadRowVec(j, k, restK, remainder);

          const auto [distXKIVecNeg, distYKIVecNeg, distZKIVecNeg, distSquaredKIVec, invR5KIVec] =
              interSoADists.loadRowVec(i, k, restK, remainder);

          const auto distXKIVec = highway::Neg(distXKIVecNeg);
          const auto distYKIVec = highway::Neg(distYKIVecNeg);
          const auto distZKIVec = highway::Neg(distZKIVecNeg);

          const auto mask = highway::FirstN(tag_long, remainder ? restK : _vecLengthDouble);
          const auto ownershipK =
              highway::MaskedLoadOr(highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::dummy)), mask,
                                    tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr2[k]));

          // calculate masks for cutoff between i<->k and j<->k
          const auto maskJK = highway::Le(distSquaredJKVec, highway::Set(tag_double, _cutoffSquared));
          const auto maskKI = highway::Le(distSquaredKIVec, highway::Set(tag_double, _cutoffSquared));
          // ownership mask
          const auto maskOwnershipK =
              highway::Ne(highway::ConvertTo(tag_double, ownershipK),
                          highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::dummy)));
          // mask for triplets between i and j and multiple k
          const auto maskK = highway::And(highway::And(maskJK, maskKI), maskOwnershipK);

          if (highway::AllFalse(tag_double, maskK)) {
            continue;
          }

          // load nus
          auto nu = highway::Set(tag_double, const_nu);
          if constexpr (useMixing) {
            HWY_ALIGN double nus[_vecLengthDouble] = {0.};

            for (size_t n = 0; n < (remainder ? restK : _vecLengthDouble); ++n) {
              nus[n] = _PPLibrary->getMixingNu(typeptr1[i], typeptr1[j], typeptr2[k + n]);
            }
            nu = highway::Load(tag_double, nus);
          }

          VectorDouble forceIX, forceIY, forceIZ;
          VectorDouble forceJX, forceJY, forceJZ;
          VectorDouble factor, allDotProducts, allDistsSquared;
          SoAKernel_HWY<true>(distXIJVec, distYIJVec, distZIJVec, distXJKVec, distYJKVec, distZJKVec, distXKIVec,
                              distYKIVec, distZKIVec, distSquaredIJVec, distSquaredJKVec, distSquaredKIVec, nu,
                              invR5IJVec, invR5JKVec, invR5KIVec, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ,
                              factor, allDotProducts, allDistsSquared, maskK);

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
            const auto forceKX = highway::Neg(forceIX + forceJX);
            const auto forceKY = highway::Neg(forceIY + forceJY);
            const auto forceKZ = highway::Neg(forceIZ + forceJZ);

            // Store force acting on particle k.
            const VectorDouble fx2 =
                remainder ? highway::LoadN(tag_double, &fxptr2[k], restK) : highway::LoadU(tag_double, &fxptr2[k]);
            const VectorDouble fy2 =
                remainder ? highway::LoadN(tag_double, &fyptr2[k], restK) : highway::LoadU(tag_double, &fyptr2[k]);
            const VectorDouble fz2 =
                remainder ? highway::LoadN(tag_double, &fzptr2[k], restK) : highway::LoadU(tag_double, &fzptr2[k]);

            const VectorDouble fx2New = fx2 + forceKX;
            const VectorDouble fy2New = fy2 + forceKY;
            const VectorDouble fz2New = fz2 + forceKZ;

            storePacked<alignedSoAView>(tag_double, &fxptr2[k], fx2New, remainder ? restK : 0);
            storePacked<alignedSoAView>(tag_double, &fyptr2[k], fy2New, remainder ? restK : 0);
            storePacked<alignedSoAView>(tag_double, &fzptr2[k], fz2New, remainder ? restK : 0);

            if constexpr (calculateGlobals) {
              const auto potentialEnergy3 = factor * (allDistsSquared - _threeDoubleVec * allDotProducts);

              // sum potential energy only on owned particles
              const auto ownedMaskK =
                  highway::Eq(highway::ConvertTo(tag_double, ownershipK),
                              highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)));

              potentialEnergySum += highway::IfThenElse(ownedMaskK, potentialEnergy3, highway::Zero(tag_double));
              // Virial for k
              const auto virialKX = forceKX * (distXJKVec - distXKIVec);
              const auto virialKY = forceKY * (distYJKVec - distYKIVec);
              const auto virialKZ = forceKZ * (distZJKVec - distZKIVec);

              const auto maskedVirialKX = highway::IfThenElse(ownedMaskK, virialKX, highway::Zero(tag_double));
              const auto maskedVirialKY = highway::IfThenElse(ownedMaskK, virialKY, highway::Zero(tag_double));
              const auto maskedVirialKZ = highway::IfThenElse(ownedMaskK, virialKZ, highway::Zero(tag_double));

              virialSumX += maskedVirialKX;
              virialSumY += maskedVirialKY;
              virialSumZ += maskedVirialKZ;
            }
          }

          if constexpr (calculateGlobals) {
            const auto potentialEnergy3 = factor * (allDistsSquared - _threeDoubleVec * allDotProducts);

            // sum potential energy only on owned particles
            const auto ownedMaskI =
                highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                            highway::Set(tag_double, static_cast<double>(ownedStateI)));
            const auto ownedMaskJ =
                highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                            highway::Set(tag_double, static_cast<double>(ownedStateJ)));

            const auto potentialEnergyI = highway::IfThenElse(ownedMaskI, potentialEnergy3, highway::Zero(tag_double));
            const auto potentialEnergyJ = highway::IfThenElse(ownedMaskJ, potentialEnergy3, highway::Zero(tag_double));

            potentialEnergySum += potentialEnergyI + potentialEnergyJ;

            // Virial for i
            const auto virialIX = forceIX * (distXKIVec - distXIJVec);
            const auto virialIY = forceIY * (distYKIVec - distYIJVec);
            const auto virialIZ = forceIZ * (distZKIVec - distZIJVec);

            const auto maskedVirialIX = highway::IfThenElse(ownedMaskI, virialIX, highway::Zero(tag_double));
            const auto maskedVirialIY = highway::IfThenElse(ownedMaskI, virialIY, highway::Zero(tag_double));
            const auto maskedVirialIZ = highway::IfThenElse(ownedMaskI, virialIZ, highway::Zero(tag_double));

            // Virial for j
            const auto virialJX = forceJX * (distXIJVec - distXJKVec);
            const auto virialJY = forceJY * (distYIJVec - distYJKVec);
            const auto virialJZ = forceJZ * (distZIJVec - distZJKVec);

            const auto maskedVirialJX = highway::IfThenElse(ownedMaskJ, virialJX, highway::Zero(tag_double));
            const auto maskedVirialJY = highway::IfThenElse(ownedMaskJ, virialJY, highway::Zero(tag_double));
            const auto maskedVirialJZ = highway::IfThenElse(ownedMaskJ, virialJZ, highway::Zero(tag_double));

            virialSumX += maskedVirialIX + maskedVirialJX;
            virialSumY += maskedVirialIY + maskedVirialJY;
            virialSumZ += maskedVirialIZ + maskedVirialJZ;

            if constexpr (countFLOPs) {
              numGlobalCalcsN3Sum += remainder ? restK : _vecLengthDouble;
            }
          }
        }
        fxptr1[j] += highway::ReduceSum(tag_double, fXAccJ);
        fyptr1[j] += highway::ReduceSum(tag_double, fYAccJ);
        fzptr1[j] += highway::ReduceSum(tag_double, fZAccJ);
      }
      fxptr1[i] += highway::ReduceSum(tag_double, fXAccI);
      fyptr1[i] += highway::ReduceSum(tag_double, fYAccI);
      fzptr1[i] += highway::ReduceSum(tag_double, fZAccI);
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
    autopas::utils::ExceptionHandler::exception(
        "AxilrodTellerMutoFunctorHWY::SoAFunctorVerletImpl() is not implemented.");
  }

  template <bool newton3, bool remainder, bool alignedSoAView>
  HWY_INLINE void handleKLoopBody(
      const size_t i, const size_t j, const size_t k, const SoAFloatPrecision const_nu, const DistanceMatrix &soaDists1,
      const DistanceMatrix &soaDists2, const autopas::OwnershipState *const __restrict ownedStatePtr,
      const size_t *const __restrict typePtr, const VectorDouble &distXIJVec, const VectorDouble &distYIJVec,
      const VectorDouble &distZIJVec, const VectorDouble &distSquaredIJVec, const VectorDouble &invR5IJVec,
      VectorDouble &fXAccI, VectorDouble &fYAccI, VectorDouble &fZAccI, VectorDouble &fXAccJ, VectorDouble &fYAccJ,
      VectorDouble &fZAccJ, const autopas::OwnershipState ownedStateI, const autopas::OwnershipState ownedStateJ,
      double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr,
      VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &potentialEnergySum,
      size_t &numKernelCallsN3Sum, size_t &numGlobalCalcsSum, size_t restK = 0) {
    const auto [distXJKVec, distYJKVec, distZJKVec, distSquaredJKVec, invR5JKVec] =
        soaDists1.loadRowVec(j, k, restK, remainder);

    const auto [distXKIVecNeg, distYKIVecNeg, distZKIVecNeg, distSquaredKIVec, invR5KIVec] =
        soaDists2.loadRowVec(i, k, restK, remainder);

    const auto distXKIVec = highway::Neg(distXKIVecNeg);
    const auto distYKIVec = highway::Neg(distYKIVecNeg);
    const auto distZKIVec = highway::Neg(distZKIVecNeg);

    const auto mask = highway::FirstN(tag_long, remainder ? restK : _vecLengthDouble);
    const auto ownershipK =
        highway::MaskedLoadOr(highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::dummy)), mask,
                              tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr[k]));

    // calculate masks for cutoff between i<->k and j<->k
    const auto maskJK = highway::Le(distSquaredJKVec, highway::Set(tag_double, _cutoffSquared));
    const auto maskKI = highway::Le(distSquaredKIVec, highway::Set(tag_double, _cutoffSquared));
    // ownership mask
    const auto maskOwnershipK =
        highway::Ne(highway::ConvertTo(tag_double, ownershipK),
                    highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::dummy)));
    // mask for triplets between i and j and multiple k
    const auto maskK = highway::And(highway::And(maskJK, maskKI), maskOwnershipK);

    if (highway::AllFalse(tag_double, maskK)) {
      return;
    }

    // load nus
    auto nu = highway::Set(tag_double, const_nu);
    if constexpr (useMixing) {
      HWY_ALIGN double nus[_vecLengthDouble] = {0.};

      for (size_t n = 0; n < (remainder ? restK : _vecLengthDouble); ++n) {
        nus[n] = _PPLibrary->getMixingNu(typePtr[i], typePtr[j], typePtr[k + n]);
      }
      nu = highway::Load(tag_double, nus);
    }

    VectorDouble forceIX, forceIY, forceIZ;
    VectorDouble forceJX, forceJY, forceJZ;
    VectorDouble factor, allDotProducts, allDistsSquared;
    SoAKernel_HWY<true>(distXIJVec, distYIJVec, distZIJVec, distXJKVec, distYJKVec, distZJKVec, distXKIVec, distYKIVec,
                        distZKIVec, distSquaredIJVec, distSquaredJKVec, distSquaredKIVec, nu, invR5IJVec, invR5JKVec,
                        invR5KIVec, forceIX, forceIY, forceIZ, forceJX, forceJY, forceJZ, factor, allDotProducts,
                        allDistsSquared, maskK);

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
    const VectorDouble fx2 =
        remainder ? highway::LoadN(tag_double, &fxPtr[k], restK) : highway::LoadU(tag_double, &fxPtr[k]);
    const VectorDouble fy2 =
        remainder ? highway::LoadN(tag_double, &fyPtr[k], restK) : highway::LoadU(tag_double, &fyPtr[k]);
    const VectorDouble fz2 =
        remainder ? highway::LoadN(tag_double, &fzPtr[k], restK) : highway::LoadU(tag_double, &fzPtr[k]);

    const VectorDouble fx2New = fx2 + forceKX;
    const VectorDouble fy2New = fy2 + forceKY;
    const VectorDouble fz2New = fz2 + forceKZ;

    storePacked<alignedSoAView>(tag_double, &fxPtr[k], fx2New, remainder ? restK : 0);
    storePacked<alignedSoAView>(tag_double, &fyPtr[k], fy2New, remainder ? restK : 0);
    storePacked<alignedSoAView>(tag_double, &fzPtr[k], fz2New, remainder ? restK : 0);

    if constexpr (countFLOPs) {
      numKernelCallsN3Sum += remainder ? restK : _vecLengthDouble;
    }

    if constexpr (calculateGlobals) {
      const auto potentialEnergy3 = factor * (allDistsSquared - _threeDoubleVec * allDotProducts);

      // sum potential energy only on owned particles
      const auto ownedMaskI = highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                                          highway::Set(tag_double, static_cast<double>(ownedStateI)));
      const auto ownedMaskJ = highway::Eq(highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)),
                                          highway::Set(tag_double, static_cast<double>(ownedStateJ)));
      const auto ownedMaskK =
          highway::Eq(highway::ConvertTo(tag_double, ownershipK),
                      highway::Set(tag_double, static_cast<double>(autopas::OwnershipState::owned)));

      // potential energy for i,j,k (masked)
      const auto potentialEnergyI = highway::IfThenElse(ownedMaskI, potentialEnergy3, highway::Zero(tag_double));
      const auto potentialEnergyJ = highway::IfThenElse(ownedMaskJ, potentialEnergy3, highway::Zero(tag_double));
      const auto potentialEnergyK = highway::IfThenElse(ownedMaskK, potentialEnergy3, highway::Zero(tag_double));

      potentialEnergySum += potentialEnergyI + potentialEnergyJ + potentialEnergyK;

      // Virial for i
      const auto virialIX = forceIX * (distXKIVec - distXIJVec);
      const auto virialIY = forceIY * (distYKIVec - distYIJVec);
      const auto virialIZ = forceIZ * (distZKIVec - distZIJVec);

      const auto maskedVirialIX = highway::IfThenElse(ownedMaskI, virialIX, highway::Zero(tag_double));
      const auto maskedVirialIY = highway::IfThenElse(ownedMaskI, virialIY, highway::Zero(tag_double));
      const auto maskedVirialIZ = highway::IfThenElse(ownedMaskI, virialIZ, highway::Zero(tag_double));

      // Virial for j
      const auto virialJX = forceJX * (distXIJVec - distXJKVec);
      const auto virialJY = forceJY * (distYIJVec - distYJKVec);
      const auto virialJZ = forceJZ * (distZIJVec - distZJKVec);

      const auto maskedVirialJX = highway::IfThenElse(ownedMaskJ, virialJX, highway::Zero(tag_double));
      const auto maskedVirialJY = highway::IfThenElse(ownedMaskJ, virialJY, highway::Zero(tag_double));
      const auto maskedVirialJZ = highway::IfThenElse(ownedMaskJ, virialJZ, highway::Zero(tag_double));

      // Virial for k
      const auto virialKX = forceKX * (distXJKVec - distXKIVec);
      const auto virialKY = forceKY * (distYJKVec - distYKIVec);
      const auto virialKZ = forceKZ * (distZJKVec - distZKIVec);

      const auto maskedVirialKX = highway::IfThenElse(ownedMaskK, virialKX, highway::Zero(tag_double));
      const auto maskedVirialKY = highway::IfThenElse(ownedMaskK, virialKY, highway::Zero(tag_double));
      const auto maskedVirialKZ = highway::IfThenElse(ownedMaskK, virialKZ, highway::Zero(tag_double));

      // Reduce the virial
      virialSumX += maskedVirialIX + maskedVirialJX + maskedVirialKX;
      virialSumY += maskedVirialIY + maskedVirialJY + maskedVirialKY;
      virialSumZ += maskedVirialIZ + maskedVirialJZ + maskedVirialKZ;

      if constexpr (countFLOPs) {
        numGlobalCalcsSum += remainder ? restK : _vecLengthDouble;
      }
    }
  }

  /**
   * Inline helper to compute force components for particle I.
   * Returns tuple of (forceIX, forceIY, forceIZ, factor, allDotProducts, allDistsSquared)
   */
  HWY_INLINE void SoAKernelNoN3(const SoAFloatPrecision &distXIJ, const SoAFloatPrecision &distYIJ,
                                const SoAFloatPrecision &distZIJ, const SoAFloatPrecision &distXJK,
                                const SoAFloatPrecision &distYJK, const SoAFloatPrecision &distZJK,
                                const SoAFloatPrecision &distXKI, const SoAFloatPrecision &distYKI,
                                const SoAFloatPrecision &distZKI, const SoAFloatPrecision &distSquaredIJ,
                                const SoAFloatPrecision &distSquaredJK, const SoAFloatPrecision &distSquaredKI,
                                const SoAFloatPrecision &nu, SoAFloatPrecision &forceIX, SoAFloatPrecision &forceIY,
                                SoAFloatPrecision &forceIZ, SoAFloatPrecision &factor,
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
  HWY_INLINE auto SoAKernelN3(const SoAFloatPrecision &distXIJ, const SoAFloatPrecision &distYIJ,
                              const SoAFloatPrecision &distZIJ, const SoAFloatPrecision &distXJK,
                              const SoAFloatPrecision &distYJK, const SoAFloatPrecision &distZJK,
                              const SoAFloatPrecision &distXKI, const SoAFloatPrecision &distYKI,
                              const SoAFloatPrecision &distZKI, const SoAFloatPrecision &distSquaredIJ,
                              const SoAFloatPrecision &distSquaredJK, const SoAFloatPrecision &distSquaredKI,
                              const SoAFloatPrecision &nu, SoAFloatPrecision &forceIX, SoAFloatPrecision &forceIY,
                              SoAFloatPrecision &forceIZ, SoAFloatPrecision &forceJX, SoAFloatPrecision &forceJY,
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
   * Inline helper to compute force components for particle I.
   * Returns tuple of (forceIX, forceIY, forceIZ, factor, allDotProducts, allDistsSquared)
   */
  template <bool newton3>
  HWY_INLINE auto SoAKernel_HWY(const VectorDouble &distXIJ, const VectorDouble &distYIJ, const VectorDouble &distZIJ,
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

    // Use precomputed 1/r^5 for each pair: allDistsTo5 = (1/r_ij^5)*(1/r_jk^5)*(1/r_ki^5)
    const VectorDouble allDistsTo5 = invR5IJ * invR5JK * invR5KI;

    // factor is three * nu * product(inv r^5)
    factor = _threeDoubleVec * nu * allDistsTo5;

    // Dot products
    const VectorDouble IJDotKI = fmaHelper(distXIJ, distXKI, distYIJ, distYKI, distZIJ, distZKI);
    const VectorDouble IJDotJK = fmaHelper(distXIJ, distXJK, distYIJ, distYJK, distZIJ, distZJK);
    const VectorDouble JKDotKI = fmaHelper(distXJK, distXKI, distYJK, distYKI, distZJK, distZKI);
    allDotProducts = IJDotKI * IJDotJK * JKDotKI;

    // Force I components
    const VectorDouble factorIDirectionJK = factor * IJDotKI * (IJDotJK - JKDotKI);
    const VectorDouble factorIDirectionIJ =
        factor * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI + _fiveDoubleVec * allDotProducts / distSquaredIJ);
    const VectorDouble factorIDirectionKI = factor * (highway::Neg(IJDotJK) * JKDotKI + distSquaredIJ * distSquaredJK -
                                                      _fiveDoubleVec * allDotProducts / distSquaredKI);

    forceIX = fmaHelper(distXJK, factorIDirectionJK, distXIJ, factorIDirectionIJ, distXKI, factorIDirectionKI);
    forceIY = fmaHelper(distYJK, factorIDirectionJK, distYIJ, factorIDirectionIJ, distYKI, factorIDirectionKI);
    forceIZ = fmaHelper(distZJK, factorIDirectionJK, distZIJ, factorIDirectionIJ, distZKI, factorIDirectionKI);

    if constexpr (newton3) {
      // Force J components
      const VectorDouble factorJDirectionKI = factor * IJDotJK * (JKDotKI - IJDotKI);
      const VectorDouble factorJDirectionIJ =
          factor * (highway::Neg(IJDotKI) * JKDotKI + distSquaredJK * distSquaredKI -
                    _fiveDoubleVec * allDotProducts / distSquaredIJ);
      const VectorDouble factorJDirectionJK = factor * (IJDotKI * JKDotKI - distSquaredIJ * distSquaredKI +
                                                        _fiveDoubleVec * allDotProducts / distSquaredJK);

      forceJX = fmaHelper(distXKI, factorJDirectionKI, distXIJ, factorJDirectionIJ, distXJK, factorJDirectionJK);
      forceJY = fmaHelper(distYKI, factorJDirectionKI, distYIJ, factorJDirectionIJ, distYJK, factorJDirectionJK);
      forceJZ = fmaHelper(distZKI, factorJDirectionKI, distZIJ, factorJDirectionIJ, distZJK, factorJDirectionJK);
    }

    // apply mask
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

  template <bool alignedSoAView, class Tag, typename T>
  HWY_INLINE static auto loadPacked(const Tag tag, const T *HWY_RESTRICT ptr, size_t numberToLoad = 0) {
    if (numberToLoad > 0) {
      // partial / remainder load
      return highway::LoadN(tag, ptr, numberToLoad);
    } else {
      // full load
      if constexpr (alignedSoAView) {
        return highway::Load(tag, ptr);
      } else {
        return highway::LoadU(tag, ptr);
      }
    }
  }

  template <bool alignedSoAView, class Tag, typename T, typename Vec>
  HWY_INLINE static void storePacked(const Tag tag, T *HWY_RESTRICT ptr, const Vec &v, size_t numberToStore = 0) {
    if (numberToStore > 0) {
      highway::StoreN(v, tag, ptr, numberToStore);
    } else {
      if constexpr (alignedSoAView) {
        highway::Store(v, tag, ptr);
      } else {
        highway::StoreU(v, tag, ptr);
      }
    }
  }

  template <class SoA, auto AttributeX, auto AttributeY, auto AttributeZ>
  HWY_INLINE bool isAligned3D(const SoA &soa, size_t alignment) {
    return soa.template isAligned<AttributeX>(alignment) && soa.template isAligned<AttributeY>(alignment) &&
           soa.template isAligned<AttributeZ>(alignment);
  }

  struct DistanceMatrix {
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> x;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> y;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> z;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> squared;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> invR5;

    size_t nRows{};
    size_t nCols{};
    bool lowerTriangle{true};

    DistanceMatrix(size_t rows, size_t cols, bool tri = true) : nRows(rows), nCols(cols), lowerTriangle(tri) {
      size_t size = (tri ? (rows * rows - rows) / 2 : rows * cols);
      x.resize(size);
      y.resize(size);
      z.resize(size);
      squared.resize(size);
      invR5.resize(size);
    }

    // Lower-triangle index
    size_t triIndex(size_t i, size_t j) const {
      assert(lowerTriangle and i > j);

      return i * (i - 1) / 2 + j;
    }

    // Fullmatrix-Index
    size_t fullIndex(size_t i, size_t j) const {
      assert(!lowerTriangle);
      return j + i * nCols;
    }

    size_t index(size_t i, size_t j) const { return lowerTriangle ? triIndex(i, j) : fullIndex(i, j); }

    void set(size_t i, size_t j, SoAFloatPrecision dx, SoAFloatPrecision dy, SoAFloatPrecision dz,
             SoAFloatPrecision dist2, SoAFloatPrecision invr5) {
      const size_t idx = index(i, j);
      x[idx] = dx;
      y[idx] = dy;
      z[idx] = dz;
      squared[idx] = dist2;
      invR5[idx] = invr5;
    }

    auto get(size_t i, size_t j) const {
      const size_t idx = index(i, j);
      return std::tuple{x[idx], y[idx], z[idx], squared[idx], invR5[idx]};
    }

    // Highway-Vector load
    auto loadRowVec(size_t i, size_t jStart, size_t width, bool remainder) const {
      size_t idx;

      if (lowerTriangle) {
        if (i > jStart) {
          idx = triIndex(i, jStart);
        } else {
          idx = triIndex(jStart, i);
        }
      } else {
        idx = fullIndex(i, jStart);
      }

      const auto dx = remainder ? highway::LoadN(tag_double, &x[idx], width) : highway::LoadU(tag_double, &x[idx]);
      const auto dy = remainder ? highway::LoadN(tag_double, &y[idx], width) : highway::LoadU(tag_double, &y[idx]);
      const auto dz = remainder ? highway::LoadN(tag_double, &z[idx], width) : highway::LoadU(tag_double, &z[idx]);
      const auto r2 =
          remainder ? highway::LoadN(tag_double, &squared[idx], width) : highway::LoadU(tag_double, &squared[idx]);
      const auto vInvR5 =
          remainder ? highway::LoadN(tag_double, &invR5[idx], width) : highway::LoadU(tag_double, &invR5[idx]);

      return std::tuple{dx, dy, dz, r2, vInvR5};
    }
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
  const VectorDouble _cutoffSquaredVec{};

  // Parameter of the Axilrod-Teller-Muto potential
  // not const because they might be reset through PPL
  double _nu = 0.0;
  VectorDouble _nuVec{highway::Zero(tag_double)};

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

  const VectorDouble _threeDoubleVec = highway::Set(tag_double, 3.0);
  const VectorDouble _fiveDoubleVec = highway::Set(tag_double, 5.0);

  // Alignment for SoAFloatPrecision vector register
  static constexpr std::size_t alignmentSoAFloatHwyVector = highway::Lanes(tag_double) * sizeof(SoAFloatPrecision);
};
}  // namespace mdLib

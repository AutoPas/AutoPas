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
 */
template <class Particle, bool useMixing = false, bool useLUT = false,autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false>
class AxilrodTellerFunctor
    : public autopas::TriwiseFunctor<Particle,
                                     AxilrodTellerFunctor<Particle, useMixing, useLUT, useNewton3, calculateGlobals>> {
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
      : autopas::TriwiseFunctor<Particle, AxilrodTellerFunctor<Particle, useMixing, useLUT, useNewton3, calculateGlobals>>(
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
    // LUT beginA

    const std::array<double, 3> forceI = {0., 0., 0.};
    std::array<double, 3> forceJ;
    std::array<double, 3> forceK;
    double factor;
    double allDistsSquared;
    double allDotProducts;
    double potentialEnergy = 0.;


    if constexpr (useLUT) {
      AutoPasLog(DEBUG, "Used LUT with {} | {} | {}", displacementIJ, displacementJK, displacementKI);
      auto values = _PPLibrary->getATLUT().retrieveValue(displacementIJ, displacementJK, displacementKI);
      forceI = values.first.at(0);
      forceJ = values.first.at(1);
      forceK = values.first.at(2);
      potentialEnergy = values.second;
    }
    else {
      const double IJDotKI = autopas::utils::ArrayMath::dot(displacementIJ, displacementKI);
      const double IJDotJK = autopas::utils::ArrayMath::dot(displacementIJ, displacementJK);
      const double JKDotKI = autopas::utils::ArrayMath::dot(displacementJK, displacementKI);
      allDotProducts = IJDotKI * IJDotJK * JKDotKI;

      allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
      const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
      factor = 3.0 * nu / allDistsTo5;

      const auto forceIDirectionJK = displacementJK * IJDotKI * (IJDotJK - JKDotKI);  // return from LUT (* factor)
      const auto forceIDirectionIJ =
          displacementIJ * (IJDotJK * JKDotKI - distSquaredJK * distSquaredKI +
                            5.0 * allDotProducts / distSquaredIJ);  // return from LUT (* factor)
      const auto forceIDirectionKI =
          displacementKI * (-IJDotJK * JKDotKI + distSquaredIJ * distSquaredJK -
                            5.0 * allDotProducts / distSquaredKI);  // return from LUT (* factor)

      // factor also needs to be included in LUT
      forceI = (forceIDirectionJK + forceIDirectionIJ + forceIDirectionKI) * factor;
      i.addF(forceI);

      // Easy part, what else do we need?
      // -> forceI, J, K
      // Directions?

      forceJ = forceI;
      forceK = forceI;
      if (newton3) {
        const auto forceJDirectionKI = displacementKI * IJDotJK * (JKDotKI - IJDotKI);
        const auto forceJDirectionIJ = displacementIJ * (-IJDotKI * JKDotKI + distSquaredJK * distSquaredKI -
                                                         5.0 * allDotProducts / distSquaredIJ);
        const auto forceJDirectionJK =
            displacementJK * (IJDotKI * JKDotKI - distSquaredIJ * distSquaredKI + 5.0 * allDotProducts / distSquaredJK);

        forceJ = (forceJDirectionKI + forceJDirectionIJ + forceJDirectionJK) * factor;
        j.addF(forceJ);

        forceK = (forceI + forceJ) * (-1.0);
        k.addF(forceK);
      }
    }

    if constexpr (calculateGlobals) {
      // Calculate third of total potential energy from 3-body interaction
      if constexpr (!useLUT) {
        potentialEnergy = factor * (allDistsSquared - 3.0 * allDotProducts) / 9.0;
      }

      // Virial is calculated as f_i * r_i
      // see Thompson et al.: https://doi.org/10.1063/1.3245303
      const auto virialI = forceI * i.getR();
      const int threadnum = autopas::autopas_get_thread_num();
      if (i.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialI;
      }
      // for non-newton3 particles j and/or k will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        const auto virialJ = forceJ * j.getR();
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialJ;
      }
      if (newton3 and k.isOwned()) {
        const auto virialK = forceK * k.getR();
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy;
        _aosThreadData[threadnum].virialSum += virialK;
      }
    }
  }

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

  // Parameter of the Axilrod-Teller potential
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

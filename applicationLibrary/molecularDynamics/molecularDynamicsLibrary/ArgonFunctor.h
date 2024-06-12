/**
* @file ArgonFunctor.h
* @author I. Angelucci
* @date 11/06/24
*/

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib {

namespace {
  size_t indexRepulsivePart(const size_t i, const size_t j, const size_t k) {
    if (i == 0 && j == 0 && k == 0) { return 0; }
    if (i == 0 && j == 0 && k == 1) { return 1; }
    if (i == 0 && j == 1 && k == 1) { return 2; }
    if (i == 1 && j == 1 && k == 1) { return 3; }
    if (i == 0 && j == 0 && k == 2) { return 4; }
    if (i == 0 && j == 1 && k == 2) { return 5; }
    if (i == 1 && j == 1 && k == 2) { return 6; }
    if (i == 0 && j == 2 && k == 2) { return 7; }
    if (i == 1 && j == 2 && k == 2) { return 8; }
    if (i == 2 && j == 2 && k == 2) { return 9; }
    if (i == 0 && j == 0 && k == 3) { return 10; }
    if (i == 0 && j == 1 && k == 3) { return 11; }
    if (i == 1 && j == 1 && k == 3) { return 12; }
    if (i == 0 && j == 2 && k == 3) { return 13; }
    if (i == 1 && j == 2 && k == 3) { return 14; }
    if (i == 0 && j == 3 && k == 3) { return 15; }
    if (i == 0 && j == 0 && k == 4) { return 16; }
    if (i == 0 && j == 1 && k == 4) { return 17; }
    if (i == 1 && j == 1 && k == 4) { return 18; }
    if (i == 0 && j == 2 && k == 4) { return 19; }
    if (i == 0 && j == 0 && k == 5) { return 20; }
    if (i == 0 && j == 1 && k == 5) { return 21; }
    if (i == 0 && j == 0 && k == 6) { return 22; }
    throw autopas::utils::ExceptionHandler::AutoPasException({});
  }

  size_t indexDispersionPart(const size_t i, const size_t j, const size_t k) {
    if (i == 1 && j == 1 && k == 1) { return 0; }
    if (i == 1 && j == 1 && k == 2) { return 1; }
    if (i == 1 && j == 2 && k == 2) { return 2; }
    if (i == 2 && j == 2 && k == 2) { return 3; }
    if (i == 1 && j == 1 && k == 3) { return 4; }
    throw autopas::utils::ExceptionHandler::AutoPasException({});
  }

  enum param{A, alpha, Z, beta};

  template<param P>
  size_t getIndex(const size_t i, const size_t j, const size_t k) {
    size_t index{};
    try {
      if (P == A || P == alpha) {
        index = indexRepulsivePart(i, j, k);
      }
      else if (P == Z || P == beta) {
        index = indexDispersionPart(i, j, k);
      }
    }
    catch (autopas::utils::ExceptionHandler::AutoPasException exception) {
      std::stringstream message;
      message << "Parameter " << P << ": " << "trying to access parameter value at index (i, j, k) = (" << i << ", "
              << j << ", " << k << ").";
      throw autopas::utils::ExceptionHandler::AutoPasException(message.str());
    }
    return index;
  }

} // namespace

template <class Particle, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false>
class ArgonFunctor
    : public autopas::TriwiseFunctor<Particle,
                                     ArgonFunctor<Particle, useMixing, useNewton3, calculateGlobals>> {

  using SoAArraysType = typename Particle::SoAArraysType;
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  ArgonFunctor() = delete;

  /**
   * Constructor for Functor with mixing disabled. When using this functor it is necessary to call
   * setParticleProperties() to set internal constants because it does not use a particle properties library.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit ArgonFunctor(double cutoff) : ArgonFunctor(cutoff, nullptr) {
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
  explicit ArgonFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : ArgonFunctor(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
  }

  std::string getName() final { return "ArgonFunctorAutoVec"; }

  bool isRelevantForTuning() final { return true; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle &i, Particle &j, Particle &k, bool newton3) {}

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
   * @note The molecule types make no difference for ArgonFunctor, but are kept to have a consistent interface
   * for other functors where they may.
   * @return the number of floating point operations
   */
   //TODO @ireneangelucci compute number of flops needed per kernel call, once Functor implementation is completed
  static unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, size_t molCType, bool newton3) {
    return 0;
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
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit ArgonFunctor(double cutoff, void * /*dummy*/)
      : autopas::TriwiseFunctor<Particle, ArgonFunctor<Particle, useMixing, useNewton3, calculateGlobals>>(
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

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    autopas::utils::ExceptionHandler::exception("ArgonFunctor::SoAFunctorVerletImpl() is not yet implemented.");
  }

  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
   //TODO @ireneangelucci Why do we want it to be k*64 Bytes by using __remainingTo64? And not just checking it has k*32 Bytes without storing __remainingTo64?
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    std::array<double, 3> virialSum;
    double potentialEnergySum;

   private:
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  const double _cutoffSquared;

  static constexpr std::array<double, 23> A {{-0.170806873130E+01, //000
                                       -0.316818997395E+02, //001
                                       -0.571545817498E+05, //011
                                       0.848780677578E+02, //111
                                       0.163923794220E+07, //002
                                       0.380809830366E+02, //012
                                       -0.217403993198E+03, //112
                                       0.244997545538E+03, //022
                                       0.128926029735E+03, //122
                                       0.815601247450E+02, //222
                                       0.409987725022E+02, //003
                                       -0.978512983041E+06, //013
                                       0.104383189893E+07, //113
                                       -0.383820796134E+02, //023
                                       0.143934125087E+03, //123
                                       0.102161665959E+04, //033
                                       -0.569593762549E+02, //004
                                       0.178356269769E+04, //014
                                       0.242202158097E+02, //114
                                       -0.279617357863E+01, //024
                                       -0.324585542907E+02, //005
                                       -0.963264559888E-01, //015
                                       -0.898942588279E+05} //006
  };

  static constexpr std::array<double, 23> alpha{{0.428132039316E+00, //000
                                                 0.503934786518E+00, //001
                                                 0.104706730543E+01, //011
                                                 0.456769339560E+00, //111
                                                 0.131047310452E+01, //002
                                                 0.444052360076E+00, //012
                                                 0.480469535570E+00, //112
                                                 0.737327026170E+00, //022
                                                 0.496177745527E+00, //122
                                                 0.424365319847E+00, //222
                                                 0.428946186456E+00, //003
                                                 0.117979281352E+01, //013
                                                 0.119534448663E+01, //113
                                                 0.416753172892E+00, //023
                                                 0.507114743788E+00, //123
                                                 0.764351644551E+00, //033
                                                 0.422619330972E+00, //004
                                                 0.757543022081E+00, //014
                                                 0.482734248672E+00, //114
                                                 0.419340374650E+00, //024
                                                 0.635761316281E+00, //005
                                                 0.375600311119E+00, //015
                                                 0.130334333132E+01} //006
  };

  static constexpr std::array<double, 5> Z{{0.273486414323E+03, //111
                                            -0.213475877256E+05, //112
                                            0.108226781130E+07, //122
                                            -0.213710093072E+07, //222
                                            0.364515182541E+06} //113
  };

  static constexpr std::array<double, 5> beta{{0.211602562917E+02, //111
                                               0.149623190559E+01, //112
                                               0.132161541056E+01, //122
                                               0.208199482789E+01, //222
                                               0.179870559008E+01} //113
  };

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

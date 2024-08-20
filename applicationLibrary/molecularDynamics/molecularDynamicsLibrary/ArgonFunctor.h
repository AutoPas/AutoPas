/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 11/06/24
 */

#pragma once

#include "ArgonInclude/DisplacementHandle.h"
#include "ArgonInclude/RepulsiveTerm.h"
#include "ArgonInclude/DispersionTerm.h"
#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib::Argon {

template <class Particle, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false>
class ArgonFunctor : public autopas::TriwiseFunctor<Particle, ArgonFunctor<Particle, useNewton3, calculateGlobals>> {
  using SoAArraysType = typename Particle::SoAArraysType;
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  ArgonFunctor() = delete;

  /**
   * Constructor of ArgonFunctor
   * @param cutoff
   */
  explicit ArgonFunctor(double cutoff) : ArgonFunctor(cutoff, nullptr) {}

  std::string getName() final { return "ArgonFunctorAutoVec"; }

  bool isRelevantForTuning() final { return true; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  /**
   *
   * @param i particle i
   * @param j particle j
   * @param k paritcle k
   * @param newton3
   */
  void AoSFunctor(Particle &i, Particle &j, Particle &k, bool newton3) {
    const auto A_000{A[index<param::A>(0, 0, 0)]};

    if (i.isDummy() or j.isDummy() or k.isDummy()) {
      return;
    }

    enum IJK { I, J, K };

    const auto displacementIJ = DisplacementHandle(i.getR(), j.getR(), I, J);
    const auto displacementJK = DisplacementHandle(j.getR(), k.getR(), J, K);
    const auto displacementKI = DisplacementHandle(k.getR(), i.getR(), K, I);

    const auto distSquaredIJ = autopas::utils::ArrayMath::dot(displacementIJ.getDisplacement(), displacementIJ.getDisplacement());
    const auto distSquaredJK = autopas::utils::ArrayMath::dot(displacementJK.getDisplacement(), displacementJK.getDisplacement());
    const auto distSquaredKI = autopas::utils::ArrayMath::dot(displacementKI.getDisplacement(), displacementKI.getDisplacement());

    if (distSquaredIJ > _cutoffSquared or distSquaredJK > _cutoffSquared or distSquaredKI > _cutoffSquared) {
      return;
    }

    const auto forceI_repulsive = autopas::utils::ArrayMath::Argon::F_repulsive<I>(A, alpha, displacementIJ, displacementJK, displacementKI);
    const auto forceI_disperisve = autopas::utils::ArrayMath::Argon::F_dispersive<I>(Z, beta, displacementIJ, displacementJK, displacementKI);
    const auto forceI = forceI_repulsive + forceI_disperisve;
    i.addF(forceI);

    const auto forceJ_repulsive = autopas::utils::ArrayMath::Argon::F_repulsive<J>(A, alpha, displacementIJ, displacementJK, displacementKI);
    const auto forceJ_disperisve = autopas::utils::ArrayMath::Argon::F_dispersive<J>(Z, beta, displacementIJ, displacementJK, displacementKI);
    const auto forceJ = forceJ_repulsive + forceJ_disperisve;
    j.addF(forceJ);

    const auto forceK_repulsive = autopas::utils::ArrayMath::Argon::F_repulsive<K>(A, alpha, displacementIJ, displacementJK, displacementKI);
    const auto forceK_disperisve = autopas::utils::ArrayMath::Argon::F_dispersive<K>(Z, beta, displacementIJ, displacementJK, displacementKI);
    const auto forceK = forceK_repulsive + forceK_disperisve;
    k.addF(forceK);
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
  // TODO @ireneangelucci compute number of flops needed per kernel call, once Functor implementation is completed
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
      : autopas::TriwiseFunctor<Particle, ArgonFunctor<Particle, useNewton3, calculateGlobals>>(cutoff),
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
  // TODO @ireneangelucci Why do we want it to be k*64 Bytes by using __remainingTo64? And not just checking it has k*32
  // Bytes without storing __remainingTo64?
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

  static constexpr std::array<double, 23> A{
      {-5.39365445993314E+05,  // 000
       -1.00043526760807E+07,  // 001
       -1.80479894697093E+10,  // 011
       2.68023739515408E+07,   // 111
       5.17630401857978E+11,   // 002
       1.20250233629457E+07,   // 012
       -6.86507513446023E+07,  // 112
       7.73641060191982E+07,   // 022
       4.07116202374599E+07,   // 122
       2.57546504143754E+07,   // 222
       1.29463884038186E+07,   // 003
       -3.08989961490876E+11,  // 013
       3.29616043775900E+11,   // 113
       -1.21201021419532E+07,  // 023
       4.54508019194995E+07,   // 123
       3.22601026022283E+08,   // 033
       -1.79863484497154E+07,  // 004
       5.63204555102674E+08,   // 014
       7.64813924806795E+06,   // 114
       -8.82961781148373E+05,  // 024
       -1.02496007862500E+07,  // 005
       -3.04174890291515E+04,  // 015
       -2.83863618111236E+10}  // 006
  };

  static constexpr std::array<double, 23> alpha{
      {8.09052299484753E+00,  // 000
       9.52298731190775E+00,  // 001
       1.97867044131258E+01,  // 011
       8.63168953894591E+00,  // 111
       2.47643526123088E+01,  // 002
       8.39137345537347E+00,  // 012
       9.07955833453440E+00,  // 112
       1.39334614374608E+01,  // 022
       9.37640048180289E+00,  // 122
       8.01934231300047E+00,  // 222
       8.10590814604500E+00,  // 003
       2.22948530135448E+01,  // 013
       2.25887370431209E+01,  // 113
       7.87549358334693E+00,  // 023
       9.58307979519088E+00,  // 123
       1.44441527110870E+01,  // 033
       7.98634790509653E+00,  // 004
       1.43154883935442E+01,  // 014
       9.12235520967071E+00,  // 114
       7.92438461086463E+00,  // 024
       1.20141476840267E+01,  // 005
       7.09781720339141E+00,  // 015
       2.46296194255218E+01}  // 006
  };

  static constexpr std::array<double, 5> Z{
      {2.81011266190959E-04,   // 111
       -6.14241347348619E-05,  // 112
       8.72021611550457E-06,   // 122
       -4.82191783511564E-08,  // 222
       2.93702828075611E-06}   // 113
  };

  static constexpr std::array<double, 5> beta{
      {3.99870891182023E+02,  // 111
       2.82746852049202E+01,  // 112
       2.49749116804324E+01,  // 122
       3.93440001759947E+01,  // 222
       3.39906094408459E+01}  // 113
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
}  // namespace mdLib::Argon

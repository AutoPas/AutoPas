/**
 * @file AsBuildPairGeneratorFunctor.h
 * @author humig
 * @date 09.07.19
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class Particle>
class VerletNeighborListAsBuild;

namespace internal {

/**
 * This functor can generate or check variable verlet lists using the typical pairwise
 * traversal.
 *
 * @tparam Particle The particle class.
 * @tparam callCheckInstead If false, generate a neighbor list. If true, check the current for validity. Checking
 * validity only works with the AoSFunctor().
 */
template <class Particle, bool callCheckInstead = false>
class AsBuildPairGeneratorFunctor
    : public autopas::Functor<Particle, AsBuildPairGeneratorFunctor<Particle, callCheckInstead>> {
 public:
  /**
   * The structure of the SoAs is defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 4>{
        Particle::AttributeNames::ptr, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 4>{
        Particle::AttributeNames::id, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() { return std::array<typename Particle::AttributeNames, 0>{}; }

  bool allowsNewton3() override { return true; }
  bool allowsNonNewton3() override { return true; }
  bool allowsMixedNewton3() override { return false; }

  /**
   * Constructor of the functor.
   * @param neighborList The neighbor list to fill.
   * @param cutoffskin The cutoff skin to use.
   */
  AsBuildPairGeneratorFunctor(VerletNeighborListAsBuild<Particle> &neighborList, double cutoffskin)
      : autopas::Functor<Particle, AsBuildPairGeneratorFunctor>(cutoffskin),
        _list(neighborList),
        _cutoffskinsquared(cutoffskin * cutoffskin) {}

  [[nodiscard]] bool isRelevantForTuning() override { return false; }

  /**
   * Adds the given pair to the neighbor list.
   * @param i The first particle of the pair.
   * @param j The second particle of the pair.
   * @param newton3 Not used!
   */
  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    using namespace autopas::utils::ArrayMath::literals;

    const auto dist = i.getR() - j.getR();
    const double distsquare = utils::ArrayMath::dot(dist, dist);
    if (distsquare < _cutoffskinsquared) {
      if (callCheckInstead) {
        _list.checkPair(&i, &j);
      } else {
        _list.addPair(&i, &j);
      }
    }
  }

  /**
   * Adds all pairs of the SoA to the neighbor list.
   * @param soa The SoA to add.
   * @param newton3 Whether to use newton 3 or not.
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
    if (soa.getNumberOfParticles() == 0) return;

    auto **const __restrict ptrptr = soa.template begin<Particle::AttributeNames::ptr>();
    double *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    size_t numPart = soa.getNumberOfParticles();
    for (unsigned int i = 0; i < numPart; ++i) {
      for (unsigned int j = i + 1; j < numPart; ++j) {
        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _cutoffskinsquared) {
          _list.addPair(ptrptr[i], ptrptr[j]);
          if (not newton3) {
            _list.addPair(ptrptr[j], ptrptr[i]);
          }
        }
      }
    }
  }

  /**
   * Adds all pairs (p1, p2), p1 element soa1, p2 element soa2 to the neighbor list.
   * @param soa1
   * @param soa2
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool /*newton3*/) override {
    if (soa1.getNumberOfParticles() == 0 || soa2.getNumberOfParticles() == 0) return;

    auto **const __restrict ptrptr1 = soa1.template begin<Particle::AttributeNames::ptr>();
    double *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    double *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    double *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();

    auto **const __restrict ptrptr2 = soa2.template begin<Particle::AttributeNames::ptr>();
    double *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    double *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    double *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    size_t numPart1 = soa1.getNumberOfParticles();
    for (unsigned int i = 0; i < numPart1; ++i) {
      size_t numPart2 = soa2.getNumberOfParticles();
      for (unsigned int j = 0; j < numPart2; ++j) {
        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _cutoffskinsquared) {
          _list.addPair(ptrptr1[i], ptrptr2[j]);
        }
      }
    }
  }

 private:
  /**
   * The neighbor list to fill.
   */
  VerletNeighborListAsBuild<Particle> &_list;
  /**
   * The squared cutoff skin to determine if a pair should be added to the list.
   */
  double _cutoffskinsquared;
};

}  // namespace internal

}  // namespace autopas

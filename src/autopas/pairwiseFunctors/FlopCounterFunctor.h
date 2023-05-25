/**
 * @file FlopCounterFunctor.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
/**
 * This class helps in getting the number of performed floating point
 * operations. It is a functor that only calculated the amount of floating point
 * operations.
 * @todo this class currently is limited to the following case:
 *  - constant cutoff radius
 *  - cases where cutoff is determined with a single distance calculation between two positions ( the notable case where
 *    this doesn't hold is in a CoM-to-site style SoA Functor for multi-site molecules )
 * @todo: we may want the possibility of doing this faster in cases where the number of flops per kernel call is
 * constant
 * @tparam Particle
 * @tparam ForceFunctorType
 */
template <class Particle, class ForceFunctorType>
class FlopCounterFunctor : public Functor<Particle, FlopCounterFunctor<Particle, ForceFunctorType>> {
 public:
  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  /**
   * constructor of FlopCounterFunctor
   * @param forceFunctor force functor whose flops are counted
   * @param cutoffRadius the cutoff radius
   */
  explicit FlopCounterFunctor<Particle, ForceFunctorType>(ForceFunctorType forceFunctor, double cutoffRadius)
      : autopas::Functor<Particle, FlopCounterFunctor<Particle, ForceFunctorType>>(cutoffRadius),
        _forceFunctor(forceFunctor),
        _cutoffSquare(cutoffRadius * cutoffRadius),
        _distanceCalculations(0ul),
        _kernelCalls(0ul),
        _kernelFlops(0ul) {}

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    const auto displacement = i.getR() - j.getR();
    const auto distanceSquared = utils::ArrayMath::dot(displacement, displacement);
    _distanceCalculations.fetch_add(1, std::memory_order_relaxed);

    if (distanceSquared <= _cutoffSquare) {
      _kernelCalls.fetch_add(1, std::memory_order_relaxed);
      _kernelFlops.fetch_add(_forceFunctor.getNumFlopsPerKernelCall(i.getTypeId(), j.getTypeId(), newton3),
                             std::memory_order_relaxed);
    }
  }

  /**
   * @copydoc Functor::SoAFunctorSingle()
   * This SoA Functor does not use any vectorization.
   */
  void SoAFunctorSingle(SoAView<typename Particle::SoAArraysType> soa, bool newton3) override {
    if (soa.getNumberOfParticles() == 0) return;

    double *const __restrict xPtr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict yPtr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict zPtr = soa.template begin<Particle::AttributeNames::posZ>();

    unsigned long *const __restrict typePtr = soa.template begin<Particle::AttributeNames::typeId>();

    for (size_t i = 0; i < soa.getNumberOfParticles(); ++i) {
      for (size_t j = i + 1; j < soa.getNumberOfParticles(); ++j) {
        ++_distanceCalculations;

        const double drx = xPtr[i] - xPtr[j];
        const double dry = yPtr[i] - yPtr[j];
        const double drz = zPtr[i] - zPtr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 <= _cutoffSquare) {
          ++_kernelCalls;
          _kernelFlops += _forceFunctor.getNumFlopsPerKernelCall(typePtr[i], typePtr[j],
                                                                 true);  // SoAFunctorSingle always uses newton3.
        }
      }
    }
  }

  /**
   * @copydoc Functor::SoAFunctorPair()
   */
  void SoAFunctorPair(SoAView<typename Particle::SoAArraysType> soa1, SoAView<typename Particle::SoAArraysType> soa2,
                      bool newton3) override {
    double *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    double *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    double *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    double *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    double *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    double *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    unsigned long *const __restrict type1Ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    unsigned long *const __restrict type2Ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    for (size_t i = 0; i < soa1.getNumberOfParticles(); ++i) {
      for (size_t j = 0; j < soa2.getNumberOfParticles(); ++j) {
        ++_distanceCalculations;

        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 <= _cutoffSquare) {
          ++_kernelCalls;
          _kernelFlops += _forceFunctor.getNumFlopsPerKernelCall(type1Ptr[i], type2Ptr[j], newton3);
        }
      }
    }
  }

  // clang-format off
  /**
   * @copydoc Functor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(SoAView<typename Particle::SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, AlignedAllocator<size_t>> &neighborList, bool newton3) override {
    const auto numParticles = soa.getNumberOfParticles();

    if (numParticles == 0) return;

    double *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    unsigned long *const __restrict typePtr = soa.template begin<Particle::AttributeNames::typeId>();

    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict currentList = neighborList.data();

    for (size_t jNeighIndex = 0; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;

      ++_distanceCalculations;

      const double drx = xptr[indexFirst] - xptr[j];
      const double dry = yptr[indexFirst] - yptr[j];
      const double drz = zptr[indexFirst] - zptr[j];

      const double drx2 = drx * drx;
      const double dry2 = dry * dry;
      const double drz2 = drz * drz;

      const double dr2 = drx2 + dry2 + drz2;

      if (dr2 <= _cutoffSquare) {
        ++_kernelCalls;
        _kernelFlops += _forceFunctor.getNumFlopsPerKernelCall(typePtr[indexFirst], typePtr[j], newton3);
      }
    }
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static std::array<typename Particle::AttributeNames, 3> getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::posX, Particle::AttributeNames::posY, Particle::AttributeNames::posZ};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static std::array<typename Particle::AttributeNames, 3> getNeededAttr(std::false_type) {
    return getNeededAttr();
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static std::array<typename Particle::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 0>{/*Nothing*/};
  }

  /**
   * get the hit rate of the pair-wise interaction, i.e. the ratio of the number
   * of kernel calls compared to the number of distance calculations
   * @return the hit rate
   */
  double getHitRate() { return static_cast<double>(_kernelCalls) / static_cast<double>(_distanceCalculations); }

  /**
   * get the total number of flops
   * @return
   */
  [[nodiscard]] size_t getFlops() const {
    const auto distFlops = numFlopsPerDistanceCalculation * _distanceCalculations;
    return distFlops + _kernelFlops;
  }

  /**
   * get the number of calculated distance operations
   * @return
   */
  [[nodiscard]] size_t getDistanceCalculations() const { return _distanceCalculations; }

  /**
   * get the number of kernel calls, i.e. the number of pairs of particles with
   * a distance not larger than the cutoff
   * @return
   */
  [[nodiscard]] size_t getKernelCalls() const { return _kernelCalls; }

  /**
   * Get the number of kernel flops, i.e. the flops done after the distance calculation and cutoff comparison.
   * @return
   */
  [[nodiscard]] size_t getKernelFlops() const { return _kernelFlops; }

  /**
   * number of flops for one distance calculation.
   * 3 sub + 3 square + 2 add
   */
  static constexpr size_t numFlopsPerDistanceCalculation = 8ul;

 private:
  double _cutoffSquare;
  std::atomic<size_t> _distanceCalculations, _kernelCalls, _kernelFlops;
  ForceFunctorType _forceFunctor;
};

}  // namespace autopas

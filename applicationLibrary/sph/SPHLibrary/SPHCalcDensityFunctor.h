/**
 * @file SPHCalcDensityFunctor.h
 * @author seckler
 * @date 19.01.18
 */

#pragma once

#include "SPHKernels.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"

namespace sphLib {
/**
 * Class that defines the density functor.
 * It is used to calculate the density based on the given SPH kernel.
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle>
class SPHCalcDensityFunctor : public autopas::Functor<Particle, SPHCalcDensityFunctor<Particle>> {
 public:
  /// soa arrays type
  using SoAArraysType = typename Particle::SoAArraysType;

  SPHCalcDensityFunctor() : autopas::Functor<Particle, SPHCalcDensityFunctor<Particle>>(0.){};

  bool isRelevantForTuning() override { return true; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  /**
   * Calculates the density contribution of the interaction of particle i and j.
   * It is not symmetric, because the smoothing lenghts of the two particles can
   * be different.
   * @param i first particle of the interaction
   * @param j second particle of the interaction
   * @param newton3 defines whether or whether not to use newton 3
   */
  inline void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }
    const std::array<double, 3> dr = j.getR() - i.getR();  // ep_j[j].pos - ep_i[i].pos;
    const double density =
        j.getMass() * SPHKernels::W(dr, i.getSmoothingLength());  // ep_j[j].mass * W(dr, ep_i[i].smth)
    i.addDensity(density);
    if (newton3) {
      // Newton 3:
      // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
      const double density2 = i.getMass() * SPHKernels::W(dr, j.getSmoothingLength());
      j.addDensity(density2);
    }
  }

  /**
   * Get the number of floating point operations used in one full kernel call
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall() {
    unsigned long flops = 0;
    flops += 3;                            // calculating dr
    flops += 2 * SPHKernels::getFlopsW();  // flops for calling W
    flops += 2 * 1;                        // calculating density
    flops += 2 * 1;                        // adding density
    return flops;
  }

  /**
   * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType>, bool)
   * This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) override {
    if (soa.getNumberOfParticles() == 0) return;

    double *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict densityptr = soa.template begin<Particle::AttributeNames::density>();
    double *const __restrict smthptr = soa.template begin<Particle::AttributeNames::smth>();
    double *const __restrict massptr = soa.template begin<Particle::AttributeNames::mass>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    size_t numParticles = soa.getNumberOfParticles();
    for (unsigned int i = 0; i < numParticles; ++i) {
      // checks whether particle i is owned.
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      double densacc = 0.;
// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : densacc)
      for (unsigned int j = i + 1; j < numParticles; ++j) {
        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        // if second particle is a dummy, we skip the interaction.
        const bool mask = ownedStatePtr[j] != autopas::OwnershipState::dummy;

        const double density = mask ? massptr[j] * SPHKernels::W(dr2, smthptr[i]) : 0.;
        densacc += density;

        // Newton 3:
        // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
        const double density2 = mask ? massptr[i] * SPHKernels::W(dr2, smthptr[j]) : 0.;
        densityptr[j] += density2;
      }

      densityptr[i] += densacc;
    }
  }

  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType>, SoAView<SoAArraysType>, bool)
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      bool newton3) override {
    if (soa1.getNumberOfParticles() == 0 || soa2.getNumberOfParticles() == 0) return;

    double *const __restrict xptr1 = soa1.template begin<Particle::AttributeNames::posX>();
    double *const __restrict yptr1 = soa1.template begin<Particle::AttributeNames::posY>();
    double *const __restrict zptr1 = soa1.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict densityptr1 = soa1.template begin<Particle::AttributeNames::density>();
    double *const __restrict smthptr1 = soa1.template begin<Particle::AttributeNames::smth>();
    double *const __restrict massptr1 = soa1.template begin<Particle::AttributeNames::mass>();

    double *const __restrict xptr2 = soa2.template begin<Particle::AttributeNames::posX>();
    double *const __restrict yptr2 = soa2.template begin<Particle::AttributeNames::posY>();
    double *const __restrict zptr2 = soa2.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict densityptr2 = soa2.template begin<Particle::AttributeNames::density>();
    double *const __restrict smthptr2 = soa2.template begin<Particle::AttributeNames::smth>();
    double *const __restrict massptr2 = soa2.template begin<Particle::AttributeNames::mass>();

    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

    size_t numParticlesi = soa1.getNumberOfParticles();
    for (unsigned int i = 0; i < numParticlesi; ++i) {
      // checks whether particle i is in the domain box, unused if calculateGlobals is false!
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
        continue;
      }

      double densacc = 0.;
      size_t numParticlesj = soa2.getNumberOfParticles();
// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : densacc)
      for (unsigned int j = 0; j < numParticlesj; ++j) {
        const double drx = xptr1[i] - xptr2[j];
        const double dry = yptr1[i] - yptr2[j];
        const double drz = zptr1[i] - zptr2[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        // if second particle is a dummy, we skip the interaction.
        const bool mask = ownedStatePtr2[j] != autopas::OwnershipState::dummy;

        const double density = mask ? massptr2[j] * SPHKernels::W(dr2, smthptr1[i]) : 0.;
        densacc += density;
        if (newton3) {
          // Newton 3:
          // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
          const double density2 = mask ? massptr1[i] * SPHKernels::W(dr2, smthptr2[j]) : 0.;
          densityptr2[j] += density2;
        }
      }

      densityptr1[i] += densacc;
    }
  }

  // clang-format off
  /**
   * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) override {
    if (soa.getNumberOfParticles() == 0) return;

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    // checks whether particle i is owned.
    if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
      return;
    }

    double *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict densityptr = soa.template begin<Particle::AttributeNames::density>();
    double *const __restrict smthptr = soa.template begin<Particle::AttributeNames::smth>();
    double *const __restrict massptr = soa.template begin<Particle::AttributeNames::mass>();

    double densacc = 0;
    const auto &currentList = neighborList;
    size_t listSize = currentList.size();
// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : densacc)
    for (unsigned int j = 0; j < listSize; ++j) {
      const double drx = xptr[indexFirst] - xptr[currentList[j]];
      const double dry = yptr[indexFirst] - yptr[currentList[j]];
      const double drz = zptr[indexFirst] - zptr[currentList[j]];

      const double drx2 = drx * drx;
      const double dry2 = dry * dry;
      const double drz2 = drz * drz;

      const double dr2 = drx2 + dry2 + drz2;

      // if second particle is a dummy, we skip the interaction.
      const bool mask = ownedStatePtr[currentList[j]] != autopas::OwnershipState::dummy;

      const double density = mask ? massptr[currentList[j]] * SPHKernels::W(dr2, smthptr[indexFirst]) : 0.;
      densacc += density;
      if (newton3) {
        // Newton 3:
        // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
        const double density2 = mask ? massptr[indexFirst] * SPHKernels::W(dr2, smthptr[currentList[j]]) : 0.;
        densityptr[currentList[j]] += density2;
      }
    }

    densityptr[indexFirst] += densacc;
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 7>{
        Particle::AttributeNames::mass,          Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,          Particle::AttributeNames::smth, Particle::AttributeNames::density,
        Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::mass, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::smth, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 1>{Particle::AttributeNames::density};
  }
};
}  // namespace sphLib

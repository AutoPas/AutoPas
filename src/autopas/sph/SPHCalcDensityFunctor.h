/**
 * @file SPHCalcDensityFunctor.h
 * @author seckler
 * @date 19.01.18
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/sph/SPHKernels.h"
#include "autopas/sph/SPHParticle.h"

namespace autopas {
namespace sph {
/**
 * Class that defines the density functor.
 * It is used to calculate the density based on the given SPH kernel.
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle, class ParticleCell>
class SPHCalcDensityFunctor : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType,
                                             SPHCalcDensityFunctor<Particle, ParticleCell>> {
 public:
  /// soa arrays type
  using SoAArraysType = typename Particle::SoAArraysType;

  SPHCalcDensityFunctor() : autopas::Functor<Particle, ParticleCell, SoAArraysType, SPHCalcDensityFunctor>(0.){};

  bool isRelevantForTuning() override { return true; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  bool isAppropriateClusterSize(unsigned int clusterSize, DataLayoutOption::Value dataLayout) const override {
    return dataLayout == DataLayoutOption::aos;  // This functor does only support clusters via aos.
  }

  /**
   * Calculates the density contribution of the interaction of particle i and j.
   * It is not symmetric, because the smoothing lenghts of the two particles can
   * be different.
   * @param i first particle of the interaction
   * @param j second particle of the interaction
   * @param newton3 defines whether or whether not to use newton 3
   */
  inline void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    const std::array<double, 3> dr = utils::ArrayMath::sub(j.getR(), i.getR());  // ep_j[j].pos - ep_i[i].pos;
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
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ densityptr = soa.template begin<Particle::AttributeNames::density>();
    double *const __restrict__ smthptr = soa.template begin<Particle::AttributeNames::smth>();
    double *const __restrict__ massptr = soa.template begin<Particle::AttributeNames::mass>();
    size_t numParticles = soa.getNumParticles();
    for (unsigned int i = 0; i < numParticles; ++i) {
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

        const double density = massptr[j] * SPHKernels::W(dr2, smthptr[i]);
        densacc += density;

        // Newton 3:
        // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
        const double density2 = massptr[i] * SPHKernels::W(dr2, smthptr[j]);
        densityptr[j] += density2;
      }

      densityptr[i] += densacc;
    }
  }

  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType>, SoAView<SoAArraysType>, bool)
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3) override {
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

    double *const __restrict__ xptr1 = soa1.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr1 = soa1.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr1 = soa1.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ densityptr1 = soa1.template begin<Particle::AttributeNames::density>();
    double *const __restrict__ smthptr1 = soa1.template begin<Particle::AttributeNames::smth>();
    double *const __restrict__ massptr1 = soa1.template begin<Particle::AttributeNames::mass>();

    double *const __restrict__ xptr2 = soa2.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr2 = soa2.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr2 = soa2.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ densityptr2 = soa2.template begin<Particle::AttributeNames::density>();
    double *const __restrict__ smthptr2 = soa2.template begin<Particle::AttributeNames::smth>();
    double *const __restrict__ massptr2 = soa2.template begin<Particle::AttributeNames::mass>();

    size_t numParticlesi = soa1.getNumParticles();
    for (unsigned int i = 0; i < numParticlesi; ++i) {
      double densacc = 0.;
      size_t numParticlesj = soa2.getNumParticles();
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

        const double density = massptr2[j] * SPHKernels::W(dr2, smthptr1[i]);
        densacc += density;
        if (newton3) {
          // Newton 3:
          // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
          const double density2 = massptr1[i] * SPHKernels::W(dr2, smthptr2[j]);
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
  void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) override {
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ densityptr = soa.template begin<Particle::AttributeNames::density>();
    double *const __restrict__ smthptr = soa.template begin<Particle::AttributeNames::smth>();
    double *const __restrict__ massptr = soa.template begin<Particle::AttributeNames::mass>();

    double densacc = 0;
    auto &currentList = neighborList;
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

      const double density = massptr[currentList[j]] * SPHKernels::W(dr2, smthptr[indexFirst]);
      densacc += density;
      if (newton3) {
        // Newton 3:
        // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
        const double density2 = massptr[indexFirst] * SPHKernels::W(dr2, smthptr[currentList[j]]);
        densityptr[currentList[j]] += density2;
      }
    }

    densityptr[indexFirst] += densacc;
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static std::array<typename SPHParticle::AttributeNames, 6> getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::mass, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::smth, Particle::AttributeNames::density};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static std::array<typename SPHParticle::AttributeNames, 5> getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 5>{
        Particle::AttributeNames::mass, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::smth};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static std::array<typename SPHParticle::AttributeNames, 1> getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 1>{Particle::AttributeNames::density};
  }
};
}  // namespace sph
}  // namespace autopas

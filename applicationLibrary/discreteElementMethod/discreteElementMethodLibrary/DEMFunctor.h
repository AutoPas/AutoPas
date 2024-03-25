/**
 * @file DEMFunctor.h
 * @author Hoppe (hoppef@hsu-hh.de)
 * @brief Functor calculating the particle interaction force.
 * @version 0.1
 * @date 2022-11-11
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <array>
#include <cmath>

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace demLib {

/**
 * @brief Functor for the force interactions between two DEM objects (particles)
 *
 * @tparam Particle
 * @tparam useNewton3
 * @tparam relevantForTuning
 */
template <class Particle, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool relevantForTuning = true>
class DEMFunctor : public autopas::Functor<Particle, DEMFunctor<Particle, useNewton3, relevantForTuning>> {
  /**
   * @brief Structure of the SoAs defined by the particle.
   *
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * @brief Precision of SoA entries.
   *
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * @brief Delete the default constructor.
   *
   */
  DEMFunctor() = delete;

  /**
   * @brief Actual constructor a new DEMFunctor object
   *
   * @param cutoff
   */
  explicit DEMFunctor(double cutoff)

      : autopas::Functor<Particle, DEMFunctor<Particle, useNewton3, relevantForTuning>>(cutoff) {}

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  /**
   * @brief DEM Functor for an AoS
   *
   * Hertz collision between two particles i and j with the contact force
   *
   * F = 4/3*E*sqrt(R)*sqrt(delta^3)
   * R = R1*R2/(R1+R2)
   * 1/E = (1-v1^2)/E1 + (1-v2^2)/E2
   *
   * with delta the penetration depth, R1 and R2 the particle radii,
   * E1 and E2 the Young's moduli and v1 and v2 the Poisson ratios of the colliding particles.
   *
   * @param i Particle i
   * @param j Particle j
   * @param newton3 Use Newton's third law
   */
  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    if (i.isDummy() or j.isDummy()) {
      return;
    }

    // Get particle radii
    const double radI = i.getRad();
    const double radJ = j.getRad();

    // distance between ParticleCenters
    const auto dr = autopas::utils::ArrayMath::sub(i.getR(), j.getR());
    const double penDepth = radI + radJ - autopas::utils::ArrayMath::L2Norm(dr);

    // relative particle velocity
    const auto relVel = autopas::utils::ArrayMath::sub(i.getV(), j.getV());

    // Particles do not intersect => return
    if (penDepth <= 0) {
      return;
    }

    // Calculate equivalent radius
    const double radEq = radI * radJ / (radI + radJ);

    // get the particle mass
    const double massI = i.getMass();
    const double massJ = j.getMass();

    // Calculate the equivalent mass
    const double massEq = massI * massJ / (massI + massJ);

    // get the Poisson ratios
    const double poissonI = i.getPoisson();
    const double poissonJ = j.getPoisson();

    // get the Young's moduli
    const double youngI = i.getYoung();
    const double youngJ = j.getYoung();

    const double effectiveYoungI = (1.0 - poissonI * poissonI) / youngI;
    const double effectiveYoungJ = (1.0 - poissonJ * poissonJ) / youngJ;

    // Calculate the equivalent Young's modulus
    const double youngEq = 1.0 / (effectiveYoungI + effectiveYoungJ);

    // Calculate the normal force factor and damping factor
    const double normalFactorForce = 4. / 3. * youngEq * sqrt(radEq * penDepth);

    const double normalFactorDamping = _dampingRatio * 2.0 * sqrt(massEq * normalFactorForce);

    // Calculate the normal force vector
    const double normalForce = normalFactorForce * penDepth;

    const auto vecNormalForce =
        autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::normalize(dr), normalForce);

    // Calculate the normal damping
    const auto vecNormalDamping = autopas::utils::ArrayMath::mulScalar(relVel, normalFactorDamping);

    // Calculate the force
    auto vecForce = autopas::utils::ArrayMath::sub(vecNormalForce, vecNormalDamping);

    i.addF(vecForce);

    if (newton3) {
      j.subF(vecForce);
    }
  }

  /**
   *
   * @copydoc autopas::Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
   *
   * DEM Functor for an SoA
   *
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    // If SoA is empty, return
    if (soa.size() == 0) return;

    // Pointers to the start of arrays containing relevant particle data
    const auto *const __restrict idptr = soa.template begin<Particle::AttributeNames::id>();

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict vxptr = soa.template begin<Particle::AttributeNames::velocityX>();
    const auto *const __restrict vyptr = soa.template begin<Particle::AttributeNames::velocityY>();
    const auto *const __restrict vzptr = soa.template begin<Particle::AttributeNames::velocityZ>();

    const auto *const __restrict massptr = soa.template begin<Particle::AttributeNames::mass>();

    const auto *const __restrict radptr = soa.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poissonptr = soa.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict youngptr = soa.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownershipPtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    // Go over all particles in the SoA
    for (unsigned int i = 0; i < soa.size(); i++) {
      // If particle is dummy, skip it
      const auto ownedStateI = ownershipPtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      // accumulating Forces
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = i + 1; j < soa.size(); j++) {
        // Calculate pair sizes, positions and distances
        const auto ownedStateJ = ownershipPtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
        const SoAFloatPrecision dr = sqrt(dr2);

        const SoAFloatPrecision radI = radptr[i];
        const SoAFloatPrecision radJ = radptr[j];

        const SoAFloatPrecision rad = radI + radJ;

        // Mask away if particles arent intersecting or if j is dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr <= rad and ownedStateJ != autopas::OwnershipState::dummy;

        // Relative particle velocity
        const SoAFloatPrecision dvx = vxptr[i] - vxptr[j];
        const SoAFloatPrecision dvy = vyptr[i] - vyptr[j];
        const SoAFloatPrecision dvz = vzptr[i] - vzptr[j];

        // Particle masses and equivalent mass
        const SoAFloatPrecision massI = massptr[i];
        const SoAFloatPrecision massJ = massptr[j];

        const SoAFloatPrecision massEq = massI * massJ / (massI + massJ);

        // Young's modulus and Poisson ratio of partice I and J
        const SoAFloatPrecision poissonI = poissonptr[i];
        const SoAFloatPrecision youngI = youngptr[i];
        const SoAFloatPrecision poissonJ = poissonptr[j];
        const SoAFloatPrecision youngJ = youngptr[j];

        const SoAFloatPrecision effectiveYoungI = (1.0 - poissonI * poissonI) / youngI;
        const SoAFloatPrecision effectiveYoungJ = (1.0 - poissonJ * poissonJ) / youngJ;

        // Calculate the equivalent Young's modulus
        const SoAFloatPrecision youngEq = 1.0 / (effectiveYoungI + effectiveYoungJ);

        // Calculate the equivalent radius and the penetration depth
        const SoAFloatPrecision radEq = radI * radJ / (radI + radJ);

        const SoAFloatPrecision penDepth = rad - dr;

        // Normal force factor
        const SoAFloatPrecision normalFactorForce = mask ? 4. / 3. * youngEq * sqrt(radEq * penDepth) : 0.;

        // Normal damping factor
        const SoAFloatPrecision normalFactorDamping = _dampingRatio * 2.0 * sqrt(massEq * normalFactorForce);

        // Calculate the normal force and the normal force vector
        const SoAFloatPrecision normalForce = mask ? normalFactorForce * penDepth / dr : 0.0;
        const SoAFloatPrecision normalDamping = mask ? normalFactorDamping : 0.0;

        const SoAFloatPrecision fx = drx * normalForce - dvx * normalDamping;
        const SoAFloatPrecision fy = dry * normalForce - dvy * normalDamping;
        const SoAFloatPrecision fz = drz * normalForce - dvz * normalDamping;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        // Newton3
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
  }

  /**
   * @copydoc autopas::Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

 private:
  /**
   * Implementation of SoAFunctorPair(soa1, soa2, newton3) for DEM
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    // If SoA is empty, return
    if (soa1.size() == 0 || soa2.size() == 0) return;

    // Pointers to the start of arrays containing relevant particle data
    const auto *const __restrict id1ptr = soa1.template begin<Particle::AttributeNames::id>();
    const auto *const __restrict id2ptr = soa2.template begin<Particle::AttributeNames::id>();

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    SoAFloatPrecision *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    SoAFloatPrecision *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict vx1ptr = soa1.template begin<Particle::AttributeNames::velocityX>();
    const auto *const __restrict vy1ptr = soa1.template begin<Particle::AttributeNames::velocityY>();
    const auto *const __restrict vz1ptr = soa1.template begin<Particle::AttributeNames::velocityZ>();
    const auto *const __restrict vx2ptr = soa2.template begin<Particle::AttributeNames::velocityX>();
    const auto *const __restrict vy2ptr = soa2.template begin<Particle::AttributeNames::velocityY>();
    const auto *const __restrict vz2ptr = soa2.template begin<Particle::AttributeNames::velocityZ>();

    const auto *const __restrict mass1ptr = soa1.template begin<Particle::AttributeNames::mass>();
    const auto *const __restrict mass2ptr = soa2.template begin<Particle::AttributeNames::mass>();

    const auto *const __restrict rad1ptr = soa1.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poisson1ptr = soa1.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict young1ptr = soa1.template begin<Particle::AttributeNames::young>();
    const auto *const __restrict rad2ptr = soa2.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poisson2ptr = soa2.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict young2ptr = soa2.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownership1Ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownership2Ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();

    // Go over all particles in the SoA
    for (unsigned int i = 0; i < soa1.size(); i++) {
      // If particle is dummy, skip it
      const auto ownedStateI = ownership1Ptr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      // accumulating Force directions
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = 0; j < soa2.size(); ++j) {
        const auto ownedStateJ = ownership2Ptr[j];

        // Calculate pair sizes, positions and distances
        const SoAFloatPrecision drx = x1ptr[i] - x2ptr[j];
        const SoAFloatPrecision dry = y1ptr[i] - y2ptr[j];
        const SoAFloatPrecision drz = z1ptr[i] - z2ptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
        const SoAFloatPrecision dr = sqrt(dr2);

        const SoAFloatPrecision radI = rad1ptr[i];
        const SoAFloatPrecision radJ = rad2ptr[j];

        const SoAFloatPrecision rad = radI + radJ;

        // Mask away if particles arent intersecting or if j is dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr <= rad and ownedStateJ != autopas::OwnershipState::dummy;

        // Relative particle velocity
        const SoAFloatPrecision dvx = vx1ptr[i] - vx2ptr[j];
        const SoAFloatPrecision dvy = vy1ptr[i] - vy2ptr[j];
        const SoAFloatPrecision dvz = vz1ptr[i] - vz2ptr[j];

        // Particle masses and equivalent mass
        const SoAFloatPrecision massI = mass1ptr[i];
        const SoAFloatPrecision massJ = mass2ptr[j];

        const SoAFloatPrecision massEq = massI * massJ / (massI + massJ);

        // Young's modulus and Poisson ratio of partice I and J
        const SoAFloatPrecision poissonI = poisson1ptr[i];
        const SoAFloatPrecision youngI = young1ptr[i];
        const SoAFloatPrecision poissonJ = poisson2ptr[j];
        const SoAFloatPrecision youngJ = young2ptr[j];

        const SoAFloatPrecision effectiveYoungI = (1.0 - poissonI * poissonI) / youngI;
        const SoAFloatPrecision effectiveYoungJ = (1.0 - poissonJ * poissonJ) / youngJ;

        // Calculate the equivalent Young's modulus
        const SoAFloatPrecision youngEq = 1.0 / (effectiveYoungI + effectiveYoungJ);

        // Calculate the equivalent radius and the penetration depth
        const SoAFloatPrecision radEq = radI * radJ / (radI + radJ);

        const SoAFloatPrecision penDepth = rad - dr;

        // Normal force factor
        const SoAFloatPrecision normalFactorForce = mask ? 4. / 3. * youngEq * sqrt(radEq * penDepth) : 0.;

        // Normal damping factor
        const SoAFloatPrecision normalFactorDamping = _dampingRatio * 2.0 * sqrt(massEq * normalFactorForce);

        // Calculate the normal force and the normal force vector
        const SoAFloatPrecision normalForce = mask ? normalFactorForce * penDepth / dr : 0.0;
        const SoAFloatPrecision normalDamping = mask ? normalFactorDamping : 0.0;

        const SoAFloatPrecision fx = (drx * normalForce - dvx * normalDamping);
        const SoAFloatPrecision fy = (dry * normalForce - dvy * normalDamping);
        const SoAFloatPrecision fz = (drz * normalForce - dvz * normalDamping);

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        // Newton3
        if (newton3) {
          fx2ptr[j] -= fx;
          fy2ptr[j] -= fy;
          fz2ptr[j] -= fz;
        }
      }

      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
  }

 public:
  // clang-format off
  /**
   * @copydoc autopas::Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   * 
   * DEM Functor for a Verlet SoA
   * 
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    // Return, if SoA or neighbor list empty
    if (soa.size() == 0 or neighborList.empty()) return;

    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 16>{
        Particle::AttributeNames::id,        Particle::AttributeNames::posX,
        Particle::AttributeNames::posY,      Particle::AttributeNames::posZ,
        Particle::AttributeNames::forceX,    Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ,    Particle::AttributeNames::velocityX,
        Particle::AttributeNames::velocityY, Particle::AttributeNames::velocityZ,
        Particle::AttributeNames::mass,      Particle::AttributeNames::rad,
        Particle::AttributeNames::poisson,   Particle::AttributeNames::young,
        Particle::AttributeNames::typeId,    Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 13>{Particle::AttributeNames::id,
                                                             Particle::AttributeNames::posX,
                                                             Particle::AttributeNames::posY,
                                                             Particle::AttributeNames::posZ,
                                                             Particle::AttributeNames::velocityX,
                                                             Particle::AttributeNames::velocityY,
                                                             Particle::AttributeNames::velocityZ,
                                                             Particle::AttributeNames::mass,
                                                             Particle::AttributeNames::rad,
                                                             Particle::AttributeNames::poisson,
                                                             Particle::AttributeNames::young,
                                                             Particle::AttributeNames::typeId,
                                                             Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

  /**
   * Get the number of flops used per kernel call for a given particle pair. This should count the
   * floating point operations needed for two particles that lie within a cutoff radius, having already calculated the
   * distance.
   * @param parAType particle A's type id
   * @param parBType particle B's type id
   * @param newton3 is newton3 applied.
   * @note parAType and parBType make no difference for DEMFunctor, but are kept to have a consistent interface for
   * other functors where they may.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall(size_t parAType, size_t parBType, bool newton3) {
    // Kernel: THIS has to be recalculated
    return newton3 ? 18ul : 15ul;
  }

  void initTraversal() final { _postProcessed = false; }

  void endTraversal(bool newton3) final {
    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
  }

 private:
  /**
   * @brief Implementation of SoAFunctorVerlet
   *
   * @tparam newton3
   * @param soa
   * @param indexFirst
   * @param neighborList
   */
  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    // Pointers to particle properties
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    auto *const __restrict vxptr = soa.template begin<Particle::AttributeNames::velocityX>();
    auto *const __restrict vyptr = soa.template begin<Particle::AttributeNames::velocityY>();
    auto *const __restrict vzptr = soa.template begin<Particle::AttributeNames::velocityZ>();

    auto *const __restrict massptr = soa.template begin<Particle::AttributeNames::mass>();

    auto *const __restrict radptr = soa.template begin<Particle::AttributeNames::rad>();
    auto *const __restrict poissonptr = soa.template begin<Particle::AttributeNames::poisson>();
    auto *const __restrict youngptr = soa.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    SoAFloatPrecision fxacc = 0;
    SoAFloatPrecision fyacc = 0;
    SoAFloatPrecision fzacc = 0;

    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    // checks whether particle i is owned.
    const auto ownedStateI = ownedStatePtr[indexFirst];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }

    // this is a magic number, that should correspond to at least
    // vectorization width*N have testet multiple sizes:
    // 4: does not give a speedup, slower than original AoSFunctor
    // 8: small speedup compared to AoS
    // 12: highest speedup compared to Aos
    // 16: smaller speedup
    // in theory this is a variable, we could auto-tune over...
#ifdef __AVX512F__
    // use a multiple of 8 for avx
    constexpr size_t vecsize = 16;
#else
    // for everything else 12 is faster
    constexpr size_t vecsize = 12;
#endif
    size_t joff = 0;

    // if the size of the verlet list is larger than the given size vecsize,
    // we will use a vectorized version.
    if (neighborListSize >= vecsize) {
      alignas(64) std::array<SoAFloatPrecision, vecsize> xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp, masstmp, radtmp,
          poissontmp, youngtmp, xArr, yArr, zArr, vxArr, vyArr, vzArr, massArr, radArr, poissonArr, youngArr, fxArr,
          fyArr, fzArr;
      alignas(64) std::array<autopas::OwnershipState, vecsize> ownedStateArr{};

      // broadcast of the position of particle i
      for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
        xtmp[tmpj] = xptr[indexFirst];
        ytmp[tmpj] = yptr[indexFirst];
        ztmp[tmpj] = zptr[indexFirst];
        vxtmp[tmpj] = vxptr[indexFirst];
        vytmp[tmpj] = vyptr[indexFirst];
        vztmp[tmpj] = vzptr[indexFirst];
        masstmp[tmpj] = massptr[indexFirst];
        radtmp[tmpj] = radptr[indexFirst];
        poissontmp[tmpj] = poissonptr[indexFirst];
        youngtmp[tmpj] = youngptr[indexFirst];
      }

      // loop over the verlet list from 0 to x*vecsize
      for (; joff < neighborListSize - vecsize + 1; joff += vecsize) {
// in each iteration we calculate the interactions of particle i with
// vecsize particles in the neighborlist of particle i starting at
// particle joff

// gather position of particle j
#pragma omp simd safelen(vecsize)
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          xArr[tmpj] = xptr[neighborListPtr[joff + tmpj]];
          yArr[tmpj] = yptr[neighborListPtr[joff + tmpj]];
          zArr[tmpj] = zptr[neighborListPtr[joff + tmpj]];
          vxArr[tmpj] = vxptr[neighborListPtr[joff + tmpj]];
          vyArr[tmpj] = vyptr[neighborListPtr[joff + tmpj]];
          vzArr[tmpj] = vzptr[neighborListPtr[joff + tmpj]];
          massArr[tmpj] = massptr[neighborListPtr[joff + tmpj]];
          radArr[tmpj] = radptr[neighborListPtr[joff + tmpj]];
          poissonArr[tmpj] = poissonptr[neighborListPtr[joff + tmpj]];
          youngArr[tmpj] = youngptr[neighborListPtr[joff + tmpj]];
          ownedStateArr[tmpj] = ownedStatePtr[neighborListPtr[joff + tmpj]];
        }

        // do omp simd with reduction of the interaction
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc) safelen(vecsize)
        for (size_t j = 0; j < vecsize; j++) {
          const auto ownedStateJ = ownedStateArr[j];

          // Calculate pair sizes, positions and distances
          const SoAFloatPrecision drx = xtmp[j] - xArr[j];
          const SoAFloatPrecision dry = ytmp[j] - yArr[j];
          const SoAFloatPrecision drz = ztmp[j] - zArr[j];

          const SoAFloatPrecision drx2 = drx * drx;
          const SoAFloatPrecision dry2 = dry * dry;
          const SoAFloatPrecision drz2 = drz * drz;

          const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
          const SoAFloatPrecision dr = sqrt(dr2);

          const SoAFloatPrecision radI = radtmp[j];
          const SoAFloatPrecision radJ = radArr[j];

          const SoAFloatPrecision rad = radI + radJ;

          // Mask away if particles arent intersecting or if j is dummy.
          // Particle ownedStateI was already checked previously.
          const bool mask = dr <= rad and ownedStateJ != autopas::OwnershipState::dummy and dr != 0;

          // Relative particle velocity
          const SoAFloatPrecision dvx = vxtmp[j] - vxArr[j];
          const SoAFloatPrecision dvy = vytmp[j] - vyArr[j];
          const SoAFloatPrecision dvz = vztmp[j] - vzArr[j];

          // Particle masses and equivalent mass
          const SoAFloatPrecision massI = masstmp[j];
          const SoAFloatPrecision massJ = massArr[j];

          const SoAFloatPrecision massEq = massI * massJ / (massI + massJ);

          // Young's modulus and Poisson ratio of partice I and J
          const SoAFloatPrecision poissonI = poissontmp[j];
          const SoAFloatPrecision youngI = youngtmp[j];
          const SoAFloatPrecision poissonJ = poissonArr[j];
          const SoAFloatPrecision youngJ = youngArr[j];

          const SoAFloatPrecision effectiveYoungI = (1.0 - poissonI * poissonI) / youngI;
          const SoAFloatPrecision effectiveYoungJ = (1.0 - poissonJ * poissonJ) / youngJ;

          // Calculate the equivalent Young's modulus
          const SoAFloatPrecision youngEq = 1.0 / (effectiveYoungI + effectiveYoungJ);

          // Calculate the equivalent radius and the penetration depth
          const SoAFloatPrecision radEq = radI * radJ / (radI + radJ);

          const SoAFloatPrecision penDepth = rad - dr;

          // Normal force factor
          const SoAFloatPrecision normalFactorForce = mask ? 4. / 3. * youngEq * sqrt(radEq * penDepth) : 0;

          // Normal damping factor
          const SoAFloatPrecision normalFactorDamping = _dampingRatio * 2.0 * sqrt(massEq * normalFactorForce);

          // Calculate the normal force and the normal force vector
          const SoAFloatPrecision normalForce = mask ? normalFactorForce * penDepth / dr : 0.0;
          const SoAFloatPrecision normalDamping = mask ? normalFactorDamping : 0.0;

          const SoAFloatPrecision fx = (drx * normalForce - dvx * normalDamping);
          const SoAFloatPrecision fy = (dry * normalForce - dvy * normalDamping);
          const SoAFloatPrecision fz = (drz * normalForce - dvz * normalDamping);

          fxacc += fx;
          fyacc += fy;
          fzacc += fz;

          if (newton3) {
            fxArr[j] = fx;
            fyArr[j] = fy;
            fzArr[j] = fz;
          }
        }
      }
    }

    // this loop goes over the remainder and uses no optimizations
    for (size_t jNeighIndex = joff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;

      // Skipp if dummy
      const auto ownedStateJ = ownedStatePtr[j];
      if (ownedStateJ == autopas::OwnershipState::dummy) {
        continue;
      }

      // Calculate pair sizes, positions and distances
      const SoAFloatPrecision drx = xptr[indexFirst] - xptr[j];
      const SoAFloatPrecision dry = yptr[indexFirst] - yptr[j];
      const SoAFloatPrecision drz = zptr[indexFirst] - zptr[j];

      const SoAFloatPrecision drx2 = drx * drx;
      const SoAFloatPrecision dry2 = dry * dry;
      const SoAFloatPrecision drz2 = drz * drz;

      const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
      const SoAFloatPrecision dr = sqrt(dr2);

      const SoAFloatPrecision radI = radptr[indexFirst];
      const SoAFloatPrecision radJ = radptr[j];

      const SoAFloatPrecision rad = radI + radJ;

      // Skip if particles do not intersect or have the same position
      if (dr >= rad || dr == 0) {
        continue;
      }

      // Relative particle velocity
      const SoAFloatPrecision dvx = vxptr[indexFirst] - vxptr[j];
      const SoAFloatPrecision dvy = vyptr[indexFirst] - vyptr[j];
      const SoAFloatPrecision dvz = vzptr[indexFirst] - vzptr[j];

      // Particle masses and equivalent mass
      const SoAFloatPrecision massI = massptr[indexFirst];
      const SoAFloatPrecision massJ = massptr[j];

      const SoAFloatPrecision massEq = massI * massJ / (massI + massJ);

      // Young's modulus and Poisson ratio of partice I and J
      const SoAFloatPrecision poissonI = poissonptr[indexFirst];
      const SoAFloatPrecision youngI = youngptr[indexFirst];
      const SoAFloatPrecision poissonJ = poissonptr[j];
      const SoAFloatPrecision youngJ = youngptr[j];

      const SoAFloatPrecision effectiveYoungI = (1.0 - poissonI * poissonI) / youngI;
      const SoAFloatPrecision effectiveYoungJ = (1.0 - poissonJ * poissonJ) / youngJ;

      // Calculate the equivalent Young's modulus
      const SoAFloatPrecision youngEq = 1.0 / (effectiveYoungI + effectiveYoungJ);

      // Calculate the equivalent radius and the penetration depth
      const SoAFloatPrecision radEq = radI * radJ / (radI + radJ);

      const SoAFloatPrecision penDepth = rad - dr;

      // Normal force factor
      const SoAFloatPrecision normalFactorForce = 4. / 3. * youngEq * sqrt(radEq * penDepth);

      // Normal damping factor
      const SoAFloatPrecision normalFactorDamping = _dampingRatio * 2.0 * sqrt(massEq * normalFactorForce);

      // Calculate the normal force and the normal force vector
      const SoAFloatPrecision normalForce = normalFactorForce * penDepth / dr;
      const SoAFloatPrecision normalDamping = normalFactorDamping;

      const SoAFloatPrecision fx = (drx * normalForce - dvx * normalDamping);
      const SoAFloatPrecision fy = (dry * normalForce - dvy * normalDamping);
      const SoAFloatPrecision fz = (drz * normalForce - dvz * normalDamping);

      fxacc += fx;
      fyacc += fy;
      fzacc += fz;

      if (newton3) {
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;
      }
    }

    if (fxacc != 0 or fyacc != 0 or fzacc != 0) {
      fxptr[indexFirst] += fxacc;
      fyptr[indexFirst] += fyacc;
      fzptr[indexFirst] += fzacc;
    }
  }

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // Damping ratio
  double _dampingRatio = 0.1;
};

// double DEMFunctor<Particle>::_factorSubtractExcessForces{1.0};

};  // namespace demLib

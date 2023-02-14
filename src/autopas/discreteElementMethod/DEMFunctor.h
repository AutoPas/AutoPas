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

namespace autopas {

/**
 * @brief Functor for the force interactions between two DEM objects (particles)
 * 
 * @tparam Particle 
 * @tparam useNewton3 
 * @tparam relevantForTuning 
 */
template <class Particle, FunctorN3Modes useNewton3 = FunctorN3Modes::Both, bool relevantForTuning = true>
class DEMFunctor
    : public Functor<Particle,
                      DEMFunctor<Particle, useNewton3, relevantForTuning>> {

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

    : Functor<Particle, DEMFunctor<Particle, useNewton3, relevantForTuning>>(cutoff) {}

  bool isRelevantForTuning() final { 
    return relevantForTuning; 
  }

  bool allowsNewton3() final { 
    return useNewton3 == FunctorN3Modes::Newton3Only or useNewton3 == FunctorN3Modes::Both; 
  }

  bool allowsNonNewton3() final {
    return useNewton3 == FunctorN3Modes::Newton3Off or useNewton3 == FunctorN3Modes::Both;
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
   * with delta the penetration depth, R1 and R2 the particle radii and
   * E1 and E2 the Young's moduli and v1 and v2 the Poisson ratios of the colliding particles.
   * 
   * @param i Particle i 
   * @param j Particle j
   * @param newton3 Use Newton's third law
   */
  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    if(i.isDummy() or j.isDummy()) {
      return;
    }

    // Get particle radii
    double radi = i.getRad();
    double radj = j.getRad();

    //distance between ParticleCenters
    auto dr = utils::ArrayMath::sub(i.getR(), j.getR()); 
    double penDepth = radi+radj-utils::ArrayMath::L2Norm(dr);

    double penDepth3 = penDepth * penDepth * penDepth;

    // Particles do not intersect => return
    if (penDepth <= 0) {
      return;
    }

    // Calculate equivalent radius
    double radEq = radi * radj / ( radi + radj );

    // get the Poisson ratios
    double poissoni = i.getPoisson();
    double poissonj = j.getPoisson();

    // get the Young's moduli
    double youngi = i.getYoung();
    double youngj = j.getYoung();

    double ei = ( 1.0 - poissoni * poissoni ) / youngi;
    double ej = ( 1.0 - poissonj * poissonj ) / youngj;

    // Calculate the equivalent Young's modulus
    double youngEq = 1.0 / ( ei + ej );

    // Calculate the normal force and the normal force vector
    double normalForce = 4./3. * youngEq * sqrt(radEq * penDepth3); 

    auto vecNormalForce = utils::ArrayMath::mulScalar (utils::ArrayMath::normalize(dr),normalForce);

    i.addF(vecNormalForce);

    if (newton3) {
      j.subF(vecNormalForce);
    }

  }

  /**
   * 
   * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
   * 
   * DEM Functor for an SoA
   * 
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) final {

    // If SoA is empty, return
    if (soa.getNumberOfParticles() == 0) return;

    // Pointers to the start of arrays containing relevant particle data
    const auto *const __restrict idptr = soa.template begin<Particle::AttributeNames::id>();

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict radptr = soa.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poissonptr = soa.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict youngptr = soa.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownershipPtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    // Go over all particles in the SoA
    for (unsigned int i = 0; i < soa.getNumberOfParticles(); i++) {

      // If particle is dummy, skip it
      const auto ownedStateI = ownershipPtr[i];
      if (ownedStateI == OwnershipState::dummy) {
        continue;
      }

      //accumulating Forces
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = i + 1; j < soa.getNumberOfParticles(); j++) {

        // Calculate pair sizes, positions and distances
        const auto ownedStateJ = ownershipPtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
        const SoAFloatPrecision dr  = sqrt( dr2 );

        const SoAFloatPrecision radI = radptr[i];
        const SoAFloatPrecision radJ = radptr[j];
        const SoAFloatPrecision rad = radI + radJ;

        // Mask away if particles arent intersecting or if j is dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr <= rad and ownedStateJ != OwnershipState::dummy;
        
        // Young's modulus and Poisson ratio of partice I and J
        const SoAFloatPrecision poissonI = poissonptr[i];
        const SoAFloatPrecision youngI = youngptr[i];
        const SoAFloatPrecision poissonJ = poissonptr[j];
        const SoAFloatPrecision youngJ = youngptr[j];

        const SoAFloatPrecision eI = ( 1.0 - poissonI * poissonI ) / youngI;
        const SoAFloatPrecision eJ = ( 1.0 - poissonJ * poissonJ ) / youngJ;

        // Calculate the equivalent Young's modulus
        const SoAFloatPrecision youngEq = 1.0 / ( eI + eJ );

        // Calculate the equivalent radius and the penetration depth
        const SoAFloatPrecision radEq = radI * radJ / ( radI + radJ );

        const SoAFloatPrecision penDepth = rad - dr;
        const SoAFloatPrecision penDepth3 = penDepth * penDepth * penDepth;

        // Calculate the normal force and the normal force vector
        const SoAFloatPrecision normalForce = mask ? 4./3. * youngEq * sqrt(radEq * penDepth3) / dr : 0.0;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

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
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, const bool newton3) final {
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
  void SoAFunctorPairImpl(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2) {

    // If SoA is empty, return
    if (soa1.getNumberOfParticles() == 0 || soa2.getNumberOfParticles() == 0) return;

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

    const auto *const __restrict rad1ptr = soa1.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poisson1ptr = soa1.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict young1ptr = soa1.template begin<Particle::AttributeNames::young>();
    const auto *const __restrict rad2ptr = soa2.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poisson2ptr = soa2.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict young2ptr = soa2.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownership1Ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownership2Ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();

    // Go over all particles in the SoA
    for (unsigned int i = 0; i < soa1.getNumberOfParticles(); i++) {
      
      // If particle is dummy, skip it
      const auto ownedStateI = ownership1Ptr[i];
      if (ownedStateI == OwnershipState::dummy) {
        continue;
      }

      //accumulating Force directions
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

#pragme omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = 0; j < soa2.getNumberOfParticles(); ++j) {

        const auto ownedStateJ = ownership2Ptr[j];

        // Calculate pair sizes, positions and distances
        const SoAFloatPrecision drx = x1ptr[i] - x2ptr[j];
        const SoAFloatPrecision dry = y1ptr[i] - y2ptr[j];
        const SoAFloatPrecision drz = z1ptr[i] - z2ptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
        const SoAFloatPrecision dr  = sqrt( dr2 );

        const SoAFloatPrecision radI = rad1ptr[i];
        const SoAFloatPrecision radJ = rad2ptr[j];
        const SoAFloatPrecision rad = radI + radJ;

        // Mask away if particles arent intersecting or if j is dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr <= rad and ownedStateJ != OwnershipState::dummy;

        // Young's modulus and Poisson ratio of partice I and J
        const SoAFloatPrecision poissonI = poisson1ptr[i];
        const SoAFloatPrecision youngI = young1ptr[i];
        const SoAFloatPrecision poissonJ = poisson1ptr[j];
        const SoAFloatPrecision youngJ = young1ptr[j];

        const SoAFloatPrecision eI = ( 1.0 - poissonI * poissonI ) / youngI;
        const SoAFloatPrecision eJ = ( 1.0 - poissonJ * poissonJ ) / youngJ;

        // Calculate the equivalent Young's modulus
        const SoAFloatPrecision youngEq = 1.0 / ( eI + eJ );

        // Calculate the equivalent radius and the penetration depth
        const SoAFloatPrecision radEq = radI * radJ / ( radI + radJ );

        const SoAFloatPrecision penDepth = rad - dr;
        const SoAFloatPrecision penDepth3 = penDepth * penDepth * penDepth;

        // Calculate the normal force and the normal force vector
        const SoAFloatPrecision normalForce = mask ? 4./3. * youngEq * sqrt(radEq * penDepth3) / dr : 0.0;

        const SoAFloatPrecision fx = drx * normalForce;
        const SoAFloatPrecision fy = dry * normalForce;
        const SoAFloatPrecision fz = drz * normalForce;

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
   * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   * 
   * DEM Functor for a Verlet SoA
   * 
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    
    // Return, if SoA or neighbor list empty
    if (soa.getNumberOfParticles() == 0 or neighborList.empty()) return;
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
    return std::array<typename Particle::AttributeNames, 12>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::rad,    Particle::AttributeNames::poisson, 
        Particle::AttributeNames::young,  Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 9>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::rad,    Particle::AttributeNames::poisson, 
        Particle::AttributeNames::young,  Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

  /**
   * @brief THIS has to be recalculated: Get the Num Flops Per Kernel Call object
   * 
   * @return unsigned long 
   */
  static unsigned long getNumFlopsPerKernelCall() {
    // Kernel: THIS has to be recalculated
    return 18ul;
  }

  void initTraversal() final {
    _postProcessed = false;
  }

  void endTraversal(bool newton3) final {
    if (_postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
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
  void SoAFunctorVerletImpl(SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {

    // Pointers to particle properties
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

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
    if (ownedStateI == OwnershipState::dummy) {
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
      alignas(64) std::array<SoAFloatPrecision, vecsize> xtmp, ytmp, ztmp, radtmp, poissontmp, youngtmp, xArr, yArr, zArr, radArr, poissonArr, youngArr, fxArr, fyArr, fzArr;
      alignas(64) std::array<OwnershipState, vecsize> ownedStateArr{};
      
      // broadcast of the position of particle i
      for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
        xtmp[tmpj] = xptr[indexFirst];
        ytmp[tmpj] = yptr[indexFirst];
        ztmp[tmpj] = zptr[indexFirst];
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
          const bool mask = dr <= rad and ownedStateJ != OwnershipState::dummy and dr != 0;

          // Young's modulus and Poisson ratio of partice I and J
          const SoAFloatPrecision poissonI = poissontmp[j];
          const SoAFloatPrecision youngI = youngtmp[j];
          const SoAFloatPrecision poissonJ = poissonArr[j];
          const SoAFloatPrecision youngJ = youngArr[j];

          const SoAFloatPrecision eI = ( 1.0 - poissonI * poissonI ) / youngI;
          const SoAFloatPrecision eJ = ( 1.0 - poissonJ * poissonJ ) / youngJ;

          // Calculate the equivalent Young's modulus
          const SoAFloatPrecision youngEq = 1.0 / ( eI + eJ );

          // Calculate the equivalent radius and the penetration depth
          const SoAFloatPrecision radEq = radI * radJ / ( radI + radJ );

          const SoAFloatPrecision penDepth = rad - dr;
          const SoAFloatPrecision penDepth3 = penDepth * penDepth * penDepth;

          // Calculate the normal force and the normal force vector
          const SoAFloatPrecision normalForce = mask ? 4./3. * youngEq * sqrt(radEq * penDepth3) / dr : 0.0;

          const SoAFloatPrecision fx = drx * normalForce;
          const SoAFloatPrecision fy = dry * normalForce;
          const SoAFloatPrecision fz = drz * normalForce;

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
      if (ownedStateJ == OwnershipState::dummy) {
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
      
      // Young's modulus and Poisson ratio of partice I and J
      const SoAFloatPrecision poissonI = poissonptr[indexFirst];
      const SoAFloatPrecision youngI = youngptr[indexFirst];
      const SoAFloatPrecision poissonJ = poissonptr[j];
      const SoAFloatPrecision youngJ = youngptr[j];

      const SoAFloatPrecision eI = ( 1.0 - poissonI * poissonI ) / youngI;
      const SoAFloatPrecision eJ = ( 1.0 - poissonJ * poissonJ ) / youngJ;

      // Calculate the equivalent Young's modulus
      const SoAFloatPrecision youngEq = 1.0 / ( eI + eJ );

      // Calculate the equivalent radius and the penetration depth
      const SoAFloatPrecision radEq = radI * radJ / ( radI + radJ );

      const SoAFloatPrecision penDepth = rad - dr;
      const SoAFloatPrecision penDepth3 = penDepth * penDepth * penDepth;

      // Calculate the normal force and the normal force vector
      const SoAFloatPrecision normalForce = 4./3. * youngEq * sqrt(radEq * penDepth3) / dr;

      const SoAFloatPrecision fx = drx * normalForce;
      const SoAFloatPrecision fy = dry * normalForce;
      const SoAFloatPrecision fz = drz * normalForce;

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

};

}; //namespace autopas
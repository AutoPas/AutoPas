/**
 * @file AxilrodTellerFunctorAVX512.h
 * @author M. Muehlhaeusser
 * @date 25/07/23
 */

#pragma once

#include <immintrin.h>

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
 * @tparam countFLOPs counts FLOPs and hitrate
 */
template <class Particle, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGlobals = false, bool countFLOPs = false>
class AxilrodTellerFunctor
    : public autopas::TriwiseFunctor<
          Particle, AxilrodTellerFunctorAVX512<Particle, useMixing, useNewton3, calculateGlobals, countFLOPS>> {
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
  AxilrodTellerFunctorAVX512() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit AxilrodTellerFunctorAVX512(double cutoff, void * /*dummy*/)
      : autopas::TriwiseFunctor<Particle,
                                AxilrodTellerFunctor<Particle, useMixing, useNewton3, calculateGlobals, countFLOPs>>(
            cutoff),
        _cutoffSquaredAoS{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadDataGlobals(),
        _postProcessed{false},
        _cutoffSquared{_mm512_set1_pd(cutoff * cutoff)} {
    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
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
  explicit AxilrodTellerFunctorAVX512(double cutoff) : AxilrodTellerFunctorAVX512(cutoff, nullptr) {
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
  explicit AxilrodTellerFunctorAVX512(double cutoff,
                                      ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : AxilrodTellerFunctorAVX512(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "AxilrodTellerFunctorAVX512"; }

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

    const auto threadnum = autopas::autopas_get_thread_num();

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
    }

    auto nu = _nuAoS;
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
    if (distSquaredIJ > _cutoffSquaredAoS or distSquaredJK > _cutoffSquaredAoS or distSquaredKI > _cutoffSquaredAoS) {
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
      const auto virialI = forceI * i.getR();
      if (i.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy3;
        _aosThreadDataGlobals[threadnum].virialSum += virialI;
      }
      // for non-newton3 particles j and/or k will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        const auto virialJ = forceJ * j.getR();
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy3;
        _aosThreadDataGlobals[threadnum].virialSum += virialJ;
      }
      if (newton3 and k.isOwned()) {
        const auto virialK = forceK * k.getR();
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

  inline bool SoAParticlesInCutoff(const double *const x1ptr, const double *const y1ptr, const double *const z1ptr,
                                   const double *const x2ptr, const double *const y2ptr, const double *const z2ptr,
                                   size_t index1, size_t index2) {
    const SoAFloatPrecision drxij = x2ptr[index2] - x1ptr[index1];
    const SoAFloatPrecision dryij = y2ptr[index2] - y1ptr[index1];
    const SoAFloatPrecision drzij = z2ptr[index2] - z1ptr[index1];

    const SoAFloatPrecision drxij2 = drxij * drxij;
    const SoAFloatPrecision dryij2 = dryij * dryij;
    const SoAFloatPrecision drzij2 = drzij * drzij;

    const SoAFloatPrecision drij2 = drxij2 + dryij2 + drzij2;
    return drij2 <= _cutoffSquaredAoS;
  }

  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImpl<true>(soa);
    } else {
      SoAFunctorSingleImpl<false>(soa);
    }
  }

  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

  /**
   * Functor for structure of arrays (SoA)
   *
   * This functor calculates the forces
   * between all particles of soa1 and soa2 and soa3.
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
   * @param soa3 Third structure of arrays.
   * @param newton3 defines whether or whether not to use newton 3
   */
  void SoAFunctorTriple(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                        autopas::SoAView<SoAArraysType> soa3, const bool newton3) {
    if (newton3) {
      SoAFunctorTripleImpl<true>(soa1, soa2, soa3);
    } else {
      SoAFunctorTripleImpl<false>(soa1, soa2, soa3);
    }
  }

 private:
  /**
   * SoA kernel using masked simd intrinsics.
   */
  template <bool newton3j, bool newton3k, bool remainderIsMasked>
  void SoAKernelMasked(const size_t k, const __m512d &xi, const __m512d &yi, const __m512d &zi, const __m512d &xj,
                       const __m512d &yj, const __m512d &zj, const double *const __restrict xkPtr,
                       const double *const __restrict ykPtr, const double *const __restrict zkPtr, const size_t typeI,
                       const size_t typeJ, const size_t *const __restrict typeKptr, const __mmask8 &ownedMaskI,
                       const __mmask8 &ownedMaskJ, const autopas::OwnershipState *const ownedStateKptr,
                       const __m512d &drxij, const __m512d &dryij, const __m512d &drzij, const __m512d &drij2,
                       __m512d &fxiacc, __m512d &fyiacc, __m512d &fziacc, __m512d &fxjacc, __m512d &fyjacc,
                       __m512d &fzjacc, double *const __restrict fxkPtr, double *const __restrict fykPtr,
                       double *const __restrict fzkPtr, __m512d &potentialEnergySum, __m512d &virialSumX,
                       __m512d &virialSumY, __m512d &virialSumZ, unsigned int rest = 0) {
    const __m512d xk = remainderIsMasked ? _mm512_maskz_loadu_pd(_masks[rest], &xkPtr[k]) : _mm512_loadu_pd(&xkPtr[k]);
    const __m512d yk = remainderIsMasked ? _mm512_maskz_loadu_pd(_masks[rest], &ykPtr[k]) : _mm512_loadu_pd(&ykPtr[k]);
    const __m512d zk = remainderIsMasked ? _mm512_maskz_loadu_pd(_masks[rest], &zkPtr[k]) : _mm512_loadu_pd(&zkPtr[k]);

    // only required for calculating globals
    const __m512i ownedStateK = remainderIsMasked ? _mm512_maskz_loadu_epi64(_masks[rest], &ownedStateKptr[k])
                                                  : _mm512_loadu_epi64(&ownedStateKptr[k]);
    const __mmask8 ownedMaskK = _mm512_cmp_epi64_mask(ownedStateK, _ownedStateOwnedMM512i, _MM_CMPINT_EQ);

    // calculate distance k-i
    const __m512d drxki = _mm512_sub_pd(xi, xk);
    const __m512d dryki = _mm512_sub_pd(yi, yk);
    const __m512d drzki = _mm512_sub_pd(zi, zk);

    const __m512d drxki2 = _mm512_mul_pd(drxki, drxki);
    const __m512d drki2PART = _mm512_fmadd_pd(dryki, dryki, drxki2);
    const __m512d drki2 = _mm512_fmadd_pd(drzki, drzki, drki2PART);

    const __mmask8 cutoffMask_ki = _mm512_cmp_pd_mask(drki2, _cutoffSquared, _CMP_LE_OS);

    // calculate distance j-k
    const __m512d drxjk = _mm512_sub_pd(xk, xj);
    const __m512d dryjk = _mm512_sub_pd(yk, yj);
    const __m512d drzjk = _mm512_sub_pd(zk, zj);

    const __m512d drxjk2 = _mm512_mul_pd(drxjk, drxjk);
    const __m512d drjk2PART = _mm512_fmadd_pd(dryjk, dryjk, drxjk2);
    const __m512d drjk2 = _mm512_fmadd_pd(drzjk, drzjk, drjk2PART);

    const __mmask8 cutoffMask_jk = _mm512_cmp_pd_mask(drjk2, _cutoffSquared, _CMP_LE_OS);

    const __mmask8 cutoffMask = _mm512_kand(cutoffMask_jk, cutoffMask_ki);

    const __mmask8 dummyMask = _mm512_cmp_epi64_mask(ownedStateK, _ownedStateDummyMM512i, _MM_CMPINT_NE);

    // mask with bits set if particle is not a dummy and within cutoff of the first two particles
    const __mmask8 mask = remainderIsMasked ? _mm512_kand(_masks[rest], _mm512_kand(cutoffMask, dummyMask))
                                            : _mm512_kand(cutoffMask, dummyMask);

    if (std::bitset<vecLength>(mask).count() == 0) {
      // early exit if all particles are ignored
      return;
    }

    // calculate force
    const __m512d drxi2 = _mm512_mul_pd(drxij, drxki);
    const __m512d dri2PART = _mm512_fmadd_pd(dryij, dryki, drxi2);
    const __m512d dri2 = _mm512_fmadd_pd(drzij, drzki, dri2PART);

    const __m512d drxj2 = _mm512_mul_pd(drxij, drxjk);
    const __m512d drj2PART = _mm512_fmadd_pd(dryij, dryjk, drxj2);
    const __m512d drj2 = _mm512_fmadd_pd(drzij, drzjk, drj2PART);

    const __m512d drxk2 = _mm512_mul_pd(drxjk, drxki);
    const __m512d drk2PART = _mm512_fmadd_pd(dryjk, dryki, drxk2);
    const __m512d drk2 = _mm512_fmadd_pd(drzjk, drzki, drk2PART);

    const __m512d drijk2PART = _mm512_mul_pd(dri2, drj2);
    const __m512d drijk2 = _mm512_mul_pd(drijk2PART, drk2);

    const __m512d dr2PART = _mm512_mul_pd(drij2, drjk2);
    const __m512d dr2 = _mm512_mul_pd(dr2PART, drki2);
    const __m512d dr2sqrt = _mm512_sqrt_pd(dr2);
    const __m512d dr5PART = _mm512_mul_pd(dr2, dr2);
    const __m512d dr5 = _mm512_mul_pd(dr5PART, dr2sqrt);

    if constexpr (useMixing) {
      _nu =
          _mm512_set_pd(not remainderIsMasked or rest > 7 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[k + 7]) : 0,
                        not remainderIsMasked or rest > 6 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[k + 6]) : 0,
                        not remainderIsMasked or rest > 5 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[k + 5]) : 0,
                        not remainderIsMasked or rest > 4 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[k + 4]) : 0,
                        not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[k + 3]) : 0,
                        not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[k + 2]) : 0,
                        not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[k + 1]) : 0,
                        _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[k]));
    }
    const __m512d invdr5 = _mm512_maskz_div_pd(mask, _nu, dr5);
    const __m512d invdr53 = _mm512_mul_pd(_three, invdr5);

    const __m512d drj2_drk2 = _mm512_sub_pd(drj2, drk2);
    const __m512d i_jkFactor = _mm512_mul_pd(dri2, drj2_drk2);
    const __m512d drj2drk2 = _mm512_mul_pd(drj2, drk2);
    const __m512d drjk2drki2 = _mm512_mul_pd(drjk2, drki2);
    const __m512d drijk2BYdrij2 = _mm512_maskz_div_pd(mask, drijk2, drij2);
    const __m512d i_ijFactorPART = _mm512_sub_pd(drj2drk2, drjk2drki2);
    const __m512d i_ijFactor = _mm512_fmadd_pd(drijk2BYdrij2, _five, i_ijFactorPART);
    const __m512d drij2drjk2 = _mm512_mul_pd(drij2, drjk2);
    const __m512d drijk2BYdrki2 = _mm512_maskz_div_pd(mask, drijk2, drki2);
    const __m512d i_kiFactorPART = _mm512_sub_pd(drj2drk2, drij2drjk2);
    const __m512d i_kiFactor = _mm512_fmadd_pd(drijk2BYdrki2, _five, i_kiFactorPART);

    const __m512d fxiPART = _mm512_mul_pd(drxjk, i_jkFactor);
    const __m512d fxiPART2 = _mm512_fmadd_pd(drxij, i_ijFactor, fxiPART);
    const __m512d fxiPART3 = _mm512_fnmadd_pd(drxki, i_kiFactor, fxiPART2);
    const __m512d fxi = _mm512_mul_pd(invdr53, fxiPART3);
    fxiacc = _mm512_mask_add_pd(fxiacc, mask, fxi, fxiacc);
    const __m512d fyiPART = _mm512_mul_pd(dryjk, i_jkFactor);
    const __m512d fyiPART2 = _mm512_fmadd_pd(dryij, i_ijFactor, fyiPART);
    const __m512d fyiPART3 = _mm512_fnmadd_pd(dryki, i_kiFactor, fyiPART2);
    const __m512d fyi = _mm512_mul_pd(invdr53, fyiPART3);
    fyiacc = _mm512_mask_add_pd(fyiacc, mask, fyi, fyiacc);
    const __m512d fziPART = _mm512_mul_pd(drzjk, i_jkFactor);
    const __m512d fziPART2 = _mm512_fmadd_pd(drzij, i_ijFactor, fziPART);
    const __m512d fziPART3 = _mm512_fnmadd_pd(drzki, i_kiFactor, fziPART2);
    const __m512d fzi = _mm512_mul_pd(invdr53, fziPART3);
    fziacc = _mm512_mask_add_pd(fziacc, mask, fzi, fziacc);

    if constexpr (calculateGlobals) {
      const __m512d virialXI = _mm512_mul_pd(fxi, xi);
      virialSumX = _mm512_mask_add_pd(virialSumX, ownedMaskI, virialSumX, virialXI);
      const __m512d virialYI = _mm512_mul_pd(fyi, yi);
      virialSumY = _mm512_mask_add_pd(virialSumY, ownedMaskI, virialSumY, virialYI);
      const __m512d virialZI = _mm512_mul_pd(fzi, zi);
      virialSumZ = _mm512_mask_add_pd(virialSumZ, ownedMaskI, virialSumZ, virialZI);
    }

    if constexpr (newton3j) {
      const __m512d drk2_dri2 = _mm512_sub_pd(drk2, dri2);
      const __m512d j_kiFactor = _mm512_mul_pd(drj2, drk2_dri2);
      const __m512d dri2drk2 = _mm512_mul_pd(dri2, drk2);
      const __m512d j_ijFactorPART = _mm512_sub_pd(dri2drk2, drjk2drki2);
      const __m512d j_ijFactor = _mm512_fmadd_pd(drijk2BYdrij2, _five, j_ijFactorPART);
      const __m512d drij2drki2 = _mm512_mul_pd(drij2, drki2);
      const __m512d drijk2BYdrjk2 = _mm512_maskz_div_pd(mask, drijk2, drjk2);
      const __m512d j_jkFactorPART = _mm512_sub_pd(dri2drk2, drij2drki2);
      const __m512d j_jkFactor = _mm512_fmadd_pd(drijk2BYdrjk2, _five, j_jkFactorPART);

      const __m512d fxjPART = _mm512_mul_pd(drxki, j_kiFactor);
      const __m512d fxjPART2 = _mm512_fnmadd_pd(drxij, j_ijFactor, fxjPART);
      const __m512d fxjPART3 = _mm512_fmadd_pd(drxjk, j_jkFactor, fxjPART2);
      const __m512d fxj = _mm512_mul_pd(invdr53, fxjPART3);
      fxjacc = _mm512_mask_add_pd(fxjacc, mask, fxj, fxjacc);
      const __m512d fyjPART = _mm512_mul_pd(dryki, j_kiFactor);
      const __m512d fyjPART2 = _mm512_fnmadd_pd(dryij, j_ijFactor, fyjPART);
      const __m512d fyjPART3 = _mm512_fmadd_pd(dryjk, j_jkFactor, fyjPART2);
      const __m512d fyj = _mm512_mul_pd(invdr53, fyjPART3);
      fyjacc = _mm512_mask_add_pd(fyjacc, mask, fyj, fyjacc);
      const __m512d fzjPART = _mm512_mul_pd(drzki, j_kiFactor);
      const __m512d fzjPART2 = _mm512_fnmadd_pd(drzij, j_ijFactor, fzjPART);
      const __m512d fzjPART3 = _mm512_fmadd_pd(drzjk, j_jkFactor, fzjPART2);
      const __m512d fzj = _mm512_mul_pd(invdr53, fzjPART3);
      fzjacc = _mm512_mask_add_pd(fzjacc, mask, fzj, fzjacc);

      if constexpr (calculateGlobals) {
        const __m512d virialXJ = _mm512_mul_pd(fxj, xj);
        virialSumX = _mm512_mask_add_pd(virialSumX, ownedMaskJ, virialSumX, virialXJ);
        const __m512d virialYJ = _mm512_mul_pd(fyj, yj);
        virialSumY = _mm512_mask_add_pd(virialSumY, ownedMaskJ, virialSumY, virialYJ);
        const __m512d virialZJ = _mm512_mul_pd(fzj, zj);
        virialSumZ = _mm512_mask_add_pd(virialSumZ, ownedMaskJ, virialSumZ, virialZJ);
      }

      if constexpr (newton3k) {
        const __m512d nfxk = _mm512_add_pd(fxi, fxj);
        const __m512d nfyk = _mm512_add_pd(fyi, fyj);
        const __m512d nfzk = _mm512_add_pd(fzi, fzj);

        const __m512d fxk_old = _mm512_maskz_loadu_pd(mask, &fxkPtr[k]);
        const __m512d fyk_old = _mm512_maskz_loadu_pd(mask, &fykPtr[k]);
        const __m512d fzk_old = _mm512_maskz_loadu_pd(mask, &fzkPtr[k]);

        const __m512d fxk_new = _mm512_sub_pd(fxk_old, nfxk);
        const __m512d fyk_new = _mm512_sub_pd(fyk_old, nfyk);
        const __m512d fzk_new = _mm512_sub_pd(fzk_old, nfzk);

        _mm512_mask_storeu_pd(&fxkPtr[k], mask, fxk_new);
        _mm512_mask_storeu_pd(&fykPtr[k], mask, fyk_new);
        _mm512_mask_storeu_pd(&fzkPtr[k], mask, fzk_new);

        if constexpr (calculateGlobals) {
          const __m512d virialXK = _mm512_mul_pd(nfxk, xk);
          virialSumX = _mm512_mask_sub_pd(virialSumX, ownedMaskK, virialSumX, virialXK);
          const __m512d virialYK = _mm512_mul_pd(nfyk, yk);
          virialSumY = _mm512_mask_sub_pd(virialSumY, ownedMaskK, virialSumY, virialYK);
          const __m512d virialZK = _mm512_mul_pd(nfzk, zk);
          virialSumZ = _mm512_mask_sub_pd(virialSumZ, ownedMaskK, virialSumZ, virialZK);
        }
      }
    }

    if constexpr (calculateGlobals) {
      // Add 3 * potential energy to every owned particle of the interaction.
      // Division to the correct value is handled in endTraversal().
      // invdr53 is zero if particles are not within cutoff, so no cutoff mask is needed here.
      const __m512d potentialEnergy3 = _mm512_mul_pd(invdr53, _mm512_fnmadd_pd(dr2, _three, drijk2));
      potentialEnergySum = _mm512_mask_add_pd(potentialEnergySum, ownedMaskI, potentialEnergySum, potentialEnergy3);
      if constexpr (newton3j) {
        potentialEnergySum = _mm512_mask_add_pd(potentialEnergySum, ownedMaskJ, potentialEnergySum, potentialEnergy3);
      }
      if constexpr (newton3k) {
        potentialEnergySum = _mm512_mask_add_pd(potentialEnergySum, ownedMaskK, potentialEnergySum, potentialEnergy3);
      }
    }
  }

  /**
   * Implementation function of SoAFunctorSingle(soa, newton3)
   *
   * @tparam newton3
   * @param soa
   */
  template <bool newton3>
  inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
    if (soa.size() == 0) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    __m512d virialSumX = _mm512_setzero_pd();
    __m512d virialSumY = _mm512_setzero_pd();
    __m512d virialSumZ = _mm512_setzero_pd();
    __m512d potentialEnergySum = _mm512_setzero_pd();

    if constexpr (not useMixing) {
      _nu = _mm512_set1_pd(_nuAoS);
    }

    for (size_t i = soa.size() - 1; static_cast<long>(i) >= 2; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        // skip dummy particles
        continue;
      }

      // only required for calculating globals
      const __mmask8 ownedMaskI = ownedStatePtr[i] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

      const __m512d xi = _mm512_set1_pd(xptr[i]);
      const __m512d yi = _mm512_set1_pd(yptr[i]);
      const __m512d zi = _mm512_set1_pd(zptr[i]);

      __m512d fxiacc = _mm512_setzero_pd();
      __m512d fyiacc = _mm512_setzero_pd();
      __m512d fziacc = _mm512_setzero_pd();

      for (size_t j = i - 1; static_cast<long>(j) >= 1; --j) {
        if (ownedStatePtr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(xptr, yptr, zptr, xptr, yptr, zptr, i, j)) {
          // skip dummy particles
          continue;
        }

        // only required for calculating globals
        const __mmask8 ownedMaskJ = ownedStatePtr[j] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

        const __m512d xj = _mm512_set1_pd(xptr[j]);
        const __m512d yj = _mm512_set1_pd(yptr[j]);
        const __m512d zj = _mm512_set1_pd(zptr[j]);

        __m512d fxjacc = _mm512_setzero_pd();
        __m512d fyjacc = _mm512_setzero_pd();
        __m512d fzjacc = _mm512_setzero_pd();

        // calculate distance i-j
        const __m512d drxij = _mm512_sub_pd(xj, xi);
        const __m512d dryij = _mm512_sub_pd(yj, yi);
        const __m512d drzij = _mm512_sub_pd(zj, zi);

        const __m512d drxij2 = _mm512_mul_pd(drxij, drxij);
        const __m512d dryij2 = _mm512_mul_pd(dryij, dryij);
        const __m512d drzij2 = _mm512_mul_pd(drzij, drzij);

        const __m512d drij2PART = _mm512_add_pd(drxij2, dryij2);
        const __m512d drij2 = _mm512_add_pd(drij2PART, drzij2);

        size_t k = 0;
        // loop up to multiple of vecLength
        // only works if vecLength is a power of 2
        for (; k < (j & ~(vecLength - 1)); k += vecLength) {
          SoAKernelMasked<true, true, false>(k, xi, yi, zi, xj, yj, zj, xptr, yptr, zptr, typeptr[i], typeptr[j],
                                             typeptr, ownedMaskI, ownedMaskJ, ownedStatePtr, drxij, dryij, drzij, drij2,
                                             fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fxptr, fyptr, fzptr,
                                             potentialEnergySum, virialSumX, virialSumY, virialSumZ);
        }
        // process remainder
        if (k < j) {
          SoAKernelMasked<true, true, true>(k, xi, yi, zi, xj, yj, zj, xptr, yptr, zptr, typeptr[i], typeptr[j],
                                            typeptr, ownedMaskI, ownedMaskJ, ownedStatePtr, drxij, dryij, drzij, drij2,
                                            fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fxptr, fyptr, fzptr,
                                            potentialEnergySum, virialSumX, virialSumY, virialSumZ, j - k);
        }

        fxptr[j] += _mm512_reduce_add_pd(fxjacc);
        fyptr[j] += _mm512_reduce_add_pd(fyjacc);
        fzptr[j] += _mm512_reduce_add_pd(fzjacc);
      }
      fxptr[i] += _mm512_reduce_add_pd(fxiacc);
      fyptr[i] += _mm512_reduce_add_pd(fyiacc);
      fzptr[i] += _mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      _aosThreadDataGlobals[threadnum].potentialEnergySum += _mm512_reduce_add_pd(potentialEnergySum);
      _aosThreadDataGlobals[threadnum].virialSum[0] += _mm512_reduce_add_pd(virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += _mm512_reduce_add_pd(virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += _mm512_reduce_add_pd(virialSumZ);
    }
  }
  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    __m512d potentialEnergySum = _mm512_setzero_pd();
    __m512d virialSumX = _mm512_setzero_pd();
    __m512d virialSumY = _mm512_setzero_pd();
    __m512d virialSumZ = _mm512_setzero_pd();

    if constexpr (not useMixing) {
      _nu = _mm512_set1_pd(_nuAoS);
    }

    // particle 1 always from soa1
    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        // skip dummy particles
        continue;
      }

      // only required for calculating globals
      const __mmask8 ownedMaskI = ownedState1ptr[i] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

      const __m512d xi = _mm512_set1_pd(x1ptr[i]);
      const __m512d yi = _mm512_set1_pd(y1ptr[i]);
      const __m512d zi = _mm512_set1_pd(z1ptr[i]);

      __m512d fxiacc = _mm512_setzero_pd();
      __m512d fyiacc = _mm512_setzero_pd();
      __m512d fziacc = _mm512_setzero_pd();

      // particle 2 from soa1 and 3 from soa2
      for (size_t j = i + 1; j < soa1.size(); ++j) {
        if (ownedState1ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x1ptr, y1ptr, z1ptr, i, j)) {
          // skip dummy particles
          continue;
        }

        // only required for calculating globals
        const __mmask8 ownedMaskJ = ownedState2ptr[j] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

        const __m512d xj = _mm512_set1_pd(x1ptr[j]);
        const __m512d yj = _mm512_set1_pd(y1ptr[j]);
        const __m512d zj = _mm512_set1_pd(z1ptr[j]);

        __m512d fxjacc = _mm512_setzero_pd();
        __m512d fyjacc = _mm512_setzero_pd();
        __m512d fzjacc = _mm512_setzero_pd();

        // calculate distance i-j
        const __m512d drxij = _mm512_sub_pd(xj, xi);
        const __m512d dryij = _mm512_sub_pd(yj, yi);
        const __m512d drzij = _mm512_sub_pd(zj, zi);

        const __m512d drxij2 = _mm512_mul_pd(drxij, drxij);
        const __m512d dryij2 = _mm512_mul_pd(dryij, dryij);
        const __m512d drzij2 = _mm512_mul_pd(drzij, drzij);

        const __m512d drij2PART = _mm512_add_pd(drxij2, dryij2);
        const __m512d drij2 = _mm512_add_pd(drij2PART, drzij2);

        // particle 3 always from soa2
        size_t k = 0;
        // loop up to multiple of vecLength
        // only works if vecLength is a power of 2
        for (; k < (soa2.size() & ~(vecLength - 1)); k += vecLength) {
          SoAKernelMasked<true, newton3, false>(
              k, xi, yi, zi, xj, yj, zj, x2ptr, y2ptr, z2ptr, type1ptr[i], type1ptr[j], type2ptr, ownedMaskI,
              ownedMaskJ, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc,
              fx2ptr, fy2ptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ);
        }
        // process remainder
        if (k < soa2.size()) {
          SoAKernelMasked<true, newton3, true>(
              k, xi, yi, zi, xj, yj, zj, x2ptr, y2ptr, z2ptr, type1ptr[i], type1ptr[j], type2ptr, ownedMaskI,
              ownedMaskJ, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc,
              fx2ptr, fy2ptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ, soa2.size() - k);
        }
        fx1ptr[j] += _mm512_reduce_add_pd(fxjacc);
        fy1ptr[j] += _mm512_reduce_add_pd(fyjacc);
        fz1ptr[j] += _mm512_reduce_add_pd(fzjacc);
      }

      // both particles 2 and 3 from soa2
      for (size_t j = soa2.size() - 1; static_cast<long>(j) >= 0; --j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, j)) {
          // skip dummy particles
          continue;
        }

        // only required for calculating globals
        const __mmask8 ownedMaskJ = ownedState2ptr[j] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

        const __m512d xj = _mm512_set1_pd(x2ptr[j]);
        const __m512d yj = _mm512_set1_pd(y2ptr[j]);
        const __m512d zj = _mm512_set1_pd(z2ptr[j]);

        __m512d fxjacc = _mm512_setzero_pd();
        __m512d fyjacc = _mm512_setzero_pd();
        __m512d fzjacc = _mm512_setzero_pd();

        // calculate distance i-j
        const __m512d drxij = _mm512_sub_pd(xj, xi);
        const __m512d dryij = _mm512_sub_pd(yj, yi);
        const __m512d drzij = _mm512_sub_pd(zj, zi);

        const __m512d drxij2 = _mm512_mul_pd(drxij, drxij);
        const __m512d dryij2 = _mm512_mul_pd(dryij, dryij);
        const __m512d drzij2 = _mm512_mul_pd(drzij, drzij);

        const __m512d drij2PART = _mm512_add_pd(drxij2, dryij2);
        const __m512d drij2 = _mm512_add_pd(drij2PART, drzij2);

        // particle 3 always from soa 2
        size_t k = 0;
        // loop up to multiple of vecLength
        // only works if vecLength is a power of 2
        for (; k < (j & ~(vecLength - 1)); k += vecLength) {
          SoAKernelMasked<newton3, newton3, false>(
              k, xi, yi, zi, xj, yj, zj, x2ptr, y2ptr, z2ptr, type1ptr[i], type2ptr[j], type2ptr, ownedMaskI,
              ownedMaskJ, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc,
              fx2ptr, fy2ptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ);
        }
        // process remainder
        if (k < j) {
          SoAKernelMasked<newton3, newton3, true>(
              k, xi, yi, zi, xj, yj, zj, x2ptr, y2ptr, z2ptr, type1ptr[i], type2ptr[j], type2ptr, ownedMaskI,
              ownedMaskJ, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc,
              fx2ptr, fy2ptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ, j - k);
        }

        if constexpr (newton3) {
          fx2ptr[j] += _mm512_reduce_add_pd(fxjacc);
          fy2ptr[j] += _mm512_reduce_add_pd(fyjacc);
          fz2ptr[j] += _mm512_reduce_add_pd(fzjacc);
        }
      }

      fx1ptr[i] += _mm512_reduce_add_pd(fxiacc);
      fy1ptr[i] += _mm512_reduce_add_pd(fyiacc);
      fz1ptr[i] += _mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      _aosThreadDataGlobals[threadnum].potentialEnergySum += _mm512_reduce_add_pd(potentialEnergySum);
      _aosThreadDataGlobals[threadnum].virialSum[0] += _mm512_reduce_add_pd(virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += _mm512_reduce_add_pd(virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += _mm512_reduce_add_pd(virialSumZ);
    }
  }

  /**
   * Implementation function of SoAFunctorTriple(soa1, soa2, soa3, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   * @param soa3
   */
  template <bool newton3>
  inline void SoAFunctorTripleImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                   autopas::SoAView<SoAArraysType> soa3) {
    if (soa1.size() == 0 or soa2.size() == 0 or soa3.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x3ptr = soa3.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y3ptr = soa3.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z3ptr = soa3.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState3ptr = soa3.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx3ptr = soa3.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy3ptr = soa3.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz3ptr = soa3.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type3ptr = soa3.template begin<Particle::AttributeNames::typeId>();

    __m512d potentialEnergySum = _mm512_setzero_pd();
    __m512d virialSumX = _mm512_setzero_pd();
    __m512d virialSumY = _mm512_setzero_pd();
    __m512d virialSumZ = _mm512_setzero_pd();

    if constexpr (not useMixing) {
      _nu = _mm512_set1_pd(_nuAoS);
    }

    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        // skip dummy particles
        continue;
      }

      // only required for calculating globals
      const __mmask8 ownedMaskI = ownedState1ptr[i] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

      const __m512d xi = _mm512_set1_pd(x1ptr[i]);
      const __m512d yi = _mm512_set1_pd(y1ptr[i]);
      const __m512d zi = _mm512_set1_pd(z1ptr[i]);

      __m512d fxiacc = _mm512_setzero_pd();
      __m512d fyiacc = _mm512_setzero_pd();
      __m512d fziacc = _mm512_setzero_pd();

      for (size_t j = 0; j < soa2.size(); ++j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, j)) {
          // skip dummy particles
          continue;
        }

        // only required for calculating globals
        const __mmask8 ownedMaskJ = ownedState2ptr[j] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

        const __m512d xj = _mm512_set1_pd(x2ptr[j]);
        const __m512d yj = _mm512_set1_pd(y2ptr[j]);
        const __m512d zj = _mm512_set1_pd(z2ptr[j]);

        __m512d fxjacc = _mm512_setzero_pd();
        __m512d fyjacc = _mm512_setzero_pd();
        __m512d fzjacc = _mm512_setzero_pd();

        // calculate distance i-j
        const __m512d drxij = _mm512_sub_pd(xj, xi);
        const __m512d dryij = _mm512_sub_pd(yj, yi);
        const __m512d drzij = _mm512_sub_pd(zj, zi);

        const __m512d drxij2 = _mm512_mul_pd(drxij, drxij);
        const __m512d dryij2 = _mm512_mul_pd(dryij, dryij);
        const __m512d drzij2 = _mm512_mul_pd(drzij, drzij);

        const __m512d drij2PART = _mm512_add_pd(drxij2, dryij2);
        const __m512d drij2 = _mm512_add_pd(drij2PART, drzij2);

        size_t k = 0;
        // loop up to multiple of vecLength
        // only works if vecLength is a power of 2
        for (; k < (soa3.size() & ~(vecLength - 1)); k += vecLength) {
          SoAKernelMasked<newton3, newton3, false>(
              k, xi, yi, zi, xj, yj, zj, x3ptr, y3ptr, z3ptr, type1ptr[i], type2ptr[j], type3ptr, ownedMaskI,
              ownedMaskJ, ownedState3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc,
              fx3ptr, fy3ptr, fz3ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ);
        }
        // process remainder
        if (k < soa3.size()) {
          SoAKernelMasked<newton3, newton3, true>(
              k, xi, yi, zi, xj, yj, zj, x3ptr, y3ptr, z3ptr, type1ptr[i], type2ptr[j], type3ptr, ownedMaskI,
              ownedMaskJ, ownedState3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc,
              fx3ptr, fy3ptr, fz3ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ, soa3.size() - k);
        }

        if constexpr (newton3) {
          fx2ptr[j] += _mm512_reduce_add_pd(fxjacc);
          fy2ptr[j] += _mm512_reduce_add_pd(fyjacc);
          fz2ptr[j] += _mm512_reduce_add_pd(fzjacc);
        }
      }

      fx1ptr[i] += _mm512_reduce_add_pd(fxiacc);
      fy1ptr[i] += _mm512_reduce_add_pd(fyiacc);
      fz1ptr[i] += _mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      _aosThreadDataGlobals[threadnum].potentialEnergySum += _mm512_reduce_add_pd(potentialEnergySum);
      _aosThreadDataGlobals[threadnum].virialSum[0] += _mm512_reduce_add_pd(virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += _mm512_reduce_add_pd(virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += _mm512_reduce_add_pd(virialSumZ);
    }
  }

  // -------------------------------------------------------------
  // functions for a compress + alignr + gather + scatter approach
  // TODO currently not used. remove completely or keep as an option?
  // -------------------------------------------------------------

  /**
   * SoA kernel using simd gather intrinsics to load predetermined batch of particles.
   */
  template <bool newton3j, bool newton3k, bool remainderIsMasked>
  void SoAKernelGatherScatter(const __m512i indicesK, const __m512d &xi, const __m512d &yi, const __m512d &zi,
                              const __m512d &xj, const __m512d &yj, const __m512d &zj,
                              const double *const __restrict xkPtr, const double *const __restrict ykPtr,
                              const double *const __restrict zkPtr, const size_t typeI, const size_t typeJ,
                              const size_t *const __restrict typeKptr, const __mmask8 &ownedMaskI,
                              const __mmask8 &ownedMaskJ,
                              const autopas::OwnershipState *const __restrict ownedStateKptr, const __m512d &drxij,
                              const __m512d &dryij, const __m512d &drzij, const __m512d &drij2, __m512d &fxiacc,
                              __m512d &fyiacc, __m512d &fziacc, __m512d &fxjacc, __m512d &fyjacc, __m512d &fzjacc,
                              double *const __restrict fxkPtr, double *const __restrict fykPtr,
                              double *const __restrict fzkPtr, __m512d &potentialEnergySum, __m512d &virialSumX,
                              __m512d &virialSumY, __m512d &virialSumZ, unsigned int rest = 0) {
    const __m512d xk = remainderIsMasked ? _mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, xkPtr, 8)
                                         : _mm512_i64gather_pd(indicesK, xkPtr, 8);
    const __m512d yk = remainderIsMasked ? _mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, ykPtr, 8)
                                         : _mm512_i64gather_pd(indicesK, ykPtr, 8);
    const __m512d zk = remainderIsMasked ? _mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, zkPtr, 8)
                                         : _mm512_i64gather_pd(indicesK, zkPtr, 8);

    // only required for calculating globals
    const __m512i ownedStateK =
        remainderIsMasked
            ? _mm512_mask_i64gather_epi64(_ownedStateDummyMM512i, _masks[rest], indicesK,
                                          reinterpret_cast<const long long *const>(ownedStateKptr), 8)
            : _mm512_i64gather_epi64(indicesK, reinterpret_cast<const long long *const>(ownedStateKptr), 8);
    const __mmask8 ownedMaskK = _mm512_cmp_epi64_mask(ownedStateK, _ownedStateOwnedMM512i, _MM_CMPINT_EQ);

    // calculate distance j-k
    const __m512d drxjk = _mm512_sub_pd(xk, xj);
    const __m512d dryjk = _mm512_sub_pd(yk, yj);
    const __m512d drzjk = _mm512_sub_pd(zk, zj);

    const __m512d drxjk2 = _mm512_mul_pd(drxjk, drxjk);
    const __m512d drjk2PART = _mm512_fmadd_pd(dryjk, dryjk, drxjk2);
    const __m512d drjk2 = _mm512_fmadd_pd(drzjk, drzjk, drjk2PART);

    // calculate distance k-i
    const __m512d drxki = _mm512_sub_pd(xi, xk);
    const __m512d dryki = _mm512_sub_pd(yi, yk);
    const __m512d drzki = _mm512_sub_pd(zi, zk);

    const __m512d drxki2 = _mm512_mul_pd(drxki, drxki);
    const __m512d drki2PART = _mm512_fmadd_pd(dryki, dryki, drxki2);
    const __m512d drki2 = _mm512_fmadd_pd(drzki, drzki, drki2PART);

    // particles are presumed to be within cutoff -> calculate force
    const __m512d drxi2 = _mm512_mul_pd(drxij, drxki);
    const __m512d dri2PART = _mm512_fmadd_pd(dryij, dryki, drxi2);
    const __m512d dri2 = _mm512_fmadd_pd(drzij, drzki, dri2PART);

    const __m512d drxj2 = _mm512_mul_pd(drxij, drxjk);
    const __m512d drj2PART = _mm512_fmadd_pd(dryij, dryjk, drxj2);
    const __m512d drj2 = _mm512_fmadd_pd(drzij, drzjk, drj2PART);

    const __m512d drxk2 = _mm512_mul_pd(drxjk, drxki);
    const __m512d drk2PART = _mm512_fmadd_pd(dryjk, dryki, drxk2);
    const __m512d drk2 = _mm512_fmadd_pd(drzjk, drzki, drk2PART);

    const __m512d drijk2PART = _mm512_mul_pd(dri2, drj2);
    const __m512d drijk2 = _mm512_mul_pd(drijk2PART, drk2);

    const __m512d dr2PART = _mm512_mul_pd(drij2, drjk2);
    const __m512d dr2 = _mm512_mul_pd(dr2PART, drki2);
    const __m512d dr2sqrt = _mm512_sqrt_pd(dr2);
    const __m512d dr5PART = _mm512_mul_pd(dr2, dr2);
    const __m512d dr5 = _mm512_mul_pd(dr5PART, dr2sqrt);

    if constexpr (useMixing) {
      _nu = _mm512_set_pd(
          not remainderIsMasked or rest > 7 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[indicesK[7]]) : 0,
          not remainderIsMasked or rest > 6 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[indicesK[6]]) : 0,
          not remainderIsMasked or rest > 5 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[indicesK[5]]) : 0,
          not remainderIsMasked or rest > 4 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[indicesK[4]]) : 0,
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[indicesK[3]]) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[indicesK[2]]) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[indicesK[1]]) : 0,
          _PPLibrary->getMixingNu(typeI, typeJ, typeKptr[indicesK[0]]));
    }
    const __m512d invdr5 = remainderIsMasked ? _mm512_maskz_div_pd(_masks[rest], _nu, dr5) : _mm512_div_pd(_nu, dr5);
    const __m512d invdr53 = _mm512_mul_pd(_three, invdr5);

    const __m512d drj2_drk2 = _mm512_sub_pd(drj2, drk2);
    const __m512d i_jkFactor = _mm512_mul_pd(dri2, drj2_drk2);
    const __m512d drj2drk2 = _mm512_mul_pd(drj2, drk2);
    const __m512d drjk2drki2 = _mm512_mul_pd(drjk2, drki2);
    const __m512d drijk2BYdrij2 =
        remainderIsMasked ? _mm512_maskz_div_pd(_masks[rest], drijk2, drij2) : _mm512_div_pd(drijk2, drij2);
    const __m512d i_ijFactorPART = _mm512_sub_pd(drj2drk2, drjk2drki2);
    const __m512d i_ijFactor = _mm512_fmadd_pd(drijk2BYdrij2, _five, i_ijFactorPART);
    const __m512d drij2drjk2 = _mm512_mul_pd(drij2, drjk2);
    const __m512d drijk2BYdrki2 =
        remainderIsMasked ? _mm512_maskz_div_pd(_masks[rest], drijk2, drki2) : _mm512_div_pd(drijk2, drki2);
    const __m512d i_kiFactorPART = _mm512_sub_pd(drj2drk2, drij2drjk2);
    const __m512d i_kiFactor = _mm512_fmadd_pd(drijk2BYdrki2, _five, i_kiFactorPART);

    const __m512d fxiPART = _mm512_mul_pd(drxjk, i_jkFactor);
    const __m512d fxiPART2 = _mm512_fmadd_pd(drxij, i_ijFactor, fxiPART);
    const __m512d fxiPART3 = _mm512_fnmadd_pd(drxki, i_kiFactor, fxiPART2);
    const __m512d fxi = _mm512_mul_pd(invdr53, fxiPART3);
    fxiacc = _mm512_add_pd(fxi, fxiacc);
    const __m512d fyiPART = _mm512_mul_pd(dryjk, i_jkFactor);
    const __m512d fyiPART2 = _mm512_fmadd_pd(dryij, i_ijFactor, fyiPART);
    const __m512d fyiPART3 = _mm512_fnmadd_pd(dryki, i_kiFactor, fyiPART2);
    const __m512d fyi = _mm512_mul_pd(invdr53, fyiPART3);
    fyiacc = _mm512_add_pd(fyi, fyiacc);
    const __m512d fziPART = _mm512_mul_pd(drzjk, i_jkFactor);
    const __m512d fziPART2 = _mm512_fmadd_pd(drzij, i_ijFactor, fziPART);
    const __m512d fziPART3 = _mm512_fnmadd_pd(drzki, i_kiFactor, fziPART2);
    const __m512d fzi = _mm512_mul_pd(invdr53, fziPART3);
    fziacc = _mm512_add_pd(fzi, fziacc);

    if constexpr (calculateGlobals) {
      const __m512d virialXI = _mm512_mul_pd(fxi, xi);
      virialSumX = _mm512_mask_add_pd(virialSumX, ownedMaskI, virialSumX, virialXI);
      const __m512d virialYI = _mm512_mul_pd(fyi, yi);
      virialSumY = _mm512_mask_add_pd(virialSumY, ownedMaskI, virialSumY, virialYI);
      const __m512d virialZI = _mm512_mul_pd(fzi, zi);
      virialSumZ = _mm512_mask_add_pd(virialSumZ, ownedMaskI, virialSumZ, virialZI);
    }

    if constexpr (newton3j) {
      const __m512d drk2_dri2 = _mm512_sub_pd(drk2, dri2);
      const __m512d j_kiFactor = _mm512_mul_pd(drj2, drk2_dri2);
      const __m512d dri2drk2 = _mm512_mul_pd(dri2, drk2);
      const __m512d j_ijFactorPART = _mm512_sub_pd(dri2drk2, drjk2drki2);
      const __m512d j_ijFactor = _mm512_fmadd_pd(drijk2BYdrij2, _five, j_ijFactorPART);
      const __m512d drij2drki2 = _mm512_mul_pd(drij2, drki2);
      const __m512d drijk2BYdrjk2 =
          remainderIsMasked ? _mm512_maskz_div_pd(_masks[rest], drijk2, drjk2) : _mm512_div_pd(drijk2, drjk2);
      const __m512d j_jkFactorPART = _mm512_sub_pd(dri2drk2, drij2drki2);
      const __m512d j_jkFactor = _mm512_fmadd_pd(drijk2BYdrjk2, _five, j_jkFactorPART);

      const __m512d fxjPART = _mm512_mul_pd(drxki, j_kiFactor);
      const __m512d fxjPART2 = _mm512_fnmadd_pd(drxij, j_ijFactor, fxjPART);
      const __m512d fxjPART3 = _mm512_fmadd_pd(drxjk, j_jkFactor, fxjPART2);
      const __m512d fxj = _mm512_mul_pd(invdr53, fxjPART3);
      fxjacc = _mm512_add_pd(fxj, fxjacc);
      const __m512d fyjPART = _mm512_mul_pd(dryki, j_kiFactor);
      const __m512d fyjPART2 = _mm512_fnmadd_pd(dryij, j_ijFactor, fyjPART);
      const __m512d fyjPART3 = _mm512_fmadd_pd(dryjk, j_jkFactor, fyjPART2);
      const __m512d fyj = _mm512_mul_pd(invdr53, fyjPART3);
      fyjacc = _mm512_add_pd(fyj, fyjacc);
      const __m512d fzjPART = _mm512_mul_pd(drzki, j_kiFactor);
      const __m512d fzjPART2 = _mm512_fnmadd_pd(drzij, j_ijFactor, fzjPART);
      const __m512d fzjPART3 = _mm512_fmadd_pd(drzjk, j_jkFactor, fzjPART2);
      const __m512d fzj = _mm512_mul_pd(invdr53, fzjPART3);
      fzjacc = _mm512_add_pd(fzj, fzjacc);

      if constexpr (calculateGlobals) {
        const __m512d virialXJ = _mm512_mul_pd(fxj, xj);
        virialSumX = _mm512_mask_add_pd(virialSumX, ownedMaskJ, virialSumX, virialXJ);
        const __m512d virialYJ = _mm512_mul_pd(fyj, yj);
        virialSumY = _mm512_mask_add_pd(virialSumY, ownedMaskJ, virialSumY, virialYJ);
        const __m512d virialZJ = _mm512_mul_pd(fzj, zj);
        virialSumZ = _mm512_mask_add_pd(virialSumZ, ownedMaskJ, virialSumZ, virialZJ);
      }

      if constexpr (newton3k) {
        const __m512d nfxk = _mm512_add_pd(fxi, fxj);
        const __m512d nfyk = _mm512_add_pd(fyi, fyj);
        const __m512d nfzk = _mm512_add_pd(fzi, fzj);

        const __m512d fxk_old = remainderIsMasked ? _mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, fxkPtr, 8)
                                                  : _mm512_i64gather_pd(indicesK, fxkPtr, 8);
        const __m512d fyk_old = remainderIsMasked ? _mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, fykPtr, 8)
                                                  : _mm512_i64gather_pd(indicesK, fykPtr, 8);
        const __m512d fzk_old = remainderIsMasked ? _mm512_mask_i64gather_pd(_zero, _masks[rest], indicesK, fzkPtr, 8)
                                                  : _mm512_i64gather_pd(indicesK, fzkPtr, 8);

        const __m512d fxk_new = _mm512_sub_pd(fxk_old, nfxk);
        const __m512d fyk_new = _mm512_sub_pd(fyk_old, nfyk);
        const __m512d fzk_new = _mm512_sub_pd(fzk_old, nfzk);

        if constexpr (remainderIsMasked) {
          _mm512_mask_i64scatter_pd(fxkPtr, _masks[rest], indicesK, fxk_new, 8);
          _mm512_mask_i64scatter_pd(fykPtr, _masks[rest], indicesK, fyk_new, 8);
          _mm512_mask_i64scatter_pd(fzkPtr, _masks[rest], indicesK, fzk_new, 8);
        } else {
          _mm512_i64scatter_pd(fxkPtr, indicesK, fxk_new, 8);
          _mm512_i64scatter_pd(fykPtr, indicesK, fyk_new, 8);
          _mm512_i64scatter_pd(fzkPtr, indicesK, fzk_new, 8);
        }
        if constexpr (calculateGlobals) {
          const __m512d virialXK = _mm512_mul_pd(nfxk, xk);
          virialSumX = _mm512_mask_sub_pd(virialSumX, ownedMaskK, virialSumX, virialXK);
          const __m512d virialYK = _mm512_mul_pd(nfyk, yk);
          virialSumY = _mm512_mask_sub_pd(virialSumY, ownedMaskK, virialSumY, virialYK);
          const __m512d virialZK = _mm512_mul_pd(nfzk, zk);
          virialSumZ = _mm512_mask_sub_pd(virialSumZ, ownedMaskK, virialSumZ, virialZK);
        }
      }
    }

    if constexpr (calculateGlobals) {
      // Add 3 * potential energy to every owned particle of the interaction.
      // Division to the correct value is handled in endTraversal().
      // invdr53 is zero if particles are not within cutoff, so no cutoff mask is needed here.
      const __m512d potentialEnergy3 = _mm512_mul_pd(invdr53, _mm512_fnmadd_pd(dr2, _three, drijk2));
      potentialEnergySum = _mm512_mask_add_pd(potentialEnergySum, ownedMaskI, potentialEnergySum, potentialEnergy3);
      if constexpr (newton3j) {
        potentialEnergySum = _mm512_mask_add_pd(potentialEnergySum, ownedMaskJ, potentialEnergySum, potentialEnergy3);
      }
      if constexpr (newton3k) {
        potentialEnergySum = _mm512_mask_add_pd(potentialEnergySum, ownedMaskK, potentialEnergySum, potentialEnergy3);
      }
    }
  }

  /**
   * SoA kernel using simd compress and alignr instructions to compute a batch of particle indices relevant for
   * force calculation, which are then passed to SoAKernelGatherScatter.
   */
  template <bool newton3j, bool newton3k, bool remainderIsMasked>
  void SoAKernelCompressAlignr(__m512i &interactionIndices, int &numAssignedRegisters, const size_t k,
                               const __m512d &xi, const __m512d &yi, const __m512d &zi, const __m512d &xj,
                               const __m512d &yj, const __m512d &zj, const double *const xkPtr,
                               const double *const ykPtr, const double *const zkPtr, const size_t typeI,
                               const size_t typeJ, const size_t *const typeKptr, const __mmask8 &ownedMaskI,
                               const __mmask8 &ownedMaskJ, const autopas::OwnershipState *const ownedStateKptr,
                               const __m512d &drxij, const __m512d &dryij, const __m512d &drzij, const __m512d &drij2,
                               __m512d &fxiacc, __m512d &fyiacc, __m512d &fziacc, __m512d &fxjacc, __m512d &fyjacc,
                               __m512d &fzjacc, double *const fxkPtr, double *const fykPtr, double *const fzkPtr,
                               __m512d &potentialEnergySum, __m512d &virialSumX, __m512d &virialSumY,
                               __m512d &virialSumZ, unsigned int rest = 0) {
    const __m512i loopIndices = _mm512_add_epi64(_mm512_set1_epi64(k), _ascendingIndices);

    const __m512d xk = remainderIsMasked ? _mm512_maskz_loadu_pd(_masks[rest], &xkPtr[k]) : _mm512_loadu_pd(&xkPtr[k]);
    const __m512d yk = remainderIsMasked ? _mm512_maskz_loadu_pd(_masks[rest], &ykPtr[k]) : _mm512_loadu_pd(&ykPtr[k]);
    const __m512d zk = remainderIsMasked ? _mm512_maskz_loadu_pd(_masks[rest], &zkPtr[k]) : _mm512_loadu_pd(&zkPtr[k]);

    // calculate distance k-i
    const __m512d drxki = _mm512_sub_pd(xi, xk);
    const __m512d dryki = _mm512_sub_pd(yi, yk);
    const __m512d drzki = _mm512_sub_pd(zi, zk);

    const __m512d drxki2 = _mm512_mul_pd(drxki, drxki);
    const __m512d drki2PART = _mm512_fmadd_pd(dryki, dryki, drxki2);
    const __m512d drki2 = _mm512_fmadd_pd(drzki, drzki, drki2PART);

    const __mmask8 cutoffMask_ki = _mm512_cmp_pd_mask(drki2, _cutoffSquared, _CMP_LE_OS);

    // calculate distance j-k
    const __m512d drxjk = _mm512_sub_pd(xk, xj);
    const __m512d dryjk = _mm512_sub_pd(yk, yj);
    const __m512d drzjk = _mm512_sub_pd(zk, zj);

    const __m512d drxjk2 = _mm512_mul_pd(drxjk, drxjk);
    const __m512d drjk2PART = _mm512_fmadd_pd(dryjk, dryjk, drxjk2);
    const __m512d drjk2 = _mm512_fmadd_pd(drzjk, drzjk, drjk2PART);

    const __mmask8 cutoffMask_jk = _mm512_cmp_pd_mask(drjk2, _cutoffSquared, _CMP_LE_OS);

    const __mmask8 cutoffMask = _mm512_kand(cutoffMask_jk, cutoffMask_ki);

    const __m512i ownedStateK = _mm512_loadu_epi64(&ownedStateKptr[k]);
    const __mmask8 dummyMask = _mm512_cmp_epi64_mask(ownedStateK, _ownedStateDummyMM512i, _MM_CMPINT_NE);

    // mask with bits set if particle is not a dummy and within cutoff of the first two particles
    const __mmask8 mask = remainderIsMasked ? _mm512_kand(_masks[rest], _mm512_kand(cutoffMask, dummyMask))
                                            : _mm512_kand(cutoffMask, dummyMask);

    const int popCountMask = std::bitset<vecLength>(mask).count();

    if (popCountMask == 0) {
      // early exit if all particles are ignored
      return;
    }

    // compress indices of relevant particles towards lower end of vector register
    const __m512i newInteractionIndices = _mm512_maskz_compress_epi64(mask, loopIndices);
    if (numAssignedRegisters + popCountMask < vecLength) {
      // concatenate new and old indices and shift to higher end of vector register
      interactionIndices = _mm512_alignr_epi64(newInteractionIndices, interactionIndices, popCountMask);
      numAssignedRegisters += popCountMask;
    } else {
      // concatenate new and old indices and shift to higher end of vector register
      interactionIndices =
          _mm512_alignr_epi64(newInteractionIndices, interactionIndices, vecLength - numAssignedRegisters);

      SoAKernelGatherScatter<newton3j, newton3k, false>(
          interactionIndices, xi, yi, zi, xj, yj, zj, xkPtr, ykPtr, zkPtr, typeI, typeJ, typeKptr, ownedMaskI,
          ownedMaskJ, ownedStateKptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc,
          fxkPtr, fykPtr, fzkPtr, potentialEnergySum, virialSumX, virialSumY, virialSumZ);

      // shift remaining new indices to higher end of vector register
      interactionIndices = _mm512_alignr_epi64(newInteractionIndices, _zeroMM512i, popCountMask);
      numAssignedRegisters += popCountMask - vecLength;
    }
  }

  /**
   * Implementation function of SoAFunctorSingle(soa, newton3)
   *
   * @tparam newton3
   * @param soa
   */
  template <bool newton3>
  void SoAFunctorSingleImplCompressAlign(autopas::SoAView<SoAArraysType> soa) {
    if (soa.size() == 0) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    __m512d virialSumX = _mm512_setzero_pd();
    __m512d virialSumY = _mm512_setzero_pd();
    __m512d virialSumZ = _mm512_setzero_pd();
    __m512d potentialEnergySum = _mm512_setzero_pd();

    if constexpr (not useMixing) {
      _nu = _mm512_set1_pd(_nuAoS);
    }

    for (size_t i = soa.size() - 1; static_cast<long>(i) >= 2; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        // skip dummy particles
        continue;
      }

      // only required for calculating globals
      const __mmask8 ownedMaskI = ownedStatePtr[i] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

      const __m512d xi = _mm512_set1_pd(xptr[i]);
      const __m512d yi = _mm512_set1_pd(yptr[i]);
      const __m512d zi = _mm512_set1_pd(zptr[i]);

      __m512d fxiacc = _mm512_setzero_pd();
      __m512d fyiacc = _mm512_setzero_pd();
      __m512d fziacc = _mm512_setzero_pd();

      for (size_t j = i - 1; static_cast<long>(j) >= 1; --j) {
        if (ownedStatePtr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(xptr, yptr, zptr, xptr, yptr, zptr, i, j)) {
          // skip dummy particles
          continue;
        }

        // only required for calculating globals
        const __mmask8 ownedMaskJ = ownedStatePtr[j] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

        const __m512d xj = _mm512_set1_pd(xptr[j]);
        const __m512d yj = _mm512_set1_pd(yptr[j]);
        const __m512d zj = _mm512_set1_pd(zptr[j]);

        __m512d fxjacc = _mm512_setzero_pd();
        __m512d fyjacc = _mm512_setzero_pd();
        __m512d fzjacc = _mm512_setzero_pd();

        // calculate distance i-j
        const __m512d drxij = _mm512_sub_pd(xj, xi);
        const __m512d dryij = _mm512_sub_pd(yj, yi);
        const __m512d drzij = _mm512_sub_pd(zj, zi);

        const __m512d drxij2 = _mm512_mul_pd(drxij, drxij);
        const __m512d dryij2 = _mm512_mul_pd(dryij, dryij);
        const __m512d drzij2 = _mm512_mul_pd(drzij, drzij);

        const __m512d drij2PART = _mm512_add_pd(drxij2, dryij2);
        const __m512d drij2 = _mm512_add_pd(drij2PART, drzij2);

        __m512i interactionIndices = _zeroMM512i;
        int numAssignedRegisters = 0;

        size_t k = 0;
        // loop up to multiple of vecLength
        // only works if vecLength is a power of 2
        for (; k < (j & ~(vecLength - 1)); k += vecLength) {
          SoAKernelCompressAlignr<true, true, false>(interactionIndices, numAssignedRegisters, k, xi, yi, zi, xj, yj,
                                                     zj, xptr, yptr, zptr, typeptr[i], typeptr[j], typeptr, ownedMaskI,
                                                     ownedMaskJ, ownedStatePtr, drxij, dryij, drzij, drij2, fxiacc,
                                                     fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fxptr, nullptr, fzptr,
                                                     potentialEnergySum, virialSumX, virialSumY, virialSumZ);
        }
        // process remainder
        if (k < j) {
          SoAKernelCompressAlignr<true, true, true>(interactionIndices, numAssignedRegisters, k, xi, yi, zi, xj, yj, zj,
                                                    xptr, yptr, zptr, typeptr[i], typeptr[j], typeptr, ownedMaskI,
                                                    ownedMaskJ, ownedStatePtr, drxij, dryij, drzij, drij2, fxiacc,
                                                    fyiacc, fziacc, fxjacc, fyjacc, fzjacc, fxptr, nullptr, fzptr,
                                                    potentialEnergySum, virialSumX, virialSumY, virialSumZ, j - k);
        }
        // process remaining interaction indices
        if (numAssignedRegisters > 0) {
          // shift remaining interaction indices to right
          interactionIndices = _mm512_alignr_epi64(_zeroMM512i, interactionIndices, 8 - numAssignedRegisters);
          SoAKernelGatherScatter<true, true, true>(
              interactionIndices, xi, yi, zi, xj, yj, zj, xptr, yptr, zptr, typeptr[i], typeptr[j], typeptr, ownedMaskI,
              ownedMaskJ, ownedStatePtr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc, fyjacc, fzjacc,
              fxptr, fyptr, fzptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numAssignedRegisters);
        }

        fxptr[j] += _mm512_reduce_add_pd(fxjacc);
        fyptr[j] += _mm512_reduce_add_pd(fyjacc);
        fzptr[j] += _mm512_reduce_add_pd(fzjacc);
      }
      fxptr[i] += _mm512_reduce_add_pd(fxiacc);
      fyptr[i] += _mm512_reduce_add_pd(fyiacc);
      fzptr[i] += _mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      _aosThreadDataGlobals[threadnum].potentialEnergySum += _mm512_reduce_add_pd(potentialEnergySum);
      _aosThreadDataGlobals[threadnum].virialSum[0] += _mm512_reduce_add_pd(virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += _mm512_reduce_add_pd(virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += _mm512_reduce_add_pd(virialSumZ);
    }
  }

  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImplCompressAlign(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    __m512d potentialEnergySum = _mm512_setzero_pd();
    __m512d virialSumX = _mm512_setzero_pd();
    __m512d virialSumY = _mm512_setzero_pd();
    __m512d virialSumZ = _mm512_setzero_pd();

    if constexpr (not useMixing) {
      _nu = _mm512_set1_pd(_nuAoS);
    }

    // particle 1 always from soa1
    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        // skip dummy particles
        continue;
      }

      // only required for calculating globals
      const __mmask8 ownedMaskI = ownedState1ptr[i] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

      const __m512d xi = _mm512_set1_pd(x1ptr[i]);
      const __m512d yi = _mm512_set1_pd(y1ptr[i]);
      const __m512d zi = _mm512_set1_pd(z1ptr[i]);

      __m512d fxiacc = _mm512_setzero_pd();
      __m512d fyiacc = _mm512_setzero_pd();
      __m512d fziacc = _mm512_setzero_pd();

      // particle 2 from soa1 and 3 from soa2
      for (size_t j = i + 1; j < soa1.size(); ++j) {
        if (ownedState1ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x1ptr, y1ptr, z1ptr, i, j)) {
          // skip dummy particles
          continue;
        }

        // only required for calculating globals
        const __mmask8 ownedMaskJ = ownedState1ptr[j] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

        const __m512d xj = _mm512_set1_pd(x1ptr[j]);
        const __m512d yj = _mm512_set1_pd(y1ptr[j]);
        const __m512d zj = _mm512_set1_pd(z1ptr[j]);

        __m512d fxjacc = _mm512_setzero_pd();
        __m512d fyjacc = _mm512_setzero_pd();
        __m512d fzjacc = _mm512_setzero_pd();

        // calculate distance i-j
        const __m512d drxij = _mm512_sub_pd(xj, xi);
        const __m512d dryij = _mm512_sub_pd(yj, yi);
        const __m512d drzij = _mm512_sub_pd(zj, zi);

        const __m512d drxij2 = _mm512_mul_pd(drxij, drxij);
        const __m512d dryij2 = _mm512_mul_pd(dryij, dryij);
        const __m512d drzij2 = _mm512_mul_pd(drzij, drzij);

        const __m512d drij2PART = _mm512_add_pd(drxij2, dryij2);
        const __m512d drij2 = _mm512_add_pd(drij2PART, drzij2);

        __m512i interactionIndices = _zeroMM512i;
        int numAssignedRegisters = 0;

        size_t k = 0;
        // loop up to multiple of vecLength
        // only works if vecLength is a power of 2
        for (; k < (soa2.size() & ~(vecLength - 1)); k += vecLength) {
          SoAKernelCompressAlignr<true, newton3, false>(
              interactionIndices, numAssignedRegisters, k, xi, yi, zi, xj, yj, zj, x2ptr, y2ptr, z2ptr, type1ptr[i],
              type1ptr[j], type2ptr, ownedMaskI, ownedMaskJ, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc,
              fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, nullptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY,
              virialSumZ);
        }
        // process remainder
        if (k < soa2.size()) {
          SoAKernelCompressAlignr<true, newton3, true>(
              interactionIndices, numAssignedRegisters, k, xi, yi, zi, xj, yj, zj, x2ptr, y2ptr, z2ptr, type1ptr[i],
              type1ptr[j], type2ptr, ownedMaskI, ownedMaskJ, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc,
              fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, nullptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY,
              virialSumZ, soa2.size() - k);
        }
        if (numAssignedRegisters > 0) {
          // shift remainder to right
          interactionIndices = _mm512_alignr_epi64(_zeroMM512i, interactionIndices, 8 - numAssignedRegisters);
          SoAKernelGatherScatter<true, newton3, true>(
              interactionIndices, xi, yi, zi, xj, yj, zj, x2ptr, y2ptr, z2ptr, type1ptr[i], type1ptr[j], type2ptr,
              ownedMaskI, ownedMaskJ, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc,
              fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ,
              numAssignedRegisters);
        }
        fx1ptr[j] += _mm512_reduce_add_pd(fxjacc);
        fy1ptr[j] += _mm512_reduce_add_pd(fyjacc);
        fz1ptr[j] += _mm512_reduce_add_pd(fzjacc);
      }

      // both particles 2 and 3 from soa2

      for (size_t j = soa2.size() - 1; static_cast<long>(j) >= 0; --j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, j)) {
          continue;
        }

        // only required for calculating globals
        const __mmask8 ownedMask2 = ownedState2ptr[j] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

        const __m512d x2 = _mm512_set1_pd(x2ptr[j]);
        const __m512d y2 = _mm512_set1_pd(y2ptr[j]);
        const __m512d z2 = _mm512_set1_pd(z2ptr[j]);

        __m512d fxjacc = _mm512_setzero_pd();
        __m512d fyjacc = _mm512_setzero_pd();
        __m512d fzjacc = _mm512_setzero_pd();

        const __m512d drxij = _mm512_sub_pd(x2, xi);
        const __m512d dryij = _mm512_sub_pd(y2, yi);
        const __m512d drzij = _mm512_sub_pd(z2, zi);

        const __m512d drxij2 = _mm512_mul_pd(drxij, drxij);
        const __m512d dryij2 = _mm512_mul_pd(dryij, dryij);
        const __m512d drzij2 = _mm512_mul_pd(drzij, drzij);

        const __m512d drij2PART = _mm512_add_pd(drxij2, dryij2);
        const __m512d drij2 = _mm512_add_pd(drij2PART, drzij2);

        // particle 3 from soa 2

        __m512i interactionIndices = _zeroMM512i;
        int numAssignedRegisters = 0;

        size_t k = 0;
        for (; k < (j & ~(vecLength - 1)); k += vecLength) {
          SoAKernelCompressAlignr<newton3, newton3, false>(
              interactionIndices, numAssignedRegisters, k, xi, yi, zi, x2, y2, z2, x2ptr, y2ptr, z2ptr, type1ptr[i],
              type2ptr[j], type2ptr, ownedMaskI, ownedMask2, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc,
              fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, nullptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY,
              virialSumZ);
        }
        if (k < j) {
          SoAKernelCompressAlignr<newton3, newton3, true>(
              interactionIndices, numAssignedRegisters, k, xi, yi, zi, x2, y2, z2, x2ptr, y2ptr, z2ptr, type1ptr[i],
              type2ptr[j], type2ptr, ownedMaskI, ownedMask2, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc,
              fziacc, fxjacc, fyjacc, fzjacc, fx2ptr, nullptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY,
              virialSumZ, j - k);
        }
        // process remaining interaction indices
        if (numAssignedRegisters > 0) {
          // shift remaining interaction indices to right
          interactionIndices = _mm512_alignr_epi64(_zeroMM512i, interactionIndices, 8 - numAssignedRegisters);
          SoAKernelGatherScatter<newton3, newton3, true>(
              interactionIndices, xi, yi, zi, x2, y2, z2, x2ptr, y2ptr, z2ptr, type1ptr[i], type2ptr[j], type2ptr,
              ownedMaskI, ownedMask2, ownedState2ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc,
              fyjacc, fzjacc, fx2ptr, fy2ptr, fz2ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ,
              numAssignedRegisters);
        }

        if constexpr (newton3) {
          fx2ptr[j] += _mm512_reduce_add_pd(fxjacc);
          fy2ptr[j] += _mm512_reduce_add_pd(fyjacc);
          fz2ptr[j] += _mm512_reduce_add_pd(fzjacc);
        }
      }

      fx1ptr[i] += _mm512_reduce_add_pd(fxiacc);
      fy1ptr[i] += _mm512_reduce_add_pd(fyiacc);
      fz1ptr[i] += _mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      _aosThreadDataGlobals[threadnum].potentialEnergySum += _mm512_reduce_add_pd(potentialEnergySum);
      _aosThreadDataGlobals[threadnum].virialSum[0] += _mm512_reduce_add_pd(virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += _mm512_reduce_add_pd(virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += _mm512_reduce_add_pd(virialSumZ);
    }
  }

  /**
   * Implementation function of SoAFunctorTriple(soa1, soa2, soa3, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   * @param soa3
   */
  template <bool newton3>
  void SoAFunctorTripleImplCompressAlign(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                         autopas::SoAView<SoAArraysType> soa3) {
    if (soa1.size() == 0 or soa2.size() == 0 or soa3.size() == 0) {
      return;
    }

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x3ptr = soa3.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y3ptr = soa3.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z3ptr = soa3.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedState1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedState3ptr = soa3.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx3ptr = soa3.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy3ptr = soa3.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz3ptr = soa3.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict type1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type2ptr = soa2.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict type3ptr = soa3.template begin<Particle::AttributeNames::typeId>();

    __m512d potentialEnergySum = _mm512_setzero_pd();
    __m512d virialSumX = _mm512_setzero_pd();
    __m512d virialSumY = _mm512_setzero_pd();
    __m512d virialSumZ = _mm512_setzero_pd();

    if constexpr (not useMixing) {
      _nu = _mm512_set1_pd(_nuAoS);
    }

    for (size_t i = 0; i < soa1.size(); ++i) {
      if (ownedState1ptr[i] == autopas::OwnershipState::dummy) {
        // skip dummy particles
        continue;
      }

      // only required for calculating globals
      const __mmask8 ownedMask1 = ownedState1ptr[i] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

      const __m512d xi = _mm512_set1_pd(x1ptr[i]);
      const __m512d yi = _mm512_set1_pd(y1ptr[i]);
      const __m512d zi = _mm512_set1_pd(z1ptr[i]);

      __m512d fxiacc = _mm512_setzero_pd();
      __m512d fyiacc = _mm512_setzero_pd();
      __m512d fziacc = _mm512_setzero_pd();

      for (size_t j = 0; j < soa2.size(); ++j) {
        if (ownedState2ptr[j] == autopas::OwnershipState::dummy or
            not SoAParticlesInCutoff(x1ptr, y1ptr, z1ptr, x2ptr, y2ptr, z2ptr, i, j)) {
          // skip dummy particles
          continue;
        }

        // only required for calculating globals
        const __mmask8 ownedMaskJ = ownedState2ptr[j] == autopas::OwnershipState::owned ? 0b11111111 : 0b00000000;

        const __m512d xj = _mm512_set1_pd(x2ptr[j]);
        const __m512d yj = _mm512_set1_pd(y2ptr[j]);
        const __m512d zj = _mm512_set1_pd(z2ptr[j]);

        __m512d fxjacc = _mm512_setzero_pd();
        __m512d fyjacc = _mm512_setzero_pd();
        __m512d fzjacc = _mm512_setzero_pd();

        // calculate distance i-j
        const __m512d drxij = _mm512_sub_pd(xj, xi);
        const __m512d dryij = _mm512_sub_pd(yj, yi);
        const __m512d drzij = _mm512_sub_pd(zj, zi);

        const __m512d drxij2 = _mm512_mul_pd(drxij, drxij);
        const __m512d dryij2 = _mm512_mul_pd(dryij, dryij);
        const __m512d drzij2 = _mm512_mul_pd(drzij, drzij);

        const __m512d drij2PART = _mm512_add_pd(drxij2, dryij2);
        const __m512d drij2 = _mm512_add_pd(drij2PART, drzij2);

        __m512i interactionIndices = _zeroMM512i;
        int numAssignedRegisters = 0;

        size_t k = 0;
        // loop up to multiple of vecLength
        // only works if vecLength is a power of 2
        for (; k < (soa3.size() & ~(vecLength - 1)); k += vecLength) {
          SoAKernelCompressAlignr<newton3, newton3, false>(
              interactionIndices, numAssignedRegisters, k, xi, yi, zi, xj, yj, zj, x3ptr, y3ptr, z3ptr, type1ptr[i],
              type2ptr[j], type3ptr, ownedMask1, ownedMaskJ, ownedState3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc,
              fziacc, fxjacc, fyjacc, fzjacc, fx3ptr, nullptr, fz3ptr, potentialEnergySum, virialSumX, virialSumY,
              virialSumZ);
        }
        // process remainder
        if (k < soa3.size()) {
          SoAKernelCompressAlignr<newton3, newton3, true>(
              interactionIndices, numAssignedRegisters, k, xi, yi, zi, xj, yj, zj, x3ptr, y3ptr, z3ptr, type1ptr[i],
              type2ptr[j], type3ptr, ownedMask1, ownedMaskJ, ownedState3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc,
              fziacc, fxjacc, fyjacc, fzjacc, fx3ptr, nullptr, fz3ptr, potentialEnergySum, virialSumX, virialSumY,
              virialSumZ, soa3.size() - k);
        }
        // process remaining interaction indices
        if (numAssignedRegisters > 0) {
          // shift remaining interaction indices to right
          interactionIndices = _mm512_alignr_epi64(_zeroMM512i, interactionIndices, 8 - numAssignedRegisters);
          SoAKernelGatherScatter<newton3, newton3, true>(
              interactionIndices, xi, yi, zi, xj, yj, zj, x3ptr, y3ptr, z3ptr, type1ptr[i], type2ptr[j], type3ptr,
              ownedMask1, ownedMaskJ, ownedState3ptr, drxij, dryij, drzij, drij2, fxiacc, fyiacc, fziacc, fxjacc,
              fyjacc, fzjacc, fx3ptr, fy3ptr, fz3ptr, potentialEnergySum, virialSumX, virialSumY, virialSumZ,
              numAssignedRegisters);
        }

        if constexpr (newton3) {
          fx2ptr[j] += _mm512_reduce_add_pd(fxjacc);
          fy2ptr[j] += _mm512_reduce_add_pd(fyjacc);
          fz2ptr[j] += _mm512_reduce_add_pd(fzjacc);
        }
      }

      fx1ptr[i] += _mm512_reduce_add_pd(fxiacc);
      fy1ptr[i] += _mm512_reduce_add_pd(fyiacc);
      fz1ptr[i] += _mm512_reduce_add_pd(fziacc);
    }
    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      _aosThreadDataGlobals[threadnum].potentialEnergySum += _mm512_reduce_add_pd(potentialEnergySum);
      _aosThreadDataGlobals[threadnum].virialSum[0] += _mm512_reduce_add_pd(virialSumX);
      _aosThreadDataGlobals[threadnum].virialSum[1] += _mm512_reduce_add_pd(virialSumY);
      _aosThreadDataGlobals[threadnum].virialSum[2] += _mm512_reduce_add_pd(virialSumZ);
    }
  }

 public:
  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param nu The Axilrod-Teller potential parameter
   */
  void setParticleProperties(SoAFloatPrecision nu) { _nuAoS = nu; }

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
   * - virial: 3 without n3, 9 with n3
   * - accumulation: 4 without n3, 12 with n3
   * - Total: 10 without n3, 24 with n3
   *
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

      constexpr size_t numFLOPsPerDistanceCall = 24;
      constexpr size_t numFLOPsPerN3KernelCall = 100;
      constexpr size_t numFLOPsPerNoN3KernelCall = 59;
      constexpr size_t numFLOPsPerN3GlobalCalc = 24;
      constexpr size_t numFLOPsPerNoN3GlobalCalc = 10;

      return numDistCallsAcc * numFLOPsPerDistanceCall + numKernelCallsN3Acc * numFLOPsPerN3KernelCall +
             numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall + numGlobalCalcsN3Acc * numFLOPsPerN3GlobalCalc +
             numGlobalCalcsNoN3Acc * numFLOPsPerNoN3GlobalCalc;
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<size_t>::max();
    }
  }

  [[nodiscard]] double getHitRate() const override {
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

      return (static_cast<double>(numKernelCallsNoN3Acc) + static_cast<double>(numKernelCallsN3Acc)) /
             (static_cast<double>(numDistCallsAcc));
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

 private:
  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    autopas::utils::ExceptionHandler::exception(
        "AxilrodTellerFunctorAVX512::SoAFunctorVerletImpl() is not implemented.");
  }

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
    double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
  };

  // make sure of the size of AoSThreadDataGlobals
  static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
  static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

  const double _cutoffSquaredAoS;
  const __m512d _cutoffSquared;

  // Parameter of the Axilrod-Teller potential
  // not const because they might be reset through PPL
  double _nuAoS = 0.0;
  __m512d _nu;

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

  const __m512d _zero{_mm512_setzero_pd()};
  const __m512i _zeroMM512i{_mm512_setzero_si512()};
  const __m512d _three{_mm512_set1_pd(3.0)};
  const __m512d _five{_mm512_set1_pd(5.0)};

  static constexpr int vecLength = 8;
  const __mmask8 _masks[vecLength]{0b00000000, 0b00000001, 0b00000011, 0b00000111,
                                   0b00001111, 0b00011111, 0b00111111, 0b01111111};

  __m512i _ascendingIndices{_mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)};

  const __m512i _ownedStateDummyMM512i{_mm512_set1_epi64(static_cast<int64_t>(autopas::OwnershipState::dummy))};
  const __m512i _ownedStateOwnedMM512i{_mm512_set1_epi64(static_cast<int64_t>(autopas::OwnershipState::owned))};
};
}  // namespace mdLib

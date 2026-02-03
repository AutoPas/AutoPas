/**
 * @file EmptyPairwiseFunctor.h
 * @author seckler
 * @date 26.03.20
 */

#pragma once

#include "autopas/baseFunctors/Functor.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

/**
 * Empty Functor, this functor is empty and can be used for testing purposes.
 * It returns that it is applicable for everything.
 */
template <class Particle_T>
class EmptyPairwiseFunctor : public autopas::PairwiseFunctor<Particle_T, EmptyPairwiseFunctor<Particle_T>> {
 private:
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Default constructor.
   */
  EmptyPairwiseFunctor() : autopas::PairwiseFunctor<Particle_T, EmptyPairwiseFunctor<Particle_T>>(0.){};

  /**
   * @copydoc autopas::PairwiseFunctor::AoSFunctor()
   */
  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) override {}

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   */
  void SoAFunctorSingle(autopas::SoAView<typename Particle_T::SoAArraysType> soa, bool newton3) override {}

  /**
   * SoAFunctor for a pair of SoAs.
   * @param soa1 An autopas::SoAView for the Functor
   * @param soa2 A second autopas::SoAView for the Functor
   * @param newton3 A boolean to indicate whether to allow newton3
   */
  void SoAFunctorPair(autopas::SoAView<typename Particle_T::SoAArraysType> soa1,
                      autopas::SoAView<typename Particle_T::SoAArraysType> soa2, bool newton3) override {}

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   */
  void SoAFunctorVerlet(autopas::SoAView<typename Particle_T::SoAArraysType> soa, size_t indexFirst,
                        const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &neighborList,
                        bool newton3) override{};

  /**
   * @copydoc autopas::Functor::allowsNewton3()
   */
  bool allowsNewton3() override { return true; }

  /**
   * @copydoc autopas::Functor::allowsNonNewton3()
   */
  bool allowsNonNewton3() override { return true; }

  /**
   * @copydoc autopas::Functor::getName()
   */
  std::string getName() override { return "EmptyPairwiseFunctor"; }

  /**
   * @copydoc autopas::Functor::isRelevantForTuning()
   */
  bool isRelevantForTuning() override { return true; }
};

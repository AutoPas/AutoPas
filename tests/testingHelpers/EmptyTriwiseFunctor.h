/**
 * @file EmptyTriwiseFunctor.h
 * @author S. J. Newcome
 * @date 08/10/2026
 */

#pragma once

#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/options/DataLayoutOption.h"

/**
 * Empty Functor, this functor is empty and can be used for testing purposes.
 * It returns that it is applicable for everything.
 */
template <class Particle_T>
class EmptyTriwiseFunctor : public autopas::TriwiseFunctor<Particle_T, EmptyTriwiseFunctor<Particle_T>> {
 private:
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Default constructor.
   */
  EmptyTriwiseFunctor() : autopas::TriwiseFunctor<Particle_T, EmptyTriwiseFunctor<Particle_T>>(0.){};

  /**
   * @copydoc autopas::TriwiseFunctor::AoSFunctor()
   */
  void AoSFunctor(Particle_T &i, Particle_T &j, Particle_T &k, bool newton3) override {}

  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorSingle()
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
   * SoAFunctor for a triple of SoAs.
   * @param soa1 An autopas::SoAView for the Functor
   * @param soa2 A second autopas::SoAView for the Functor
   * @param soa3 A third autopas::SoAView for the Functor
   * @param newton3 A boolean to indicate whether to allow newton3
   */
  void SoAFunctorTriple(autopas::SoAView<typename Particle_T::SoAArraysType> soa1,
                        autopas::SoAView<typename Particle_T::SoAArraysType> soa2,
                        autopas::SoAView<typename Particle_T::SoAArraysType> soa3, bool newton3) override {}

  /**
   * @copydoc autopas::TriwiseFunctor::SoAFunctorVerlet()
   */
  void SoAFunctorVerlet(autopas::SoAView<typename Particle_T::SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
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
  std::string getName() override { return "EmptyTriwiseFunctor"; }

  /**
   * @copydoc autopas::Functor::isRelevantForTuning()
   */
  bool isRelevantForTuning() override { return true; }
};
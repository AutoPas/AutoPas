/**
 * @file checkFunctorType.h
 *
 * @date 28 Aug 2023
 * @author muehlhaeusser
 */

#pragma once

#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"

namespace {
/**
 * Returns false for any types that are not a PairwiseFunctor
 * @return false
 */
std::false_type isPairwiseFunctorImpl(...);
/**
 * Return true for a PairwiseFunctor type
 * @tparam Particle_T
 * @tparam FunctorT
 * @return true
 */
template <typename Particle_T, typename FunctorT>
std::true_type isPairwiseFunctorImpl(const autopas::PairwiseFunctor<Particle_T, FunctorT> &);

/**
 * Returns false for any types that are not a TriwiseFunctor
 * @return false
 */
std::false_type isTriwiseFunctorImpl(...);
/**
 * Return true for a TriwiseFunctor type
 * @tparam Particle_T
 * @tparam FunctorT
 * @return true
 */
template <typename Particle_T, typename FunctorT>
std::true_type isTriwiseFunctorImpl(const autopas::TriwiseFunctor<Particle_T, FunctorT> &);
}  // namespace

namespace autopas::utils {
/**
 * Check whether a Functor Type is inheriting from PairwiseFunctor
 * @tparam FunctorT
 */
template <typename FunctorT>
using isPairwiseFunctor = decltype(isPairwiseFunctorImpl(std::declval<FunctorT>()));

/**
 * Check whether a Functor Type is inheriting from TriwiseFunctor
 * @tparam FunctorT
 */
template <typename FunctorT>
using isTriwiseFunctor = decltype(isTriwiseFunctorImpl(std::declval<FunctorT>()));

}  // namespace autopas::utils

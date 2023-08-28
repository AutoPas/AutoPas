/**
 * @file checkFunctorType.h
 *
 * @date 28 Aug 2023
 * @author muehlhaeusser
 */

#pragma once

#include "autopas/pairwiseFunctors/PairwiseFunctor.h"
#include "autopas/pairwiseFunctors/TriwiseFunctor.h"


namespace {
std::false_type isPairwiseFunctorImpl(...);
template <typename T, typename U>
std::true_type isPairwiseFunctorImpl(autopas::PairwiseFunctor<T,U> const volatile&);

std::false_type isTriwiseFunctorImpl(...);
template <typename T, typename U>
std::true_type isTriwiseFunctorImpl(autopas::TriwiseFunctor<T,U> const volatile&);
}

namespace autopas::utils {

template <typename T>
using isPairwiseFunctor = decltype(isPairwiseFunctorImpl(std::declval<T&>()));

template <typename T>
using isTriwiseFunctor = decltype(isTriwiseFunctorImpl(std::declval<T&>()));

}  // namespace autopas::utils



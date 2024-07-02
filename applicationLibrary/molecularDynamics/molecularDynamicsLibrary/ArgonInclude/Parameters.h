/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 13/06/24
 */

#pragma once

#include "autopas/utils/WrapOpenMP.h"

namespace mdLib::Argon {

/**
 *
 * @param i first index
 * @param j second index
 * @param k third index
 * @return returns the position of the element in the parameter array corresponding to the triplet (i, j, k), to be
 * used for parameters A and alpha (repulsive part)
 */
constexpr size_t indexRepulsivePart(const size_t i, const size_t j, const size_t k) {
  if (i == 0 && j == 0 && k == 0) {
    return 0;
  }
  if (i == 0 && j == 0 && k == 1) {
    return 1;
  }
  if (i == 0 && j == 1 && k == 1) {
    return 2;
  }
  if (i == 1 && j == 1 && k == 1) {
    return 3;
  }
  if (i == 0 && j == 0 && k == 2) {
    return 4;
  }
  if (i == 0 && j == 1 && k == 2) {
    return 5;
  }
  if (i == 1 && j == 1 && k == 2) {
    return 6;
  }
  if (i == 0 && j == 2 && k == 2) {
    return 7;
  }
  if (i == 1 && j == 2 && k == 2) {
    return 8;
  }
  if (i == 2 && j == 2 && k == 2) {
    return 9;
  }
  if (i == 0 && j == 0 && k == 3) {
    return 10;
  }
  if (i == 0 && j == 1 && k == 3) {
    return 11;
  }
  if (i == 1 && j == 1 && k == 3) {
    return 12;
  }
  if (i == 0 && j == 2 && k == 3) {
    return 13;
  }
  if (i == 1 && j == 2 && k == 3) {
    return 14;
  }
  if (i == 0 && j == 3 && k == 3) {
    return 15;
  }
  if (i == 0 && j == 0 && k == 4) {
    return 16;
  }
  if (i == 0 && j == 1 && k == 4) {
    return 17;
  }
  if (i == 1 && j == 1 && k == 4) {
    return 18;
  }
  if (i == 0 && j == 2 && k == 4) {
    return 19;
  }
  if (i == 0 && j == 0 && k == 5) {
    return 20;
  }
  if (i == 0 && j == 1 && k == 5) {
    return 21;
  }
  if (i == 0 && j == 0 && k == 6) {
    return 22;
  }
  throw autopas::utils::ExceptionHandler::AutoPasException("Parameter cannot be accessed");
}

/**
 *
 * @param i first index
 * @param j second index
 * @param k third index
 * @return returns the position of the element in the parameter array corresponding to the triplet (i, j, k), to be
 * used for parameters Z and beta (dispersion part)
 */
constexpr size_t indexDispersionPart(const size_t i, const size_t j, const size_t k) {
  if (i == 1 && j == 1 && k == 1) {
    return 0;
  }
  if (i == 1 && j == 1 && k == 2) {
    return 1;
  }
  if (i == 1 && j == 2 && k == 2) {
    return 2;
  }
  if (i == 2 && j == 2 && k == 2) {
    return 3;
  }
  if (i == 1 && j == 1 && k == 3) {
    return 4;
  }
  throw autopas::utils::ExceptionHandler::AutoPasException("Parameter cannot be accessed");
}

enum param { A, alpha, Z, beta };

/**
 *
 * @tparam P the parameter : either A, alpha, Z or beta
 * @param i first index
 * @param j second index
 * @param k third index
 * @return returns the position of the element in the parameter array corresponding to the triplet (i, j, k)
 */
template <param P>
constexpr size_t index(const size_t i, const size_t j, const size_t k);

template <>
constexpr size_t index<A>(const size_t i, const size_t j, const size_t k) {
  return indexRepulsivePart(i, j, k);
}

template <>
constexpr size_t index<alpha>(const size_t i, const size_t j, const size_t k) {
  return indexRepulsivePart(i, j, k);
}

template <>
constexpr size_t index<Z>(const size_t i, const size_t j, const size_t k) {
  return indexDispersionPart(i, j, k);
}

template <>
constexpr size_t index<beta>(const size_t i, const size_t j, const size_t k) {
  return indexDispersionPart(i, j, k);
}

}  // namespace mdLib::Argon
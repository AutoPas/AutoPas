/**
 * @file MoleculeLJ.cpp
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */
#include "autopas/particles/MoleculeLJ.h"

namespace autopas {
template <typename floatType>
double MoleculeLJBase<floatType>::getEpsilon() {
  return EPSILON;
}

template <typename floatType>
void MoleculeLJBase<floatType>::setEpsilon(double epsilon) {
  EPSILON = epsilon;
}

template <typename floatType>
double MoleculeLJBase<floatType>::getSigma() {
  return SIGMA;
}

template <typename floatType>
void MoleculeLJBase<floatType>::setSigma(double sigma) {
  SIGMA = sigma;
}

template <typename floatType>
double MoleculeLJBase<floatType>::EPSILON;
template <typename floatType>
double MoleculeLJBase<floatType>::SIGMA;

template class MoleculeLJBase<float>;
template class MoleculeLJBase<double>;
}  // namespace autopas

/**
 * @file MoleculeLJ.cpp
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */
#include "autopas/particles/MoleculeLJ.h"

namespace autopas {

template <typename floatType>
floatType MoleculeLJBase<floatType>::EPSILON;
template <typename floatType>
floatType MoleculeLJBase<floatType>::SIGMA;

template class MoleculeLJBase<float>;
template class MoleculeLJBase<double>;
}  // namespace autopas

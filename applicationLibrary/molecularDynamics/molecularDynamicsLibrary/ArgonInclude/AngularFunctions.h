/**
* @file ArgonFunctor.h
* @author I. Angelucci
* @date 19/06/24
 */

#pragma once

#include "DisplacementHandle.h"
#include "CosineHandle.h"

namespace autopas::utils::ArrayMath::Argon {

  template <size_t a, size_t b, size_t c>
  [[nodiscard]] double AngularTerm(const DisplacementHandle& displacementIJ, const DisplacementHandle& displacementJK, const DisplacementHandle& displacementKI,
                             const CosineHandle& cosineI, const CosineHandle& cosineJ, const CosineHandle& cosineK);

  template <size_t a, size_t b, size_t c>
  [[nodiscard]] double AngularTerm_derive_wrt(const DisplacementHandle& displacementIJ, const DisplacementHandle& displacementJK, const DisplacementHandle& displacementKI,
                         const CosineHandle& cosineI, const CosineHandle& cosineJ, const CosineHandle& cosineK);

} //  namespace autopas::utils::ArrayMath::Argon
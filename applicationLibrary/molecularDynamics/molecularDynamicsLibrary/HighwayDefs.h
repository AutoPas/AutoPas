/**
 * @file HighwayDefs.h
 * @author D. Martin
 * @date 13/11/25
 */

#pragma once
#include "autopas/options/VectorizationPatternOption.h"
#include "hwy/highway.h"

namespace mdLib {

namespace highway = hwy::HWY_NAMESPACE;
/** Higwhay tag for full double register */
inline const highway::ScalableTag<double> tag_double;
/** Highway tag for full long register */
inline const highway::ScalableTag<int64_t> tag_long;
/** Highway tag for full index register */
inline const highway::ScalableTag<size_t> tag_size_t;
/** Nuber of double values in a full register */
inline const size_t _vecLengthDouble{highway::Lanes(tag_double)};
/** Type for a Double vector register */
using VectorDouble = decltype(highway::Zero(tag_double));
/** Type for a Long vector register */
using VectorLong = decltype(highway::Zero(tag_long));
/** Type for a size_t vector register */
using VectorSizeT = decltype(highway::Zero(tag_size_t));
/** Highway tag for a half-filled double register */
inline const highway::Half<highway::DFromV<VectorDouble>> tag_double_half;
/** Type for a Double Mask */
using MaskDouble = decltype(highway::FirstN(tag_double, 1));
/** Type for a Long Mask */
using MaskLong = decltype(highway::FirstN(tag_long, 2));
/** Vectorization Pattern Type */
using VectorizationPattern = autopas::VectorizationPatternOption::Value;

}  // namespace mdLib

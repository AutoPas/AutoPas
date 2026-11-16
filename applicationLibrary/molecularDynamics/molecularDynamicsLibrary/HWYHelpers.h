/**
 * @file HWYHelpers.h
 *
 * Shared Google Highway tag aliases and horizontal reduction helpers for molecularDynamicsLibrary functors.
 *
 * @date 15.11.2026
 * @author muehlhaeusser
 */

#pragma once

#include <hwy/highway.h>

#include <cstddef>
#include <cstdint>

namespace mdLib {

namespace highway = hwy::HWY_NAMESPACE;
/** Highway tag for full double register */
inline constexpr highway::ScalableTag<double> tag_double;
/** Highway tag for full long register */
inline constexpr highway::ScalableTag<int64_t> tag_long;
/**
 * Number of double values in a full register. Not constexpr on all architectures, so cannot be assumed, but, when
 * possible (e.g. on x86), we provide this constexpr so compilers may make optimizations.
 */
HWY_LANES_CONSTEXPR inline size_t _vecLengthDouble{highway::Lanes(tag_double)};
/**
 * Upper bound on the number of double values in a full register. Always constexpr.
 */
inline constexpr size_t _maxVecLengthDouble{highway::MaxLanes(tag_double)};
/** Type for a Double vector register */
using VectorDouble = decltype(highway::Zero(tag_double));
/** Type for a Long vector register */
using VectorLong = decltype(highway::Zero(tag_long));
/** Highway tag for a half-filled double register */
inline constexpr highway::Half<highway::DFromV<VectorDouble>> tag_double_half;
/** Type for a Double Mask */
using MaskDouble = decltype(highway::FirstN(tag_double, 1));
/** Type for a Long Mask */
using MaskLong = decltype(highway::FirstN(tag_long, 2));

namespace hwy_helpers {
/**
 * Selects the tuned fixed-width halving/interleave horizontal reduction path over the generic highway::ReduceSum
 * fallback.
 *
 * Disabled on:
 * - Scalable targets (HWY_HAVE_SCALABLE != 0, e.g. Arm SVE/SVE2, RISC-V RVV), where runtime vector lengths can exceed
 *   8 lanes and runtime 1-lane half-vectors cannot use MaxLanes()-only compile-time branching.
 * - Fixed SVE targets (HWY_SVE_256, HWY_SVE2_128), where DFromV<svfloat64_t> resolves back to ScalableTag<double>
 *   rather than a 2-lane vector type, LowerHalf is a no-op, and ReduceSum lowers directly to predicated svaddv.
 * - Targets with MaxLanes == 1 or > 8 double lanes.
 *
 * @tparam D Highway Simd tag for the vector to reduce.
 * @tparam forceGenericFallback Test hook to force the generic ReduceSum fallback on fixed-width targets.
 */
template <class D, bool forceGenericFallback = false>
constexpr bool useFixedHalvingReduction() {
  if constexpr (forceGenericFallback) {
    return false;
  } else {
    constexpr size_t maxLanes = highway::MaxLanes(D{});
    constexpr bool isSveTarget = (HWY_TARGET & HWY_ALL_SVE) != 0;
    return !HWY_HAVE_SCALABLE && !isSveTarget && (maxLanes == 2 || maxLanes == 4 || maxLanes == 8);
  }
}

/**
 * Recursively folds a fixed 4- or 8-lane double vector down to 2 lanes by adding its lower and upper halves.
 */
template <class D>
HWY_INLINE auto foldTo2Lanes(D /*d*/, const highway::VFromD<D> &v) {
  constexpr size_t maxLanes = highway::MaxLanes(D{});
  static_assert(maxLanes <= 8 && !HWY_HAVE_SCALABLE && ((HWY_TARGET & HWY_ALL_SVE) == 0),
                "foldTo2Lanes is only valid for fixed non-SVE vectors with at most 8 double lanes.");
  if constexpr (maxLanes == 8) {
    const highway::Half<D> d4;
    const auto v4 = highway::Add(highway::LowerHalf(d4, v), highway::UpperHalf(d4, v));
    const highway::Half<decltype(d4)> d2;
    return highway::Add(highway::LowerHalf(d2, v4), highway::UpperHalf(d2, v4));
  } else if constexpr (maxLanes == 4) {
    const highway::Half<D> d2;
    return highway::Add(highway::LowerHalf(d2, v), highway::UpperHalf(d2, v));
  } else {
    return v;
  }
}

/**
 * Reduces a single vector to a scalar sum, using tree halving on fixed 4/8-lane AVX targets and highway::ReduceSum
 * elsewhere (scalable SVE/RVV, fixed SVE, 128-bit, >8-lane, and scalar).
 */
template <bool forceGenericFallback = false, class D>
HWY_INLINE double reduceSumHalving([[maybe_unused]] D d, const highway::VFromD<D> &v) {
  if constexpr (useFixedHalvingReduction<D, forceGenericFallback>() && highway::MaxLanes(D{}) > 2) {
    const auto v2 = foldTo2Lanes(d, v);
    const highway::FixedTag<double, 2> d2;
    return highway::ReduceSum(d2, v2);
  } else if constexpr (highway::MaxLanes(D{}) == 1) {
    return highway::GetLane(v);
  } else {
    if (highway::Lanes(d) == 1) {
      return highway::GetLane(v);
    }
    return highway::ReduceSum(d, v);
  }
}

/**
 * Shared implementation for reduceAdd2_InterleaveFirst and reduceSub2_InterleaveFirst.
 */
template <bool subtract, bool forceGenericFallback = false, class D>
HWY_INLINE void reduce2_InterleaveFirst([[maybe_unused]] D d, const highway::VFromD<D> &va,
                                        const highway::VFromD<D> &vb, double &outA, double &outB) {
  double sumA{};
  double sumB{};
  if constexpr (useFixedHalvingReduction<D, forceGenericFallback>()) {
    const auto ab = highway::Add(highway::InterleaveLower(d, va, vb), highway::InterleaveUpper(d, va, vb));
    const auto ab2 = foldTo2Lanes(d, ab);
    sumA = highway::GetLane(ab2);
    sumB = highway::ExtractLane(ab2, 1);
  } else if constexpr (highway::MaxLanes(D{}) == 1) {
    sumA = highway::GetLane(va);
    sumB = highway::GetLane(vb);
  } else {
    if (highway::Lanes(d) == 1) {
      sumA = highway::GetLane(va);
      sumB = highway::GetLane(vb);
    } else {
      sumA = highway::ReduceSum(d, va);
      sumB = highway::ReduceSum(d, vb);
    }
  }
  if constexpr (subtract) {
    outA -= sumA;
    outB -= sumB;
  } else {
    outA += sumA;
    outB += sumB;
  }
}

/**
 * Reduces two vectors simultaneously and adds the results to output scalar accumulators.
 */
template <bool forceGenericFallback = false, class D>
HWY_INLINE void reduceAdd2_InterleaveFirst(D d, const highway::VFromD<D> &va, const highway::VFromD<D> &vb,
                                           double &outA, double &outB) {
  reduce2_InterleaveFirst<false, forceGenericFallback>(d, va, vb, outA, outB);
}

/**
 * Reduces two vectors simultaneously and subtracts the results from output scalar accumulators.
 */
template <bool forceGenericFallback = false, class D>
HWY_INLINE void reduceSub2_InterleaveFirst(D d, const highway::VFromD<D> &va, const highway::VFromD<D> &vb,
                                           double &outA, double &outB) {
  reduce2_InterleaveFirst<true, forceGenericFallback>(d, va, vb, outA, outB);
}

/**
 * Shared implementation for reduceAdd3_InterleaveFirst and reduceSub3_InterleaveFirst.
 */
template <bool subtract, bool forceGenericFallback = false, class D>
HWY_INLINE void reduce3_InterleaveFirst(D d, const highway::VFromD<D> &fx, const highway::VFromD<D> &fy,
                                        const highway::VFromD<D> &fz, double &outX, double &outY, double &outZ) {
  reduce2_InterleaveFirst<subtract, forceGenericFallback>(d, fx, fy, outX, outY);
  const double sumZ = reduceSumHalving<forceGenericFallback>(d, fz);
  if constexpr (subtract) {
    outZ -= sumZ;
  } else {
    outZ += sumZ;
  }
}

/**
 * Reduces 3 vectors (x, y, z) into scalar accumulators by interleaving (x, y) and halving/reducing z.
 */
template <bool forceGenericFallback = false, class D>
HWY_INLINE void reduceAdd3_InterleaveFirst(D d, const highway::VFromD<D> &fx, const highway::VFromD<D> &fy,
                                           const highway::VFromD<D> &fz, double &outX, double &outY, double &outZ) {
  reduce3_InterleaveFirst<false, forceGenericFallback>(d, fx, fy, fz, outX, outY, outZ);
}

/**
 * Reduces 3 vectors (x, y, z) and subtracts the results from output scalar accumulators.
 */
template <bool forceGenericFallback = false, class D>
HWY_INLINE void reduceSub3_InterleaveFirst(D d, const highway::VFromD<D> &fx, const highway::VFromD<D> &fy,
                                           const highway::VFromD<D> &fz, double &outX, double &outY, double &outZ) {
  reduce3_InterleaveFirst<true, forceGenericFallback>(d, fx, fy, fz, outX, outY, outZ);
}

/**
 * Reduces 4 vectors into scalar accumulators using two parallel interleaved pairs (or fallback ReduceSum).
 */
template <bool forceGenericFallback = false, class D>
HWY_INLINE void reduceAdd4_Interleaved(D d, const highway::VFromD<D> &v0, const highway::VFromD<D> &v1,
                                       const highway::VFromD<D> &v2, const highway::VFromD<D> &v3, double &out0,
                                       double &out1, double &out2, double &out3) {
  reduce2_InterleaveFirst<false, forceGenericFallback>(d, v0, v1, out0, out1);
  reduce2_InterleaveFirst<false, forceGenericFallback>(d, v2, v3, out2, out3);
}

/**
 * Shared implementation for reduceAdd6Half_Vec2Store and reduceSub6Half_Vec2Store.
 *
 * Updates strictly 2 contiguous scalar positions per coordinate (ptr[0] from lf*, ptr[1] from uf*). On fixed
 * 4/8-lane AVX targets this uses 2-lane vector memory ops via FixedTag<double, 2>; on scalable/SVE/scalar targets it
 * accumulates via GetLane or highway::ReduceSum and scalar updates so vector-length mismatches or SVE DFromV typing
 * never over-read or overwrite adjacent elements.
 */
template <bool subtract, bool forceGenericFallback = false, class DHalf>
HWY_INLINE void reduce6Half_Vec2Store([[maybe_unused]] DHalf dh, const highway::VFromD<DHalf> &lfx,
                                      const highway::VFromD<DHalf> &lfy, const highway::VFromD<DHalf> &lfz,
                                      const highway::VFromD<DHalf> &ufx, const highway::VFromD<DHalf> &ufy,
                                      const highway::VFromD<DHalf> &ufz, double *const __restrict fxPtr,
                                      double *const __restrict fyPtr, double *const __restrict fzPtr) {
  if constexpr (useFixedHalvingReduction<DHalf, forceGenericFallback>()) {
    const auto fxLU = highway::Add(highway::InterleaveLower(dh, lfx, ufx), highway::InterleaveUpper(dh, lfx, ufx));
    const auto fyLU = highway::Add(highway::InterleaveLower(dh, lfy, ufy), highway::InterleaveUpper(dh, lfy, ufy));
    const auto fzLU = highway::Add(highway::InterleaveLower(dh, lfz, ufz), highway::InterleaveUpper(dh, lfz, ufz));
    const auto fx2 = foldTo2Lanes(dh, fxLU);
    const auto fy2 = foldTo2Lanes(dh, fyLU);
    const auto fz2 = foldTo2Lanes(dh, fzLU);
    const highway::FixedTag<double, 2> d2;

    if constexpr (subtract) {
      highway::StoreU(highway::Sub(highway::LoadU(d2, fxPtr), fx2), d2, fxPtr);
      highway::StoreU(highway::Sub(highway::LoadU(d2, fyPtr), fy2), d2, fyPtr);
      highway::StoreU(highway::Sub(highway::LoadU(d2, fzPtr), fz2), d2, fzPtr);
    } else {
      highway::StoreU(highway::Add(highway::LoadU(d2, fxPtr), fx2), d2, fxPtr);
      highway::StoreU(highway::Add(highway::LoadU(d2, fyPtr), fy2), d2, fyPtr);
      highway::StoreU(highway::Add(highway::LoadU(d2, fzPtr), fz2), d2, fzPtr);
    }
  } else {
    double lfxSum{}, lfySum{}, lfzSum{}, ufxSum{}, ufySum{}, ufzSum{};
    auto assignGetLane = [&]() {
      lfxSum = highway::GetLane(lfx);
      lfySum = highway::GetLane(lfy);
      lfzSum = highway::GetLane(lfz);
      ufxSum = highway::GetLane(ufx);
      ufySum = highway::GetLane(ufy);
      ufzSum = highway::GetLane(ufz);
    };
    if constexpr (highway::MaxLanes(DHalf{}) == 1) {
      assignGetLane();
    } else if (highway::Lanes(dh) == 1) {
      assignGetLane();
    } else {
      lfxSum = highway::ReduceSum(dh, lfx);
      lfySum = highway::ReduceSum(dh, lfy);
      lfzSum = highway::ReduceSum(dh, lfz);
      ufxSum = highway::ReduceSum(dh, ufx);
      ufySum = highway::ReduceSum(dh, ufy);
      ufzSum = highway::ReduceSum(dh, ufz);
    }
    if constexpr (subtract) {
      fxPtr[0] -= lfxSum;
      fyPtr[0] -= lfySum;
      fzPtr[0] -= lfzSum;
      fxPtr[1] -= ufxSum;
      fyPtr[1] -= ufySum;
      fzPtr[1] -= ufzSum;
    } else {
      fxPtr[0] += lfxSum;
      fyPtr[0] += lfySum;
      fzPtr[0] += lfzSum;
      fxPtr[1] += ufxSum;
      fyPtr[1] += ufySum;
      fzPtr[1] += ufzSum;
    }
  }
}

/**
 * Reduces 6 half-vectors (lower and upper halves of 3 force components) and adds them to 2 contiguous particles.
 */
template <bool forceGenericFallback = false, class DHalf>
HWY_INLINE void reduceAdd6Half_Vec2Store(DHalf dh, const highway::VFromD<DHalf> &lfx, const highway::VFromD<DHalf> &lfy,
                                         const highway::VFromD<DHalf> &lfz, const highway::VFromD<DHalf> &ufx,
                                         const highway::VFromD<DHalf> &ufy, const highway::VFromD<DHalf> &ufz,
                                         double *const __restrict fxPtr, double *const __restrict fyPtr,
                                         double *const __restrict fzPtr) {
  reduce6Half_Vec2Store<false, forceGenericFallback>(dh, lfx, lfy, lfz, ufx, ufy, ufz, fxPtr, fyPtr, fzPtr);
}

/**
 * Reduces 6 half-vectors (lower and upper halves of 3 force components) and subtracts them from 2 contiguous
 * particles.
 */
template <bool forceGenericFallback = false, class DHalf>
HWY_INLINE void reduceSub6Half_Vec2Store(DHalf dh, const highway::VFromD<DHalf> &lfx, const highway::VFromD<DHalf> &lfy,
                                         const highway::VFromD<DHalf> &lfz, const highway::VFromD<DHalf> &ufx,
                                         const highway::VFromD<DHalf> &ufy, const highway::VFromD<DHalf> &ufz,
                                         double *const __restrict fxPtr, double *const __restrict fyPtr,
                                         double *const __restrict fzPtr) {
  reduce6Half_Vec2Store<true, forceGenericFallback>(dh, lfx, lfy, lfz, ufx, ufy, ufz, fxPtr, fyPtr, fzPtr);
}

}  // namespace hwy_helpers

}  // namespace mdLib

/**
 * @file FunctorBenchmarkTraits.h
 * @author H. Meyran
 */

#pragma once

namespace autopas {

/**
 * Trait struct declaring which startup benchmarks a functor supports.
 * Default: no benchmarks supported. Specialize for functors that opt in.
 * Add new boolean fields here as new benchmark types are introduced.
 * @tparam Functor_T The functor type to query.
 */
template <typename Functor_T>
struct FunctorBenchmarkTraits {
  /**
   * Whether the functor supports PatternBenchmark (vectorization pattern lookup table selection via SoAFunctorPair
   * micro-benchmark). A functor that sets this to true must provide a real setVecPattern() implementation and a
   * SoAFunctorPair() that dispatches on the currently set pattern.
   */
  static constexpr bool supportsPatternBenchmark = false;
};

}  // namespace autopas

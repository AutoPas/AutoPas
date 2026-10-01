/**
 * @file HashCombineTest.cpp
 * @date 14/09/2026
 * @author S. J. Newcome
 */

#include <gtest/gtest.h>

#include <cstddef>
#include <functional>
#include <unordered_set>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/utils/HashCombine.h"

using namespace autopas;

namespace {

/**
 * Minimal stand-in for an AutoPas Option: implicitly convertible to std::size_t but with no std::hash specialization.
 *
 * This exercises the static_cast fallback path of hashValue(). It is also used for the reference value test below,
 * because its hashValue() is exactly the wrapped number, whereas std::hash of a built-in integer is only required to
 * be consistent, not to have any particular value, and so differs between standard library implementations.
 */
struct SizeTLike {
  /**
   * The wrapped value.
   */
  std::size_t value;

  /**
   * Implicit conversion, mirroring Option::operator Value().
   * @return The wrapped value.
   */
  constexpr operator std::size_t() const { return value; }
};

}  // namespace

/**
 * Checks hashCombine() against hard coded reference values, to pin down the mixing constants and the order of
 * operations. Without this, hashCombine() could be silently altered (e.g. a mistyped constant or a wrong shift width)
 * and every other test here would still pass, because they only compare hashes against each other.
 *
 * The reference values were not taken from this implementation but computed independently from the definition of
 * boost::hash_combine (Boost >= 1.81), i.e. by repeatedly applying, in 64 bit unsigned arithmetic:
 *
 *   seed = mix(seed + 0x9e3779b9 + value), starting from seed = 0
 *   mix(x): x ^= x >> 32; x *= 0xe9846af9b1a615d; x ^= x >> 32; x *= 0xe9846af9b1a615d; x ^= x >> 28
 *
 * Recomputing them from this formula is enough to regenerate the expectations if the algorithm is ever changed.
 */
TEST(HashCombineTest, testReferenceValues) {
  EXPECT_EQ(utils::hashCombine(SizeTLike{42}), 0x393c360f4e323eaeULL);
  EXPECT_EQ(utils::hashCombine(SizeTLike{1}, SizeTLike{2}), 0x30b3fc98529bf99eULL);
  EXPECT_EQ(utils::hashCombine(SizeTLike{2}, SizeTLike{1}), 0x31854bc10639eee4ULL);
  EXPECT_EQ(utils::hashCombine(SizeTLike{7}, SizeTLike{7}), 0x65c22b14c97b54bfULL);
}

/**
 * hashValue() should defer to std::hash where it exists and fall back to a static_cast otherwise.
 */
TEST(HashCombineTest, testHashValueDispatch) {
  EXPECT_EQ(utils::hashValue(123), std::hash<int>{}(123));
  EXPECT_EQ(utils::hashValue(1.5), std::hash<double>{}(1.5));
  EXPECT_EQ(utils::hashValue(SizeTLike{5}), 5ul);
  EXPECT_EQ(utils::hashValue(ContainerOption{ContainerOption::linkedCells}),
            static_cast<std::size_t>(ContainerOption::linkedCells));
}

/**
 * Equal inputs must always produce equal hashes.
 */
TEST(HashCombineTest, testDeterministic) {
  EXPECT_EQ(utils::hashCombine(1, 2.5, ContainerOption{ContainerOption::linkedCells}),
            utils::hashCombine(1, 2.5, ContainerOption{ContainerOption::linkedCells}));
  EXPECT_EQ(utils::hashCombine(), utils::hashCombine());
  EXPECT_EQ(utils::hashCombine(SizeTLike{13}, SizeTLike{13}), utils::hashCombine(SizeTLike{13}, SizeTLike{13}));
  EXPECT_EQ(utils::hashCombine(Newton3Option{Newton3Option::enabled}, -0.5, SizeTLike{0}),
            utils::hashCombine(Newton3Option{Newton3Option::enabled}, -0.5, SizeTLike{0}));
}

/**
 * Unlike XOR-ing or summing hashes, hashCombine() must be order dependent. Otherwise configurations that merely
 * permute the same option values (e.g. a container and a traversal that happen to have swapped indices) would share a
 * hash.
 */
TEST(HashCombineTest, testOrderDependent) {
  EXPECT_NE(utils::hashCombine(1, 2), utils::hashCombine(2, 1));
  EXPECT_NE(utils::hashCombine(1, 2, 3), utils::hashCombine(3, 2, 1));
  EXPECT_NE(utils::hashCombine(SizeTLike{4}, SizeTLike{9}), utils::hashCombine(SizeTLike{9}, SizeTLike{4}));
  EXPECT_NE(utils::hashCombine(ContainerOption{ContainerOption::verletLists}, Newton3Option{Newton3Option::enabled}),
            utils::hashCombine(Newton3Option{Newton3Option::enabled}, ContainerOption{ContainerOption::verletLists}));
}

/**
 * Combining a value twice must differ from combining it once and from combining it not at all, which is a failure
 * point in other hash combine functions (e.g. XOR-ing hashes)
 */
TEST(HashCombineTest, testRepeatedValuesDoNotCancel) {
  EXPECT_NE(utils::hashCombine(7, 7), 0ul);
  EXPECT_NE(utils::hashCombine(7, 7), utils::hashCombine());
  EXPECT_NE(utils::hashCombine(7, 7), utils::hashCombine(7));
  EXPECT_NE(utils::hashCombine(7, 7, 3), utils::hashCombine(3));
}

/**
 * The number of combined values must matter, even if the additional ones are zero.
 */
TEST(HashCombineTest, testArityMatters) {
  EXPECT_NE(utils::hashCombine(0), utils::hashCombine(0, 0));
  EXPECT_NE(utils::hashCombine(0, 0), utils::hashCombine(0, 0, 0));
}

/**
 * Every argument must influence the result, so that no member of a hashed struct is silently ignored.
 */
TEST(HashCombineTest, testEveryArgumentContributes) {
  const auto reference = utils::hashCombine(1, 2, 3, 4);
  EXPECT_NE(utils::hashCombine(9, 2, 3, 4), reference);
  EXPECT_NE(utils::hashCombine(1, 9, 3, 4), reference);
  EXPECT_NE(utils::hashCombine(1, 2, 9, 4), reference);
  EXPECT_NE(utils::hashCombine(1, 2, 3, 9), reference);
}

/**
 * Doubles must be hashed as doubles and not truncated to integers by an accidental static_cast.
 */
TEST(HashCombineTest, testDoublesAreNotTruncated) {
  EXPECT_NE(utils::hashCombine(1.0), utils::hashCombine(1.5));
  EXPECT_NE(utils::hashCombine(1.0, 2.0), utils::hashCombine(1.2, 2.3));
}
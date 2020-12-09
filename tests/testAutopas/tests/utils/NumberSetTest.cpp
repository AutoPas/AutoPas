/**
 * @file NumberSetTest.cpp
 * @author Jan Nguyen
 * @date 09.06.19
 */

#include <gtest/gtest.h>

#include <map>

#include "autopas/utils/NumberInterval.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/NumberSetFinite.h"

namespace NumberSetTest {

using namespace autopas;

TEST(NumberSetTest, testFiniteSet) {
  std::array<double, 3> a = {1., 2., 3.};
  NumberSetFinite<double> fSet({1., 2., 3.});
  NumberSet<double> &set = fSet;

  // check min and max
  ASSERT_EQ(set.getMin(), 1.);
  ASSERT_EQ(set.getMax(), 3.);

  // set should be finite by definition
  ASSERT_TRUE(set.isFinite());

  // set should keep all initialized values
  std::set<double> values = set.getAll();
  ASSERT_EQ(values.size(), a.size());
  ASSERT_EQ(set.size(), a.size());
  auto it = values.begin();
  for (unsigned indexFirst = 0; indexFirst < a.size(); ++indexFirst, ++it) {
    ASSERT_EQ(*it, a[indexFirst]);
  }
}

TEST(NumberSetTest, testFiniteSetSample) {
  Random rng(32);
  unsigned numSamples = 32;

  std::array<double, 3> a = {1., 2., 3.};
  NumberSetFinite<double> fSet({1., 2., 3.});
  NumberSet<double> &set = fSet;

  std::map<double, unsigned> occurences;
  for (auto val : a) {
    occurences[val] = 0;
  }
  std::vector<double> samples = set.uniformSample(numSamples, rng);
  ASSERT_EQ(samples.size(), numSamples);
  for (auto sample : samples) {
    // sample should be in the set
    ASSERT_EQ(occurences.count(sample), 1);

    // count occurences of the samples
    ++occurences[sample];
  }

  // get the min and max occurences
  unsigned minOccur = numSamples;
  unsigned maxOccur = 0;
  for (auto val : a) {
    if (occurences[val] < minOccur) {
      minOccur = occurences[val];
    }

    if (occurences[val] > maxOccur) {
      maxOccur = occurences[val];
    }
  }

  // occurences should be equally distributed
  ASSERT_LE(maxOccur - minOccur, 1);
}

TEST(NumberSetTest, testInterval) {
  double min = -3.;
  double max = 2.;
  NumberInterval<double> iSet(min, max);
  NumberSet<double> &set = iSet;

  // set should be infinite
  ASSERT_FALSE(set.isFinite());

  // cant get all objects or size of an infinite set
  EXPECT_THROW(set.getAll(), autopas::utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(set.size(), autopas::utils::ExceptionHandler::AutoPasException);
}

TEST(NumberSetTest, testIntervalFinite) {
  double min = 2.;
  double max = 2.;
  NumberInterval<double> iSet(min, max);
  NumberSet<double> &set = iSet;

  // set should be finite
  ASSERT_TRUE(set.isFinite());

  // set should contain exactly one value
  ASSERT_EQ(set.size(), 1);
  std::set<double> values = set.getAll();
  ASSERT_EQ(values.size(), 1);
  ASSERT_EQ(*values.begin(), max);
}

TEST(NumberSetTest, testIntervalSample) {
  Random rng(32);
  unsigned numSamples = 32;

  double min = -3.;
  double max = 2.;
  NumberInterval<double> iSet(min, max);
  NumberSet<double> &set = iSet;

  std::vector<double> samples = set.uniformSample(numSamples, rng);
  ASSERT_EQ(samples.size(), numSamples);
  for (auto sample : samples) {
    // sample should be in the set
    ASSERT_GE(sample, min);
    ASSERT_LE(sample, max);
  }
}

TEST(NumberSetTest, testClone) {
  double min = -3.;
  double max = 2.;
  NumberInterval<double> iSet(min, max);
  NumberSet<double> &set = iSet;

  auto clone = set.clone();

  // clone should be newly allocated
  ASSERT_NE(clone.get(), &set);

  // min and max should be same
  ASSERT_EQ(clone->getMin(), set.getMin());
  ASSERT_EQ(clone->getMax(), set.getMax());
}

TEST(NumberSetTest, testStringRepresentation) {
  // Empty set
  NumberSetFinite<double> emptyfSet;
  EXPECT_EQ("[]", emptyfSet.to_string());

  // Filled set
  NumberSetFinite<double> fSet({1., 2., 3.});
  EXPECT_EQ("[1, 2, 3]", fSet.to_string());
}

}  // end namespace NumberSetTest

/**
 * @file StringUtilsTest.cpp
 * @author F. Gratl
 * @date 1/15/19
 */

#include "StringUtilsTest.h"

// ToDo this needs to be generalized now that these parsers contain any type

TEST(StringUtilsTest, parseDoublesTest) {
  auto parsedOptions = autopas::utils::StringUtils::parseNumbers<double>("1.,1.5, 2,3.00,2e1");

  EXPECT_THAT(parsedOptions, ::testing::ElementsAreArray({1., 1.5, 2., 3., 20.}));
}

TEST(StringUtilsTest, parseNumberSetFiniteTest) {
  auto numberSet = autopas::utils::StringUtils::parseNumberSet<double>("1.,1.5, 2,3.00,2e1, 1.2e-2");
  auto numberSetFinite = dynamic_cast<autopas::NumberSetFinite<double> *>(numberSet.get());
  ASSERT_NE(numberSetFinite, nullptr);  // if this is null numberSet was parsed to a NumberInterval
  EXPECT_THAT(numberSetFinite->getAll(), ::testing::UnorderedElementsAre(1., 1.5, 2., 3., 20., 0.012));
}

/**
 * Checks that parseNumber() only accepts strings that are entirely a number of the requested type.
 */
TEST(StringUtilsTest, parseNumberTest) {
  using autopas::utils::StringUtils::parseNumber;

  EXPECT_EQ(parseNumber<size_t>("42"), 42);
  EXPECT_EQ(parseNumber<double>("1.5"), 1.5);
  EXPECT_EQ(parseNumber<double>("2e1"), 20.);

  EXPECT_ANY_THROW(parseNumber<size_t>("1.5"));
  EXPECT_ANY_THROW(parseNumber<double>("1.5x"));
  EXPECT_ANY_THROW(parseNumber<double>(""));
  EXPECT_ANY_THROW(parseNumber<float>("1e100"));
}

TEST(StringUtilsTest, parseNumberIntervalTest) {
  auto numberSet = autopas::utils::StringUtils::parseNumberSet<double>("1.-2e1");
  auto numberInterval = dynamic_cast<autopas::NumberInterval<double> *>(numberSet.get());
  ASSERT_NE(numberInterval, nullptr);  // if this is null numberSet was parsed to a NumberSetFinite
  EXPECT_EQ(numberInterval->getMin(), 1.);
  EXPECT_EQ(numberInterval->getMax(), 2e1);
}

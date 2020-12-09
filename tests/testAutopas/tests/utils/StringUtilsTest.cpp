/**
 * @file StringUtilsTest.cpp
 * @author F. Gratl
 * @date 1/15/19
 */

#include "StringUtilsTest.h"

namespace StringUtilsTest {

TEST(StringUtilsTest, parseDoublesTest) {
  auto parsedOptions = autopas::utils::StringUtils::parseDoubles("1.,1.5, 2,3.00,2e1");

  EXPECT_THAT(parsedOptions, ::testing::ElementsAreArray({1., 1.5, 2., 3., 20.}));
}

TEST(StringUtilsTest, parseNumberSetFiniteTest) {
  auto numberSet = autopas::utils::StringUtils::parseNumberSet("1.,1.5, 2,3.00,2e1, 1.2e-2");
  auto numberSetFinite = dynamic_cast<autopas::NumberSetFinite<double> *>(numberSet.get());
  ASSERT_NE(numberSetFinite, nullptr);  // if this is null numberSet was parsed to a NumberInterval
  EXPECT_THAT(numberSetFinite->getAll(), ::testing::UnorderedElementsAre(1., 1.5, 2., 3., 20., 0.012));
}

TEST(StringUtilsTest, parseNumberIntervalTest) {
  auto numberSet = autopas::utils::StringUtils::parseNumberSet("1.-2e1");
  auto numberInterval = dynamic_cast<autopas::NumberInterval<double> *>(numberSet.get());
  ASSERT_NE(numberInterval, nullptr);  // if this is null numberSet was parsed to a NumberSetFinite
  EXPECT_EQ(numberInterval->getMin(), 1.);
  EXPECT_EQ(numberInterval->getMax(), 2e1);
}

}  // end namespace StringUtilsTest

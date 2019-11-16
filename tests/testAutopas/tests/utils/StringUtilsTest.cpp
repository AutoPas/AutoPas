/**
 * @file StringUtilsTest.cpp
 * @author F. Gratl
 * @date 1/15/19
 */

#include "StringUtilsTest.h"

TEST(StringUtilsTest, parseDoublesTest) {
  auto parsedOptions = autopas::utils::StringUtils::parseDoubles("1.,1.5, 2,3.00,2e1");

  EXPECT_THAT(parsedOptions, ::testing::ElementsAreArray({1., 1.5, 2., 3., 20.}));
}

TEST(StringUtilsTest, parseNumberSetTest) {
  EXPECT_EQ(autopas::utils::StringUtils::parseNumberSet("1.,1.5, 2,3.00,2e1")->getAll(),
            std::set<double>({1., 1.5, 2., 3., 20.}));

  auto numberSet = autopas::utils::StringUtils::parseNumberSet("[1.,2e1]");
  auto *numberInterval = dynamic_cast<autopas::NumberInterval<double> *>(numberSet.get());
  ASSERT_NE(numberInterval, nullptr);
  EXPECT_EQ(numberInterval->getMin(), 1.);
  EXPECT_EQ(numberInterval->getMax(), 2e1);
}

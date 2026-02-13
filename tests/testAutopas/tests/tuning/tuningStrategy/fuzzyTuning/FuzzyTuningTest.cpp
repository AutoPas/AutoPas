/**
 * @file FuzzyTuningTest.cpp
 * @author Manuel Lerchner
 * @date 30.06.2024
 */

// Disable these tests if fuzzy tuning is not enabled
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
#include "FuzzyTuningTest.h"

#include <filesystem>
#include <fstream>

#include "autopas/tuning/tuningStrategy/fuzzyTuning/FuzzyTuning.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzySetFactory.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/LinguisticVariable.h"

using namespace autopas;
using namespace autopas::FuzzyLogic;

/**
 * Tests whether the triangle fuzzy set has the correct membership values at certain key points.
 */
TEST(FuzzyTuningTest, testTriangleFuzzySet) {
  // creates the (continuous) crisp set over which the fuzzy set is defined. It has the name "x" and the range [-5, 30]
  std::shared_ptr<CrispSet> crispSet = std::make_shared<CrispSet>("x", std::pair(-5, 30));

  // creates a fuzzy set with a triangle membership function which is 0 for values smaller than 0, increases
  // linearly from 0 to 1 for values between 0 and 5 and decreases linearly from 1 to 0 for values between 5 and 10
  auto triangle = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 5, 10});
  triangle->setCrispSet(crispSet);

  const auto epsilon = 1e-10;

  // test the membership function
  const std::map<std::string, double> xMinus2 = {{"x", -2}};
  EXPECT_NEAR(triangle->evaluate_membership(xMinus2), 0, epsilon);

  const std::map<std::string, double> x0 = {{"x", 0}};
  EXPECT_NEAR(triangle->evaluate_membership(x0), 0, epsilon);

  const std::map<std::string, double> x2 = {{"x", 2}};

  // as x = 2 is in the increasing part of the triangle, and the triangle increases linearly from 0 to 1 between 0 and 5
  // we can calculate the membership value at x = 2 as follows:
  // y = delta_y / delta_x * (x - 0) + 0 = (1.0 - 0) / (5 - 0) * (2 - 0) + 0 = 0.4
  EXPECT_NEAR(triangle->evaluate_membership(x2), 0.4, epsilon);

  const std::map<std::string, double> x5 = {{"x", 5}};
  EXPECT_NEAR(triangle->evaluate_membership(x5), 1, epsilon);

  const std::map<std::string, double> x7 = {{"x", 7}};

  // as x = 7 is in the decreasing part of the triangle, and the triangle decreases linearly from 1 to 0 between 5 and
  // 10 we can calculate the membership value at x = 7 as:
  // y = delta_y / delta_x * (x - 5) +1  = (0 - 1) / (10 - 5) * (7 - 5) + 1 = -0.4 + 1 = 0.6
  EXPECT_NEAR(triangle->evaluate_membership(x7), 0.6, epsilon);

  const std::map<std::string, double> x10 = {{"x", 10}};
  EXPECT_NEAR(triangle->evaluate_membership(x10), 0, epsilon);

  const std::map<std::string, double> x12 = {{"x", 12}};
  EXPECT_NEAR(triangle->evaluate_membership(x12), 0, epsilon);
}

/**
 * Tests whether the gaussian fuzzy set takes the correct values at certain key points.
 */
TEST(FuzzyTuningTest, testGaussianFuzzySet) {
  // creates the (continuous) crisp set over which the fuzzy set is defined. It has the name "x" and the range [-5, 30]
  std::shared_ptr<CrispSet> crispSet = std::make_shared<CrispSet>("x", std::pair(-5, 30));

  // creates a fuzzy set with a gaussian membership function with mean 1 and standard deviation 2
  auto gaussian = FuzzySetFactory::makeFuzzySet("low", "Gaussian", {1, 2});
  gaussian->setCrispSet(crispSet);

  const auto epsilon = 1e-10;

  // test the membership function
  const std::map<std::string, double> xMinus10 = {{"x", -100}};
  EXPECT_NEAR(gaussian->evaluate_membership(xMinus10), 0, epsilon);

  const std::map<std::string, double> x0 = {{"x", 0}};

  // the gaussian function is defined as f(x) = exp(-0.5 * ((x - mean) / std_dev)^2)
  // for the given parameters mean = 1 and std_dev = 2 the membership value at x = 0 should be:
  // f(0) = exp(-0.5 * ((0 - 1) / 2)^2) = exp(-0.5 * (-0.5)^2) = exp(-0.5 * 0.25) = exp(-0.125) ≈ 0.88249690258
  EXPECT_NEAR(gaussian->evaluate_membership(x0), 0.88249690258, epsilon);

  const std::map<std::string, double> x1 = {{"x", 1}};
  EXPECT_NEAR(gaussian->evaluate_membership(x1), 1, epsilon);

  const std::map<std::string, double> x2 = {{"x", 2}};

  // f(2) = exp(-0.5 * ((2 - 1) / 2)^2) = exp(-0.5 * 0.5^2) = exp(-0.5 * 0.25) = exp(-0.125) ≈ 0.88249690258
  EXPECT_NEAR(gaussian->evaluate_membership(x2), 0.88249690258, epsilon);

  const std::map<std::string, double> x10 = {{"x", 10}};
  // f(10) = exp(-0.5 * ((10 - 1) / 2)^2) = exp(-0.5 * 4.5^2) = exp(-0.5 * 20.25) = exp(-10.125) ≈ 0.000040065297
  EXPECT_NEAR(gaussian->evaluate_membership(x10), 0.000040065297, epsilon);
}

/**
 * Tests whether the sigmoid finite fuzzy set (increasing) takes the correct values at certain key points.
 * Since the sigmoid finite fuzzy set is setup with the left<mid<right its membership function should be increasing.
 */
TEST(FuzzyTuningTest, testSigmoidFiniteFuzzySetIncreasing) {
  // creates the (continuous) crisp set over which the fuzzy set is defined. It has the name "x" and the range [-5, 30]
  std::shared_ptr<CrispSet> crispSet = std::make_shared<CrispSet>("x", std::pair(-5, 30));

  // creates a fuzzy set with a sigmoid finite membership function inbetween 0 and 10 with a center at 5
  auto sigmoidFinite = FuzzySetFactory::makeFuzzySet("low", "SigmoidFinite", {0, 5, 10});
  sigmoidFinite->setCrispSet(crispSet);

  const auto epsilon = 1e-10;

  // test the membership function
  EXPECT_NEAR(sigmoidFinite->evaluate_membership({{"x", -10}}), 0, epsilon);

  EXPECT_NEAR(sigmoidFinite->evaluate_membership({{"x", 0}}), 0, epsilon);

  EXPECT_NEAR(sigmoidFinite->evaluate_membership({{"x", 5}}), 0.5, epsilon);

  EXPECT_NEAR(sigmoidFinite->evaluate_membership({{"x", 10}}), 1, epsilon);
}

/**
 * Tests whether the sigmoid fuzzy set (decreasing) takes the correct values at certain key points.
 * Since the sigmoid finite fuzzy set is setup with the left>mid>right its membership function should be decreasing.
 */
TEST(FuzzyTuningTest, testSigmoidFiniteFuzzySetDecreasing) {
  // creates the (continuous) crisp set over which the fuzzy set is defined. It has the name "x" and the range [-5, 30]
  std::shared_ptr<CrispSet> crispSet = std::make_shared<CrispSet>("x", std::pair(-5, 30));

  // creates a fuzzy set with a sigmoid finite membership function with a center at 5 compared to the one before the
  // left and right border are switched
  auto sigmoidFinite = FuzzySetFactory::makeFuzzySet("low", "SigmoidFinite", {10, 5, 0});
  sigmoidFinite->setCrispSet(crispSet);

  const auto epsilon = 1e-10;

  // test the membership function
  EXPECT_NEAR(sigmoidFinite->evaluate_membership({{"x", -10}}), 1, epsilon);

  EXPECT_NEAR(sigmoidFinite->evaluate_membership({{"x", 0}}), 1, epsilon);

  EXPECT_NEAR(sigmoidFinite->evaluate_membership({{"x", 5}}), 0.5, epsilon);

  EXPECT_NEAR(sigmoidFinite->evaluate_membership({{"x", 10}}), 0, epsilon);
}

/**
 * Tests whether the negation operation on a fuzzy set calculates the correct membership values.
 * The negation of a fuzzy set is defined as 1 - membership.
 */
TEST(FuzzyTuningTest, testNegationOfFuzzySet) {
  // creates the (continuous) crisp set over which the fuzzy set is defined. It has the name "x" and the range [-5, 30]
  std::shared_ptr<CrispSet> crispSet = std::make_shared<CrispSet>("x", std::pair(-5, 30));

  auto triangle = FuzzySetFactory::makeFuzzySet("low", "Gaussian", {17, 10});
  triangle->setCrispSet(crispSet);

  auto negatedTriangle = !triangle;

  for (double x = -5; x <= 30; x += 5) {
    EXPECT_NEAR(negatedTriangle->evaluate_membership({{"x", x}}), 1 - triangle->evaluate_membership({{"x", x}}), 1e-10);
  }
}

/**
 * Tests whether the intersection operation on two fuzzy calculates the correct membership values.
 * The intersection of two fuzzy sets is defined as the minimum of the membership values of the two fuzzy sets.
 */
TEST(FuzzyTuningTest, testIntersectionOfFuzzySets) {
  // creates the (continuous) crisp set over which the fuzzy set is defined. It has the name "x" and the range [-5, 30]
  std::shared_ptr<CrispSet> crispSet = std::make_shared<CrispSet>("x", std::pair(-5, 30));

  auto triangle1 = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 5, 10});
  triangle1->setCrispSet(crispSet);

  auto triangle2 = FuzzySetFactory::makeFuzzySet("high", "Triangle", {5, 10, 15});
  triangle2->setCrispSet(crispSet);

  auto intersection = triangle1 && triangle2;

  for (double x = -5; x <= 30; x += 5) {
    EXPECT_NEAR(intersection->evaluate_membership({{"x", x}}),
                std::min(triangle1->evaluate_membership({{"x", x}}), triangle2->evaluate_membership({{"x", x}})),
                1e-10);
  }
}

/**
 * Tests whether the union operation on two fuzzy sets calculates the correct membership values.
 * The union of two fuzzy sets is defined as the maximum of the membership values of the two fuzzy sets.
 */
TEST(FuzzyTuningTest, testUnionOfFuzzySets) {
  // creates the (continuous) crisp set over which the fuzzy set is defined. It has the name "x" and the range [-5, 30]
  std::shared_ptr<CrispSet> crispSet = std::make_shared<CrispSet>("x", std::pair(-5, 30));

  auto triangle1 = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 5, 10});
  triangle1->setCrispSet(crispSet);

  auto triangle2 = FuzzySetFactory::makeFuzzySet("high", "Triangle", {5, 10, 15});
  triangle2->setCrispSet(crispSet);

  auto unionSet = triangle1 || triangle2;

  for (double x = -5; x <= 30; x += 5) {
    EXPECT_NEAR(unionSet->evaluate_membership({{"x", x}}),
                std::max(triangle1->evaluate_membership({{"x", x}}), triangle2->evaluate_membership({{"x", x}})),
                1e-10);
  }
}

/**
 * Tests whether applying a fuzzy rule correctly applies the cut operator on the consequent fuzzy set.
 * Evaluating the antecedent fuzzy set at x = 5 should yield 0.5 and therfore the consequent fuzzy set should also be
 * cut at 0.5.
 */
TEST(FuzzyTuningTest, testFuzzyRule) {
  std::shared_ptr<CrispSet> crispSetAntecedent = std::make_shared<CrispSet>("x", std::pair(-5, 30));
  std::shared_ptr<CrispSet> crispSetConsequent = std::make_shared<CrispSet>("y", std::pair(-5, 30));

  auto antecedent = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  antecedent->setCrispSet(crispSetAntecedent);

  auto consequent = FuzzySetFactory::makeFuzzySet("high", "Gaussian", {5, 10});
  consequent->setCrispSet(crispSetConsequent);

  auto fuzzyRule = FuzzyRule(antecedent, consequent);

  const auto epsilon = 1e-10;

  // Evaluating the antecedent fuzzy set at x = 5 should yield 0.5
  const std::map<std::string, double> x5 = {{"x", 5}};
  EXPECT_NEAR(antecedent->evaluate_membership(x5), 0.5, epsilon);

  // Applying the fuzzy rule should cut the consequent fuzzy set at 0.5
  auto cutConsequent = fuzzyRule.apply(x5);

  // the resulting set should behave identical to the consequent set, but with limited membership values
  for (double x = -5; x <= 30; x += 5) {
    double wouldBeMembership = consequent->evaluate_membership({{"y", x}});
    double actualMembership = cutConsequent->evaluate_membership({{"y", x}});

    EXPECT_NEAR(actualMembership, std::min(wouldBeMembership, 0.5), epsilon);
  }
}

/**
 * Tests whether the operator== of the LinguisticVariable class and the operator&& of the FuzzySet class work as
 * expected. The operator== of the LinguisticVariable should select the fuzzy set with the given name and the operator&&
 * of the FuzzySet class should return the intersection of the two fuzzy sets.
 */
TEST(FuzzyTuningTest, testLinguisticVariable) {
  auto X1 = LinguisticVariable("x1", std::pair(-5, 30));
  auto x1Low = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto x1High = FuzzySetFactory::makeFuzzySet("high", "Triangle", {15, 25, 30});

  X1.addLinguisticTerm(x1Low);
  X1.addLinguisticTerm(x1High);

  auto X2 = LinguisticVariable("x2", std::pair(-5, 30));
  auto x2Low = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto x2High = FuzzySetFactory::makeFuzzySet("high", "Triangle", {15, 25, 30});

  X2.addLinguisticTerm(x2Low);
  X2.addLinguisticTerm(x2High);

  auto combinedMembership = X1 == "low" && X2 == "high";

  const auto epsilon = 1e-10;

  for (double x1 = -5; x1 <= 30; x1 += 5) {
    for (double x2 = -5; x2 <= 30; x2 += 5) {
      double expectedMembership =
          std::min(x1Low->evaluate_membership({{"x1", x1}}), x2High->evaluate_membership({{"x2", x2}}));
      double actualMembership = combinedMembership->evaluate_membership({{"x1", x1}, {"x2", x2}});
      EXPECT_NEAR(actualMembership, expectedMembership, epsilon);
    }
  }
}

/**
 * This test is similar to the `testFuzzyRule` test but uses linguistic variables instead of fuzzy sets.
 * The test checks whether the overloaded operators behave as expected.
 * The example uses the example rule:  IF X1 == "low" && X2 == "high" THEN Y == "low"
 */
TEST(FuzzyTuningTest, testFuzzyRuleWithLinguisticVariables) {
  auto X1 = LinguisticVariable("x1", std::pair(-5, 30));
  auto x1Low = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto x1High = FuzzySetFactory::makeFuzzySet("high", "Triangle", {15, 25, 30});

  X1.addLinguisticTerm(x1Low);
  X1.addLinguisticTerm(x1High);

  auto X2 = LinguisticVariable("x2", std::pair(-5, 30));
  auto x2Low = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto x2High = FuzzySetFactory::makeFuzzySet("high", "Triangle", {15, 25, 30});

  X2.addLinguisticTerm(x2Low);
  X2.addLinguisticTerm(x2High);

  auto Y = LinguisticVariable("y", std::pair(-5, 30));
  auto yLow = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto yHigh = FuzzySetFactory::makeFuzzySet("high", "Triangle", {15, 25, 30});

  Y.addLinguisticTerm(yLow);
  Y.addLinguisticTerm(yHigh);

  auto rule1 = FuzzyRule(X1 == "low" && X2 == "high", Y == "low");

  const auto epsilon = 1e-10;

  for (double x1 = -5; x1 <= 30; x1 += 5) {
    for (double x2 = -5; x2 <= 30; x2 += 5) {
      for (double y = -5; y <= 30; y += 5) {
        double antecedentMembership =
            std::min(x1Low->evaluate_membership({{"x1", x1}}), x2High->evaluate_membership({{"x2", x2}}));

        double expectedMembership = std::min(antecedentMembership, yLow->evaluate_membership({{"y", y}}));

        auto cutConsequent = rule1.apply({{"x1", x1}, {"x2", x2}});

        double actualMembership = cutConsequent->evaluate_membership({{"y", y}});
        EXPECT_NEAR(actualMembership, expectedMembership, epsilon);
      }
    }
  }
}

/**
 * Tests whether the fuzzy controller correctly applies the fuzzy rules.
 * The test uses two fuzzy rules:
 * 1. IF X1 == "low" && X2 == "high" THEN Y == "low"
 * 2. IF X1 == "high" && X2 == "low" THEN Y == "high"
 *
 * Both rules get combined in a fuzzy controller and the controller is tested for different input values.
 * Both antecedents are evaluated and the consequent is cut at the minimum of the antecedent memberships.
 *
 * Then the consequents are combined using the Fuzzy-OR operator and the resulting fuzzy set should therefore correspond
 * to the maximum of the two cut consequents when evualated.
 */
TEST(FuzzyTuningTest, testFuzzyController) {
  auto X1 = LinguisticVariable("x1", std::pair(-5, 30));
  auto x1Low = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto x1High = FuzzySetFactory::makeFuzzySet("high", "Triangle", {15, 25, 30});

  X1.addLinguisticTerm(x1Low);
  X1.addLinguisticTerm(x1High);

  auto X2 = LinguisticVariable("x2", std::pair(-5, 30));
  auto x2Low = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto x2High = FuzzySetFactory::makeFuzzySet("high", "Triangle", {15, 25, 30});

  X2.addLinguisticTerm(x2Low);
  X2.addLinguisticTerm(x2High);

  auto Y = LinguisticVariable("y", std::pair(-5, 30));
  auto yLow = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto yHigh = FuzzySetFactory::makeFuzzySet("high", "Triangle", {15, 25, 30});

  Y.addLinguisticTerm(yLow);
  Y.addLinguisticTerm(yHigh);

  auto rule1 = FuzzyRule(X1 == "low" && X2 == "high", Y == "low");
  auto rule2 = FuzzyRule(X1 == "high" && X2 == "low", Y == "high");

  auto fuzzyController = FuzzyControlSystem({});

  fuzzyController.addRule(rule1);
  fuzzyController.addRule(rule2);

  const auto epsilon = 1e-10;

  for (double x1 = -5; x1 <= 30; x1 += 5) {
    for (double x2 = -5; x2 <= 30; x2 += 5) {
      for (double y = -5; y <= 30; y += 5) {
        double antecedent1Membership =
            std::min(x1Low->evaluate_membership({{"x1", x1}}), x2High->evaluate_membership({{"x2", x2}}));
        double antecedent2Membership =
            std::min(x1High->evaluate_membership({{"x1", x1}}), x2Low->evaluate_membership({{"x2", x2}}));

        double expectedMembership1 = std::min(antecedent1Membership, yLow->evaluate_membership({{"y", y}}));
        double expectedMembership2 = std::min(antecedent2Membership, yHigh->evaluate_membership({{"y", y}}));

        auto cutConsequent = fuzzyController.applyRules({{"x1", x1}, {"x2", x2}});
        double actualMembership = cutConsequent->evaluate_membership({{"y", y}});

        EXPECT_NEAR(actualMembership, std::max(expectedMembership1, expectedMembership2), epsilon);
      }
    }
  }
}

/**
 * Tests if the center of gravity defuzzification method works as expected.
 * The test uses different triangular and gaussian fuzzy sets and checks if the center of gravity of combinations of
 * these fuzzy sets is calculated correctly.
 */
TEST(FuzzyTuningTest, testCOGDefuzzification) {
  auto crispSet = std::make_shared<CrispSet>("x", std::pair(0, 40));

  auto t1 = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto t2 = FuzzySetFactory::makeFuzzySet("mid", "Triangle", {10, 20, 30});
  auto t3 = FuzzySetFactory::makeFuzzySet("high", "Triangle", {20, 30, 40});

  auto g1 = FuzzySetFactory::makeFuzzySet("high", "Gaussian", {15, 5});

  // add crisp set to fuzzy sets
  t1->setCrispSet(crispSet);
  t2->setCrispSet(crispSet);
  t3->setCrispSet(crispSet);
  g1->setCrispSet(crispSet);

  size_t numSamples = 1000;
  // the center position of a triangle should be below the peak if the triangle is symmetric
  EXPECT_NEAR(t1->defuzzify(DefuzzificationMethodOption::CoG, numSamples), 10, 1e-1);

  // the center position of a gaussian should be at the peak
  EXPECT_NEAR(g1->defuzzify(DefuzzificationMethodOption::CoG, numSamples), 15, 1e-1);

  // center of gravity of the union of the low and the high fuzzy set should be at 20 since its the perfect center
  auto unionSet = t1 || t3;
  EXPECT_NEAR(unionSet->defuzzify(DefuzzificationMethodOption::CoG, numSamples), 20, 1e-1);
}

/**
 * Tests if the mean of maximum defuzzification method works as expected.
 * The test uses different triangular and gaussian fuzzy sets and checks if the mean of maximum of combinations of
 * these fuzzy sets is calculated correctly.
 */
TEST(FuzzyTuningTest, testMaxDefuzzification) {
  auto crispSet = std::make_shared<CrispSet>("x", std::pair(0, 40));

  auto t1 = FuzzySetFactory::makeFuzzySet("low", "Triangle", {0, 10, 20});
  auto t2 = FuzzySetFactory::makeFuzzySet("mid", "Triangle", {10, 20, 30});
  auto t3 = FuzzySetFactory::makeFuzzySet("high", "Triangle", {20, 30, 40});

  auto g1 = FuzzySetFactory::makeFuzzySet("high", "Gaussian", {15, 5});

  // add crisp set to fuzzy sets
  t1->setCrispSet(crispSet);
  t2->setCrispSet(crispSet);
  t3->setCrispSet(crispSet);
  g1->setCrispSet(crispSet);

  size_t numSamples = 1000;
  // the max position of a triangle should be at the peak
  EXPECT_NEAR(t1->defuzzify(DefuzzificationMethodOption::MoM, numSamples), 10, 1e-1);

  // the max position of a gaussian should be at the peak
  EXPECT_NEAR(g1->defuzzify(DefuzzificationMethodOption::MoM, numSamples), 15, 1e-1);

  // cut t2
  auto rule1 = FuzzyRule(t1, t2);
  auto cutT2 = rule1.apply({{"x", 5}});

  // cut t3
  auto rule2 = FuzzyRule(t2, t3);
  auto cutT3 = rule2.apply({{"x", 15}});

  // the max position of the union of the low and the cut of the other fuzzy sets should be close to the peak of t1
  auto unionSet = t1 || cutT2 || cutT3;
  EXPECT_NEAR(unionSet->defuzzify(DefuzzificationMethodOption::MoM, numSamples), 10, 1e-1);
}

/**
 * Tests if the fuzzy tuning class can parse a rule file and create the correct fuzzy control system.
 * The test uses a simple rule file with one rule and checks if the rule file can be parsed and interpreted correctly.
 */
TEST(FuzzyTuningTest, testParseRuleFile) {
  const std::string fileContent = R"(
# Define the settings of the fuzzy control system
FuzzySystemSettings:
     defuzzificationMethod: "meanOfMaximum"
     interpretOutputAs: "IndividualSystems"


# Define all of the linguistic variables together with their linguistic terms
FuzzyVariable: domain: "homogeneity" range: (-0.009, 0.1486)
     "lower than 0.049":     SigmoidFinite(0.0914, 0.049, 0.0065)
     "lower than 0.041":     SigmoidFinite(0.0834, 0.041, -0.001)
     "higher than 0.049":    SigmoidFinite(0.0065, 0.049, 0.0914)
     "higher than 0.041":    SigmoidFinite(-0.001, 0.041, 0.0834)

FuzzyVariable: domain: "particlesPerCellStdDev" range: (-0.017, 0.072)
     "lower than 0.013":     SigmoidFinite(0.0639, 0.038,  0.012)
     "lower than 0.014":     SigmoidFinite(0.0399, 0.014, -0.011)
     "higher than 0.013":    SigmoidFinite(0.012,  0.013,  0.0639)
     "higher than 0.014":    SigmoidFinite(-0.011, 0.014,  0.0399)

FuzzyVariable: domain: "Newton 3" range: (0, 1)
      "disabled, enabled":   Gaussian(0.3333, 0.1667)
      "enabled":             Gaussian(0.6667, 0.1667)

# Define how the result of the output variables should be interpreted in the context of autopas
OutputMapping:
 "Newton 3":
     0.333 => [newton3="disabled"], [newton3="enabled"]
     0.666 => [newton3="enabled"]

# Define a bunch of rules connecting the input variables to the output variables
if ("homogeneity" == "lower than 0.041") && ("particlesPerCellStdDev" == "lower than 0.014")
   then ("Newton 3" == "enabled")
)";

  // make temporary file in /tmp
  const std::string fileName = std::filesystem::temp_directory_path().string() + "/fuzzyRules.frule";
  std::ofstream file(fileName);
  file << fileContent;
  // make sure the file is written before we use it.
  file.close();

  // parse the file
  autopas::FuzzyTuning tuner{fileName};

  auto settings = tuner.getFuzzyControlSettings();
  auto outputMappings = tuner.getOutputMappings();
  auto fuzzySystems = tuner.getFuzzyControlSystems();

  // check settings
  EXPECT_EQ(settings->size(), 2);
  EXPECT_EQ(settings->at("defuzzificationMethod"), "meanOfMaximum");
  EXPECT_EQ(settings->at("interpretOutputAs"), "IndividualSystems");

  auto newtonMapper = outputMappings.at("Newton 3");

  // check output mappings
  EXPECT_EQ(outputMappings.size(), 1);
  EXPECT_EQ(newtonMapper->getOutputDomain(), "Newton 3");

  auto entries = newtonMapper->getMappings();
  EXPECT_EQ(entries.size(), 2);
  EXPECT_NEAR(entries[0].first, 0.333, 1e-10);
  EXPECT_EQ(entries[0].second.size(), 2);

  auto configsAt0333 = newtonMapper->getClosestConfigurationPatterns(0.333);
  EXPECT_EQ(configsAt0333.size(), 2);

  EXPECT_EQ(configsAt0333.size(), 2);
  auto firstEntry = configsAt0333[0];

  EXPECT_EQ(firstEntry._newton3Options.count(autopas::Newton3Option(autopas::Newton3Option::Value::disabled)), 1);

  auto secondEntry = configsAt0333[1];

  EXPECT_EQ(secondEntry._newton3Options.count(autopas::Newton3Option(autopas::Newton3Option::Value::enabled)), 1);

  auto configsAt0666 = newtonMapper->getClosestConfigurationPatterns(0.666);
  EXPECT_EQ(configsAt0666.size(), 1);

  auto firstEntryAt0666 = configsAt0666[0];

  EXPECT_EQ(firstEntryAt0666._newton3Options.count(autopas::Newton3Option(autopas::Newton3Option::Value::enabled)), 1);

  // check fuzzy systems
  EXPECT_EQ(fuzzySystems.size(), 1);

  auto newtonSystem = fuzzySystems.at("Newton 3");

  auto expectedSystemString =
      "FuzzyControlSystem: \"Newton 3\"\n"
      "\tFuzzyRule: (\"homogeneity\" == \"lower than 0.041\" && \"particlesPerCellStdDev\" == "
      "\"lower than 0.014\")  ==>  \"Newton 3\" == \"enabled\"\n";

  EXPECT_EQ(std::string(*newtonSystem), expectedSystemString);

  // Predict some data
  std::map<std::string, double> data = {{"homogeneity", 0.04}, {"particlesPerCellStdDev", 0.013}};

  auto result = newtonSystem->predict(data);
  EXPECT_NEAR(result, 0.6667, 1e-2);
}
#endif  // AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING

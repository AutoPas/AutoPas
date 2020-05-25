/**
 * @file GaussianProcessTest.cpp
 * @author Jan Nguyen
 * @date 12.06.19
 */

#include "GaussianProcessTest.h"

using namespace autopas;

TEST_F(GaussianProcessTest, wrongDimension) {
  Random rng(32);
  GaussianProcess gp(2, 0.001, rng);

  Eigen::VectorXd f1 = Eigen::VectorXd::Ones(1);
  Eigen::VectorXd f2 = Eigen::VectorXd::Zero(3);

  EXPECT_THROW(gp.addEvidence(f1, 0, true), utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(gp.addEvidence(f2, 1, true), utils::ExceptionHandler::AutoPasException);
}

TEST_F(GaussianProcessTest, noEvidence) {
  Random rng(32);
  double epsilon = 0.05;  // allowed error for tests

  double sigma = 0.001;  // noise
  GaussianProcess gp(1, sigma, rng);

  auto f1 = autopas::utils::Math::makeVectorXd({0.});

  // predict without information -> should return default values
  EXPECT_NEAR(gp.predictMean(f1), 0., epsilon);
  EXPECT_NEAR(gp.predictVar(f1), 1., epsilon);
}

TEST_F(GaussianProcessTest, oneEvidence) {
  Random rng(32);

  double epsilon = 0.05;
  double sigma = 0.001;
  GaussianProcess gp(1, sigma, rng);

  auto f1 = autopas::utils::Math::makeVectorXd({0.});
  double out1 = 42.;

  auto f2 = autopas::utils::Math::makeVectorXd({1000.});

  gp.addEvidence(f1, out1, true);

  // predicting point same as evidence -> should expect same output as evidence
  EXPECT_NEAR(gp.predictMean(f1), out1, epsilon);
  EXPECT_NEAR(gp.predictVar(f1), 0., epsilon);

  // not enough information, so other inputs lead to same mean
  EXPECT_NEAR(gp.predictMean(f2), out1, epsilon);
}

TEST_F(GaussianProcessTest, twoEvidence) {
  Random rng;

  double epsilon = 0.05;  // allowed error for tests
  double sigma = 0.001;   // noise
  GaussianProcess gp(1, sigma, rng);

  auto f1 = autopas::utils::Math::makeVectorXd({-100.});
  double out1 = 42.;

  auto f2 = autopas::utils::Math::makeVectorXd({100.});
  double out2 = -3.;

  gp.addEvidence(f1, out1, false);
  gp.addEvidence(f2, out2, true);

  // predicting point same as evidence
  // should expect same output as evidence because great distance between inputs
  EXPECT_NEAR(gp.predictMean(f1), out1, epsilon);
  EXPECT_NEAR(gp.predictVar(f1), 0., epsilon);

  EXPECT_NEAR(gp.predictMean(f2), out2, epsilon);
  EXPECT_NEAR(gp.predictVar(f2), 0., epsilon);
}

TEST_F(GaussianProcessTest, clear) {
  Random rng;

  double epsilon = 0.05;  // allowed error for tests
  double sigma = 0.001;   // noise
  GaussianProcess gp(1, sigma, rng);

  auto f1 = autopas::utils::Math::makeVectorXd({-100.});
  double out1 = 42.;

  auto f2 = autopas::utils::Math::makeVectorXd({100.});
  double out2 = -3.;

  gp.addEvidence(f1, out1, false);
  gp.addEvidence(f2, out2, true);
  gp.clear();

  // predicting points as deleted evidence
  // they should not have effect on the prediction anymore
  EXPECT_NEAR(gp.predictMean(f1), 0., epsilon);
  EXPECT_NEAR(gp.predictVar(f1), 1., epsilon);

  EXPECT_NEAR(gp.predictMean(f2), 0., epsilon);
  EXPECT_NEAR(gp.predictVar(f2), 1., epsilon);

  // test new evidence at deleted point
  out2 = 10.;
  gp.addEvidence(f2, out2, true);

  EXPECT_NEAR(gp.predictMean(f1), out2, epsilon);

  EXPECT_NEAR(gp.predictMean(f2), out2, epsilon);
  EXPECT_NEAR(gp.predictVar(f2), 0., epsilon);
}

TEST_F(GaussianProcessTest, sine) {
  Random rng(42);

  // gp should try to approximate the sine as blackbox function
  auto functor = [](double input) { return std::sin(input); };
  constexpr double epsilon = 0.05;        // allowed error
  constexpr size_t numEvidence = 6;       // number of evidence to provide
  constexpr size_t numPredictions = 500;  // number of predictions to make
  constexpr double domainStart = 0.;      // start of tested domain
  constexpr double domainEnd = 2 * M_PI;  // end of tested domain

  GaussianProcess gp(1, 0.001, rng);

  // create equidistant evidence over the domain
  double evidenceStep = (domainEnd - domainStart) / (numEvidence - 1);
  for (unsigned indexFirst = 0; indexFirst < numEvidence; ++indexFirst) {
    double input = domainStart + evidenceStep * indexFirst;
    auto f = autopas::utils::Math::makeVectorXd({input});
    double output = functor(input);

    gp.addEvidence(f, output, true);
  }

  // make equidistant prediction over the domain
  double predictStep = (domainEnd - domainStart) / (numPredictions - 1);
  for (unsigned indexFirst = 0; indexFirst < numPredictions; ++indexFirst) {
    double input = domainStart + predictStep * indexFirst;
    auto f = autopas::utils::Math::makeVectorXd({input});
    double output = functor(input);

    EXPECT_NEAR(gp.predictMean(f), output, epsilon);
  }
}

TEST_F(GaussianProcessTest, 2dMax) {
  // try to find the max of -(i1 + 1)^2 - (i2 - 1)^2
  auto functor = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 1, 2); };
  std::pair domain{NumberInterval<double>(-2, 2), NumberInterval<double>(-2, 2)};
  constexpr double maxError = 0.2;

  // max of function
  auto max = autopas::utils::Math::makeVectorXd({-1, 1});

  test2DFunction(functor, max, maxError, domain, AcquisitionFunctionOption::upperConfidenceBound, false);
}

TEST_F(GaussianProcessTest, 2dMaxGrid) {
  // functor to find max of
  auto functor = [](double i1, double i2) { return -2 * std::pow(i1 - 1, 2) - 3 * std::pow(i2 - 1, 2); };
  // discrete domain of function
  int domHalf = 10;
  std::set<double> domSet;
  for (int i = -domHalf; i <= domHalf; ++i) {
    domSet.insert(i * 2. / domHalf);
  }
  std::pair domain{NumberSetFinite<double>(domSet), NumberSetFinite<double>(domSet)};
  constexpr double maxError = 0.5;

  // max of function
  auto max = autopas::utils::Math::makeVectorXd({1, 1});

  test2DFunction(functor, max, maxError, domain, AcquisitionFunctionOption::expectedImprovement, false);
}

TEST_F(GaussianProcessTest, 2dMaxGridBig) {
  // functor to find max of
  auto functor = [](double i1, double i2) {
    return std::pow(-(std::abs(i1 - 1) + 1), 5) - std::pow((std::abs(i2 - 1) + 1), 5);
  };
  // discrete domain of function
  int domHalf = 10;
  std::set<double> domSet;
  for (int i = -domHalf; i <= domHalf; ++i) {
    domSet.insert(i * 2. / domHalf);
  }
  std::pair domain{NumberSetFinite<double>(domSet), NumberSetFinite<double>(domSet)};
  constexpr double maxError = 10;

  // max of function
  auto max = autopas::utils::Math::makeVectorXd({1, 1});

  test2DFunction(functor, max, maxError, domain, AcquisitionFunctionOption::upperConfidenceBound, false);
}

void GaussianProcessTest::printMap(int xChunks, int yChunks, const autopas::NumberSet<double> &domainX,
                                   const autopas::NumberSet<double> &domainY, const autopas::GaussianProcess &gp,
                                   autopas::AcquisitionFunctionOption af, double colorFactor) {
  // get distance between chunks
  double xSpace = (domainX.getMax() - domainX.getMin()) / (xChunks - 1);
  double ySpace = (domainY.getMax() - domainY.getMin()) / (yChunks - 1);

  // precalculate acqMap
  std::vector<std::vector<double>> acqMap;
  double acqMin = std::numeric_limits<double>::max();
  double acqMax = std::numeric_limits<double>::lowest();
  for (int y = 0; y < xChunks; ++y) {
    // calculate a row
    std::vector<double> row;
    for (int x = 0; x < xChunks; ++x) {
      // calculate value of chunk
      Eigen::Vector2d sample(x * xSpace + domainX.getMin(), (y * ySpace + domainY.getMin()));
      double val = gp.calcAcquisition(af, sample);

      row.push_back(val);

      // keep track min and max value
      acqMin = std::min(val, acqMin);
      acqMax = std::max(val, acqMax);
    }

    acqMap.push_back(std::move(row));
  }

  // get scaling such that acqMax=1 and acqMin=0
  double acqScale = 1 / (acqMax - acqMin);

  // print map
  for (int y = yChunks - 1; y >= 0; --y) {
    // acquisiton row
    for (int x = 0; x < xChunks; ++x) {
      // map value between 0 to 1 and square for clearer differences of high values
      double val = std::pow((acqMap[y][x] - acqMin) * acqScale, 2);

      // map value to color
      int color = static_cast<int>(232 + 23 * val);
      color = std::clamp(color, 232, 255);

      // print two spaces of that color
      std::cout << "\033[48;5;" << color << "m  ";
    }
    // reset color, print a space
    std::cout << "\033[0m ";

    // mean row
    for (int x = 0; x < xChunks; ++x) {
      Eigen::Vector2d sample(x * xSpace + domainX.getMin(), (y * ySpace + domainY.getMin()));
      double val = gp.predictMean(sample) * colorFactor;

      // map value to color
      int color = static_cast<int>(255 + val);
      color = std::clamp(color, 232, 255);

      // print two spaces of that color
      std::cout << "\033[48;5;" << color << "m  ";
    }
    // reset color, print a space
    std::cout << "\033[0m" << std::endl;
  }
}

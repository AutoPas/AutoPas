/**
 * @file GaussianProcessTest.cpp
 * @author Jan Nguyen
 * @date 12.06.19
 */

#include "GaussianProcessTest.h"

using namespace autopas;

TEST(GaussianProcessTest, noEvidence) {
  double epsilon = 0.05;  // allowed error for tests

  double theta = 2.;     // default variance
  double sigma = 0.001;  // noise
  GaussianProcess<Eigen::VectorXd> gp(theta, {1.}, sigma);

  Eigen::VectorXd f1(1);
  f1 << 0.;

  // predict without information -> should return default values
  ASSERT_NEAR(gp.predictMean(f1), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f1), theta, epsilon);
}

TEST(GaussianProcessTest, oneEvidence) {
  double epsilon = 0.05;

  double theta = 4.;
  double sigma = 0.001;
  GaussianProcess<Eigen::VectorXd> gp(theta, {1.}, sigma);

  Eigen::VectorXd f1(1);
  f1 << 0.;
  double out1 = 42.;

  Eigen::VectorXd f2(1);
  f2 << 1000.;

  gp.addEvidence(f1, out1);

  // predicting point same as evidence -> should expect same output as evidence
  ASSERT_NEAR(gp.predictMean(f1), out1, epsilon);
  ASSERT_NEAR(gp.predictVar(f1), 0., epsilon);

  // prediction far away from evidence should expect default.
  ASSERT_NEAR(gp.predictMean(f2), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f2), theta, epsilon);
}

TEST(GaussianProcessTest, twoEvidence) {
  double epsilon = 0.05;  // allowed error for tests

  double theta = 5.;     // default variance
  double sigma = 0.001;  // noise
  GaussianProcess<Eigen::VectorXd> gp(theta, {1.}, sigma);

  Eigen::VectorXd f1(1);
  f1 << -100.;
  double out1 = 42.;

  Eigen::VectorXd f2(1);
  f2 << 100.;
  double out2 = -3.;

  Eigen::VectorXd f3(1);
  f3 << 0.;

  gp.addEvidence(f1, out1);
  gp.addEvidence(f2, out2);

  // predicting point same as evidence
  // should expect same output as evidence because greate distance between inputs
  ASSERT_NEAR(gp.predictMean(f1), out1, epsilon);
  ASSERT_NEAR(gp.predictVar(f1), 0., epsilon);

  ASSERT_NEAR(gp.predictMean(f2), out2, epsilon);
  ASSERT_NEAR(gp.predictVar(f2), 0., epsilon);

  // prediction far away from evidence should expect default.
  ASSERT_NEAR(gp.predictMean(f3), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f3), theta, epsilon);
}

TEST(GaussianProcessTest, clear) {
  double epsilon = 0.05;  // allowed error for tests

  double theta = 5.;     // default variance
  double sigma = 0.001;  // noise
  GaussianProcess<Eigen::VectorXd> gp(theta, {1.}, sigma);

  Eigen::VectorXd f1(1);
  f1 << -100.;
  double out1 = 42.;

  Eigen::VectorXd f2(1);
  f2 << 100.;
  double out2 = -3.;

  gp.addEvidence(f1, out1);
  gp.addEvidence(f2, out2);
  gp.clear();

  // predicting points as deleted evidence
  // they should not have effect on the prediction anymore
  ASSERT_NEAR(gp.predictMean(f1), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f1), theta, epsilon);

  ASSERT_NEAR(gp.predictMean(f2), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f2), theta, epsilon);

  // test new evidence at deleted point
  out2 = 10.;
  gp.addEvidence(f2, out2);

  ASSERT_NEAR(gp.predictMean(f1), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f1), theta, epsilon);

  ASSERT_NEAR(gp.predictMean(f2), out2, epsilon);
  ASSERT_NEAR(gp.predictVar(f2), 0., epsilon);
}

TEST(GaussianProcessTest, sine) {
  // gp should try to approximate the sine as blackbox function
  auto functor = [](double input) { return std::sin(input); };
  double epsilon = 0.02;           // allowed error
  unsigned numEvidence = 9;        // number of evidence to provide
  unsigned numPredictions = 1000;  // number of predictions to make
  double domainStart = 0.;         // start of tested domain
  double domainEnd = 2 * M_PI;     // end of tested domain

  GaussianProcess<Eigen::VectorXd> gp(0.5, {0.2}, 0.001);

  // create equidistant evidence over the domain
  double evidenceStep = (domainEnd - domainStart) / (numEvidence - 1);
  for (unsigned i = 0; i < numEvidence; ++i) {
    double input = domainStart + evidenceStep * i;
    Eigen::VectorXd f(1);
    f << input;
    double output = functor(input);

    gp.addEvidence(f, output);
  }

  // make equidistant prediction over the domain
  double predictStep = (domainEnd - domainStart) / (numPredictions - 1);
  for (unsigned i = 0; i < numPredictions; ++i) {
    double input = domainStart + predictStep * i;
    Eigen::VectorXd f(1);
    f << input;
    double output = functor(input);

    ASSERT_NEAR(gp.predictMean(f), output, epsilon);
  }
}

TEST(GaussianProcessTest, 2dMax) {
  Random rng(72);  // random generator

  // try to find the max of -(i1 + 1)^2 - (i2 - 1)^2
  auto functor = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 1, 2); };
  double epsilon = 0.05;  // allowed error
  std::vector<NumberInterval<double>> domain{NumberInterval<double>(-2, 2),
                                             NumberInterval<double>(-2, 2)};  // domain of function

  // max of function
  Eigen::VectorXd max(2);
  max << -1, 1;
  unsigned numEvidence = 40;      // number of samples allowed to make
  unsigned lhsNumSamples = 1000;  // number of sample to find max of acquisition function
  AcquisitionFunctionOption af = AcquisitionFunctionOption::ucb;      // use upper confidence bound as af
  AcquisitionFunctionOption lastAf = AcquisitionFunctionOption::lcb;  // use lower confidence bound for final prediction

  GaussianProcess<Eigen::VectorXd> gp(6, {0.2, 0.2}, 0.001);

  // add first evidence
  Eigen::VectorXd first(2);
  first << 0, 0;
  gp.addEvidence(first, functor(0, 0));

  for (unsigned i = 1; i < numEvidence; ++i) {
    // create lhs samples
    std::vector<Eigen::VectorXd> lhsSamples;
    lhsSamples.reserve(lhsNumSamples);

    auto xSamples = domain[0].uniformSample(lhsNumSamples, rng);
    auto ySamples = domain[1].uniformSample(lhsNumSamples, rng);
    for (size_t i = 0; i < lhsNumSamples; ++i) {
      Eigen::VectorXd sample(2);
      sample << xSamples[i], ySamples[i];
      lhsSamples.push_back(sample);
    }

    // sample max of acquisition function
    Eigen::VectorXd am = gp.sampleAquisitionMax(af, lhsSamples);
    double amOut = functor(am[0], am[1]);

    gp.addEvidence(am, amOut);
  }

  // last lhs sample for predicted max
  auto xSamples = domain[0].uniformSample(lhsNumSamples, rng);
  auto ySamples = domain[1].uniformSample(lhsNumSamples, rng);
  std::vector<Eigen::VectorXd> lhsSamples;
  lhsSamples.reserve(lhsNumSamples);
  for (size_t i = 0; i < lhsNumSamples; ++i) {
    Eigen::VectorXd sample(2);
    sample << xSamples[i], ySamples[i];
    lhsSamples.push_back(sample);
  }

  // final sample
  Eigen::VectorXd am = gp.sampleAquisitionMax(lastAf, lhsSamples);

  // check if predicted max is near real max
  double predMax = functor(am[0], am[1]);
  double realMax = functor(max[0], max[1]);
  ASSERT_NEAR(predMax, realMax, epsilon);
}

TEST(GaussianProcessTest, 2dMin) {
  Random rng(73);  // random generator

  // try to find the min of (i1 - 1)^2 + (i2 - 1)^2
  auto functor = [](double i1, double i2) { return std::pow(i1 - 1, 2) + std::pow(i2 - 1, 2); };
  double epsilon = 0.05;  // allowed error
  std::vector<NumberInterval<double>> domain{NumberInterval<double>(-2, 2),
                                             NumberInterval<double>(-2, 2)};  // domain of function

  // min of function
  Eigen::VectorXd min(2);
  min << 1, 1;
  unsigned numEvidence = 40;      // number of samples allowed to make
  unsigned lhsNumSamples = 1000;  // number of sample to find min of acquisition function
  AcquisitionFunctionOption af = AcquisitionFunctionOption::lcb;      // use lower confidence bound as af
  AcquisitionFunctionOption lastAf = AcquisitionFunctionOption::ucb;  // use upper confidence bound for final prediction

  GaussianProcess<Eigen::VectorXd> gp(6, {0.2, 0.2}, 0.001);

  // add first evidence
  Eigen::VectorXd first(2);
  first << 0, 0;
  gp.addEvidence(first, functor(0, 0));

  for (unsigned i = 1; i < numEvidence; ++i) {
    // create lhs samples
    std::vector<Eigen::VectorXd> lhsSamples;
    lhsSamples.reserve(lhsNumSamples);

    auto xSamples = domain[0].uniformSample(lhsNumSamples, rng);
    auto ySamples = domain[1].uniformSample(lhsNumSamples, rng);
    for (size_t i = 0; i < lhsNumSamples; ++i) {
      Eigen::VectorXd sample(2);
      sample << xSamples[i], ySamples[i];
      lhsSamples.push_back(sample);
    }

    // sample min of acquisition function
    Eigen::VectorXd am = gp.sampleAquisitionMin(af, lhsSamples);
    double amOut = functor(am[0], am[1]);

    gp.addEvidence(am, amOut);
  }

  // last lhs sample for predicted min
  auto xSamples = domain[0].uniformSample(lhsNumSamples, rng);
  auto ySamples = domain[1].uniformSample(lhsNumSamples, rng);
  std::vector<Eigen::VectorXd> lhsSamples;
  lhsSamples.reserve(lhsNumSamples);
  for (size_t i = 0; i < lhsNumSamples; ++i) {
    Eigen::VectorXd sample(2);
    sample << xSamples[i], ySamples[i];
    lhsSamples.push_back(sample);
  }

  // final sample
  Eigen::VectorXd am = gp.sampleAquisitionMin(lastAf, lhsSamples);

  // check if predicted min is near real min
  double predMin = functor(am[0], am[1]);
  double realMin = functor(min[0], min[1]);
  ASSERT_NEAR(predMin, realMin, epsilon);
}

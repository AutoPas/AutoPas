/**
 * @file ArrayMathTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include "autopas/selectors/GaussianProcess.h"
#include "autopas/utils/DoubleSet.h"

using namespace autopas;

TEST(GaussianProcess, noEvidence) {
  double epsilon = 0.05;  // allowed error for tests

  double theta = 2.;     // default variance
  double sigma = 0.001;  // noise
  GaussianProcess gp(theta, {1., 1.}, sigma);

  FeatureVector f1;

  // predict without information -> should return default values
  ASSERT_NEAR(gp.predictMean(f1), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f1), theta, epsilon);
}

TEST(GaussianProcess, oneEvidence) {
  double epsilon = 0.05;

  double theta = 4.;
  double sigma = 0.001;
  GaussianProcess gp(theta, {1., 1.}, sigma);

  FeatureVector f1;
  f1.addFeature(1.);
  f1.addFeature(1.);
  double out1 = 42.;

  FeatureVector f2;
  f2.addFeature(100.);
  f2.addFeature(100.);

  gp.addEvidence(f1, out1);

  // predicting point same as evidence -> should expect same output as evidence
  ASSERT_NEAR(gp.predictMean(f1), out1, epsilon);
  ASSERT_NEAR(gp.predictVar(f1), 0., epsilon);

  // prediction far away from evidence should expect default.
  ASSERT_NEAR(gp.predictMean(f2), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f2), theta, epsilon);
}

TEST(GaussianProcess, twoEvidences) {
  double epsilon = 0.05;  // allowed error for tests

  double theta = 5.;     // default variance
  double sigma = 0.001;  // noise
  GaussianProcess gp(theta, {1., 1.}, sigma);

  FeatureVector f1;
  f1.addFeature(-100.);
  f1.addFeature(-100.);
  double out1 = 42.;

  FeatureVector f2;
  f2.addFeature(100.);
  f2.addFeature(100.);
  double out2 = -3.;

  FeatureVector f3;
  f3.addFeature(0.);
  f3.addFeature(0.);

  gp.addEvidence(f1, out1);
  gp.addEvidence(f2, out2);

  // predicting point same as evidences
  // should expect same output as evidence because greate distance between inputs
  ASSERT_NEAR(gp.predictMean(f1), out1, epsilon);
  ASSERT_NEAR(gp.predictVar(f1), 0., epsilon);

  ASSERT_NEAR(gp.predictMean(f2), out2, epsilon);
  ASSERT_NEAR(gp.predictVar(f2), 0., epsilon);

  // prediction far away from evidence should expect default.
  ASSERT_NEAR(gp.predictMean(f3), 0., epsilon);
  ASSERT_NEAR(gp.predictVar(f3), theta, epsilon);
}

TEST(GaussianProcess, clear) {
  double epsilon = 0.05;  // allowed error for tests

  double theta = 5.;     // default variance
  double sigma = 0.001;  // noise
  GaussianProcess gp(theta, {1., 1.}, sigma);

  FeatureVector f1;
  f1.addFeature(-100.);
  f1.addFeature(-100.);
  double out1 = 42.;

  FeatureVector f2;
  f2.addFeature(100.);
  f2.addFeature(100.);
  double out2 = -3.;

  gp.addEvidence(f1, out1);
  gp.addEvidence(f2, out2);
  gp.clear();

  // predicting points as deleted evidences
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

TEST(GaussianProcess, sine) {
  // gp should try to approximate the sine as blackbox function
  auto functor = [](double input) { return std::sin(input); };
  double epsilon = 0.02;           // allowed error
  unsigned numEvidences = 9;       // number of evidences to provide
  unsigned numPredictions = 1000;  // number of predictions to make
  double domainStart = 0.;         // start of tested domain
  double domainEnd = 2 * M_PI;     // end of tested domain

  GaussianProcess gp(0.5, {0.2}, 0.001);

  // create equidistant evidence over the domain
  double evidenceStep = (domainEnd - domainStart) / (numEvidences - 1);
  for (unsigned i = 0; i < numEvidences; ++i) {
    double input = domainStart + evidenceStep * i;
    FeatureVector f({input});
    double output = functor(input);

    gp.addEvidence(f, output);
  }

  // make equidistant prediction over the domain
  double predictStep = (domainEnd - domainStart) / (numPredictions - 1);
  for (unsigned i = 0; i < numPredictions; ++i) {
    double input = domainStart + predictStep * i;
    FeatureVector f({input});
    double output = functor(input);

    ASSERT_NEAR(gp.predictMean(f), output, epsilon);
  }
}

TEST(GaussianProcess, 2dMax) {
  std::default_random_engine rng(72);         // random generator

  // try to find the max of -(i1 + 1)^2 - (i2 - 1)^2
  auto functor = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 1, 2); };
  double epsilon = 0.1;                                                     // allowed error
  std::vector<DoubleInterval> domain{DoubleInterval(-2, 2), DoubleInterval(-2, 2)};  // domain of function
  FeatureVector max({-1, 1});                                               // max of function
  unsigned numEvidences = 20;                                               // number of samples allowed to make
  unsigned lhsNumSamples = 1000;                      // number of sample to find max of acquisition function
  AcquisitionFunction af = AcquisitionFunction::ucb;  // use upper confidence bound as af

  GaussianProcess gp(6, {0.2, 0.2}, 0.001);

  // add first evidence
  gp.addEvidence({0, 0}, functor(0, 0));

  for (unsigned i = 1; i < numEvidences; ++i) {
    // create lhs samples
    std::vector<FeatureVector> lhsSamples(lhsNumSamples);
    for (auto& d : domain) FeatureVector::lhsAddFeature(lhsSamples, d, rng);

    // sample max of acquisition function
    FeatureVector am = gp.sampleAquisitionMax(af, lhsSamples);
    double i1 = static_cast<DoubleFeature&>(am[0]).getValue();
    double i2 = static_cast<DoubleFeature&>(am[1]).getValue();
    double amOut = functor(i1, i2);

    gp.addEvidence(am, amOut);
  }

  // last lhs sample for predicted max
  std::vector<FeatureVector> lhsSamples(lhsNumSamples);
  for (auto& d : domain) FeatureVector::lhsAddFeature(lhsSamples, d, rng);

  // sample max of mean function
  FeatureVector am = gp.sampleAquisitionMax(AcquisitionFunction::mean, lhsSamples);

  // check if predicted max is near real max
  auto diff = am - max;
  ASSERT_NEAR(diff[0], 0, epsilon);
  ASSERT_NEAR(diff[1], 0, epsilon);
}

/**
 * @file GaussianProcessTest.cpp
 * @author Jan Nguyen
 * @date 12.06.19
 */

#include "GaussianProcessTest.h"
#include "autopas/utils/Random.h"

using namespace autopas;

TEST(GaussianProcessTest, wrongDimension) {
  Random rng(32);
  GaussianProcess<Eigen::VectorXd> gp(2, 0.001, rng);

  Eigen::VectorXd f1 = Eigen::VectorXd::Ones(1);
  Eigen::VectorXd f2 = Eigen::VectorXd::Zero(3);

  EXPECT_THROW(gp.addEvidence(f1, 0), utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(gp.addEvidence(f2, 1), utils::ExceptionHandler::AutoPasException);
}

TEST(GaussianProcessTest, distanceTest) {
  autopas::FeatureVector f1(ContainerOption::linkedCells, 1., TraversalOption::c01, DataLayoutOption::aos,
                            Newton3Option::enabled);
  autopas::FeatureVector f2(ContainerOption::linkedCells, 1., TraversalOption::c08, DataLayoutOption::aos,
                            Newton3Option::enabled);
  autopas::FeatureVector f3(ContainerOption::linkedCells, 1., TraversalOption::c08, DataLayoutOption::soa,
                            Newton3Option::enabled);
  autopas::FeatureVector f4(ContainerOption::linkedCells, 1., TraversalOption::c08, DataLayoutOption::soa,
                            Newton3Option::disabled);

  EXPECT_EQ(static_cast<Eigen::VectorXd>(f1 - f1).squaredNorm(), 0);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f2 - f2).squaredNorm(), 0);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f3 - f3).squaredNorm(), 0);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f4 - f4).squaredNorm(), 0);

  EXPECT_EQ(static_cast<Eigen::VectorXd>(f1 - f2).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f2 - f3).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f3 - f4).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f4 - f3).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f3 - f2).squaredNorm(), 1);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f2 - f1).squaredNorm(), 1);

  EXPECT_EQ(static_cast<Eigen::VectorXd>(f1 - f3).squaredNorm(), 2);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f2 - f4).squaredNorm(), 2);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f4 - f2).squaredNorm(), 2);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f3 - f1).squaredNorm(), 2);

  EXPECT_EQ(static_cast<Eigen::VectorXd>(f1 - f4).squaredNorm(), 3);
  EXPECT_EQ(static_cast<Eigen::VectorXd>(f4 - f1).squaredNorm(), 3);
}

TEST(GaussianProcessTest, noEvidence) {
  Random rng(32);

  double epsilon = 0.05;  // allowed error for tests

  double sigma = 0.001;  // noise
  GaussianProcess<Eigen::VectorXd> gp(1, sigma, rng);

  Eigen::VectorXd f1(1);
  f1 << 0.;

  // predict without information -> should return default values
  EXPECT_NEAR(gp.predictMean(f1), 0., epsilon);
  EXPECT_NEAR(gp.predictVar(f1), 1., epsilon);
}

TEST(GaussianProcessTest, oneEvidence) {
  Random rng(32);

  double epsilon = 0.05;
  double sigma = 0.001;
  GaussianProcess<Eigen::VectorXd> gp(1, sigma, rng);

  Eigen::VectorXd f1(1);
  f1 << 0.;
  double out1 = 42.;

  Eigen::VectorXd f2(1);
  f2 << 1000.;

  gp.addEvidence(f1, out1);

  // predicting point same as evidence -> should expect same output as evidence
  EXPECT_NEAR(gp.predictMean(f1), out1, epsilon);
  EXPECT_NEAR(gp.predictVar(f1), 0., epsilon);

  // not enough information, so other inputs lead to same mean
  EXPECT_NEAR(gp.predictMean(f2), out1, epsilon);
}

TEST(GaussianProcessTest, twoEvidence) {
  Random rng;

  double epsilon = 0.05;  // allowed error for tests
  double sigma = 0.001;   // noise
  GaussianProcess<Eigen::VectorXd> gp(1, sigma, rng);

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
  // should expect same output as evidence because great distance between inputs
  EXPECT_NEAR(gp.predictMean(f1), out1, epsilon);
  EXPECT_NEAR(gp.predictVar(f1), 0., epsilon);

  EXPECT_NEAR(gp.predictMean(f2), out2, epsilon);
  EXPECT_NEAR(gp.predictVar(f2), 0., epsilon);
}

TEST(GaussianProcessTest, clear) {
  Random rng;

  double epsilon = 0.05;  // allowed error for tests
  double sigma = 0.001;   // noise
  GaussianProcess<Eigen::VectorXd> gp(1, sigma, rng);

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
  EXPECT_NEAR(gp.predictMean(f1), 0., epsilon);
  EXPECT_NEAR(gp.predictVar(f1), 1., epsilon);

  EXPECT_NEAR(gp.predictMean(f2), 0., epsilon);
  EXPECT_NEAR(gp.predictVar(f2), 1., epsilon);

  // test new evidence at deleted point
  out2 = 10.;
  gp.addEvidence(f2, out2);

  EXPECT_NEAR(gp.predictMean(f1), out2, epsilon);

  EXPECT_NEAR(gp.predictMean(f2), out2, epsilon);
  EXPECT_NEAR(gp.predictVar(f2), 0., epsilon);
}

TEST(GaussianProcessTest, sine) {
  Random rng(42);

  // gp should try to approximate the sine as blackbox function
  auto functor = [](double input) { return std::sin(input); };
  double epsilon = 0.05;           // allowed error
  unsigned numEvidence = 9;        // number of evidence to provide
  unsigned numPredictions = 1000;  // number of predictions to make
  double domainStart = 0.;         // start of tested domain
  double domainEnd = 2 * M_PI;     // end of tested domain

  GaussianProcess<Eigen::VectorXd> gp(1, 0.001, rng);

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

    EXPECT_NEAR(gp.predictMean(f), output, epsilon);
  }
}

TEST(GaussianProcessTest, 2dMax) {
  Random rng(72);  // random generator

  // try to find the max of -(i1 + 1)^2 - (i2 - 1)^2
  auto functor = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 1, 2); };
  double epsilon = 0.1;  // allowed error
  std::vector<NumberInterval<double>> domain{NumberInterval<double>(-2, 2),
                                             NumberInterval<double>(-2, 2)};  // domain of function

  // max of function
  Eigen::VectorXd max(2);
  max << -1, 1;
  unsigned numEvidence = 10;      // number of samples allowed to make
  unsigned lhsNumSamples = 1000;  // number of sample to find max of acquisition function
  AcquisitionFunctionOption af = AcquisitionFunctionOption::upperConfidenceBound;

  GaussianProcess<Eigen::VectorXd> gp(2, 0.001, rng);

  // add first evidence
  Eigen::VectorXd first(2);
  first << 0, 0;
  gp.addEvidence(first, functor(0, 0));
  std::cout << "Origin: " << functor(0, 0) << std::endl;

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

    // print acquisition and mean map
    int xChunks = 20;
    int yChunks = 20;
    double xSpace = 4. / (xChunks - 1);
    double ySpace = 4. / (yChunks - 1);
    for (int y = yChunks - 1; y >= 0; --y) {
      // acquisition row
      for (int x = 0; x < xChunks; ++x) {
        Eigen::VectorXd sample(2);
        sample << (x * xSpace - 2), (y * ySpace - 2);
        double val = gp.calcAcquisition(af, sample);
        int color = static_cast<int>(val * 2.5 + 250);
        color = std::clamp(color, 232, 255);

        std::cout << "\033[48;5;" << color << "m  ";
      }
      std::cout << "\033[0m ";

      // mean row
      for (int x = 0; x < xChunks; ++x) {
        Eigen::VectorXd sample(2);
        sample << (x * xSpace - 2), (y * ySpace - 2);
        double val = gp.predictMean(sample);
        int color = static_cast<int>(val * 2.5 + 250);
        color = std::clamp(color, 232, 255);

        std::cout << "\033[48;5;" << color << "m  ";
      }
      std::cout << "\033[0m" << std::endl;
    }
    std::cout << "Acq max: " << std::endl << am << std::endl;
    std::cout << "Got: " << amOut << std::endl;

    gp.addEvidence(am, amOut);
  }

  // get max
  Eigen::VectorXd am = gp.getEvidenceMax();

  // check if max is near real max
  double predMax = functor(am[0], am[1]);
  double realMax = functor(max[0], max[1]);
  EXPECT_NEAR(predMax, realMax, epsilon);
}

TEST(GaussianProcessTest, 2dMin) {
  Random rng(73);  // random generator

  // try to find the min of (i1 - 1)^2 + (i2 - 1)^2
  auto functor = [](double i1, double i2) { return std::pow(i1 - 1, 2) + std::pow(i2 - 1, 2); };
  double epsilon = 0.1;  // allowed error
  std::vector<NumberInterval<double>> domain{NumberInterval<double>(-2, 2),
                                             NumberInterval<double>(-2, 2)};  // domain of function

  // min of function
  Eigen::VectorXd min(2);
  min << 1, 1;
  unsigned numEvidence = 10;      // number of samples allowed to make
  unsigned lhsNumSamples = 1000;  // number of sample to find min of acquisition function
  AcquisitionFunctionOption af = AcquisitionFunctionOption::lowerConfidenceBound;

  GaussianProcess<Eigen::VectorXd> gp(2, 0.001, rng);

  // add first evidence
  Eigen::VectorXd first(2);
  first << 0, 0;
  gp.addEvidence(first, functor(0, 0));
  std::cout << "Origin: " << functor(0, 0) << std::endl;

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

    // print acquisition map
    int xChunks = 20;
    int yChunks = 20;
    double xSpace = 4. / (xChunks - 1);
    double ySpace = 4. / (yChunks - 1);
    for (int y = yChunks - 1; y >= 0; --y) {
      for (int x = 0; x < xChunks; ++x) {
        Eigen::VectorXd sample(2);
        sample << (x * xSpace - 2), (y * ySpace - 2);
        double val = gp.calcAcquisition(af, sample);
        int color = static_cast<int>(val * 2 + 232);
        color = std::clamp(color, 232, 255);

        std::cout << "\033[48;5;" << color << "m  ";
      }
      std::cout << "\033[0m ";

      // mean row
      for (int x = 0; x < xChunks; ++x) {
        Eigen::VectorXd sample(2);
        sample << (x * xSpace - 2), (y * ySpace - 2);
        double val = gp.predictMean(sample);
        int color = static_cast<int>(val * 2 + 232);
        color = std::clamp(color, 232, 255);

        std::cout << "\033[48;5;" << color << "m  ";
      }
      std::cout << "\033[0m" << std::endl;
    }
    std::cout << "Acq min: " << std::endl << am << std::endl;
    std::cout << "Got: " << amOut << std::endl;

    gp.addEvidence(am, amOut);
  }

  // get min
  Eigen::VectorXd am = gp.getEvidenceMin();

  // check if min is near real min
  double predMin = functor(am[0], am[1]);
  double realMin = functor(min[0], min[1]);
  EXPECT_NEAR(predMin, realMin, epsilon);
}

TEST(GaussianProcessTest, 2dMinGrid) {
  Random rng(73);  // random generator

  // try to find the min of (i1 - 1)^2 + (i2 - 1)^2
  auto functor = [](double i1, double i2) { return std::pow(i1 - 1, 2) + std::pow(i2 - 1, 2); };
  double epsilon = 0.1;  // allowed error

  // domain of function
  int domHalf = 10;
  std::set<double> domSet;
  for (int i = -domHalf; i <= domHalf; ++i) {
    domSet.insert(i * 2. / domHalf);
  }
  std::vector<NumberSetFinite<double>> domain{NumberSetFinite<double>(domSet), NumberSetFinite<double>(domSet)};

  // min of function
  Eigen::VectorXd min(2);
  min << 1, 1;
  unsigned numEvidence = 10;      // number of samples allowed to make
  unsigned lhsNumSamples = 1000;  // number of sample to find min of acquisition function
  AcquisitionFunctionOption af = AcquisitionFunctionOption::lowerConfidenceBound;

  GaussianProcess<Eigen::VectorXd> gp(2, 0.001, rng);

  // add first evidence
  Eigen::VectorXd first(2);
  first << 0, 0;
  gp.addEvidence(first, functor(0, 0));
  std::cout << "Origin: " << functor(0, 0) << std::endl;

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

    // print acquisition and mean map
    int xChunks = 20;
    int yChunks = 20;
    double xSpace = 4. / (xChunks - 1);
    double ySpace = 4. / (yChunks - 1);
    for (int y = yChunks - 1; y >= 0; --y) {
      // acquisition row
      for (int x = 0; x < xChunks; ++x) {
        Eigen::VectorXd sample(2);
        sample << (x * xSpace - 2), (y * ySpace - 2);
        double val = gp.calcAcquisition(af, sample);
        int color = static_cast<int>(val * 2 + 232);
        color = std::clamp(color, 232, 255);

        std::cout << "\033[48;5;" << color << "m  ";
      }
      std::cout << "\033[0m ";

      // mean row
      for (int x = 0; x < xChunks; ++x) {
        Eigen::VectorXd sample(2);
        sample << (x * xSpace - 2), (y * ySpace - 2);
        double val = gp.predictMean(sample);
        int color = static_cast<int>(val * 2 + 232);
        color = std::clamp(color, 232, 255);

        std::cout << "\033[48;5;" << color << "m  ";
      }
      std::cout << "\033[0m" << std::endl;
    }
    std::cout << "Acq min: " << std::endl << am << std::endl;
    std::cout << "Got: " << amOut << std::endl;

    gp.addEvidence(am, amOut);
  }

  // get min
  Eigen::VectorXd am = gp.getEvidenceMin();

  // check if min is near real min
  double predMin = functor(am[0], am[1]);
  double realMin = functor(min[0], min[1]);
  EXPECT_NEAR(predMin, realMin, epsilon);
}

TEST(GaussianProcessTest, 2dMinGridBig) {
  Random rng(73);  // random generator

  // functor to find min of
  auto functor = [](double i1, double i2) {
    return std::pow((std::abs(i1 - 1) + 1), 5) + std::pow((std::abs(i2 - 1) + 1), 5);
  };
  double epsilon = 10;  // allowed error

  // domain of function
  int domHalf = 10;
  std::set<double> domSet;
  for (int i = -domHalf; i <= domHalf; ++i) {
    domSet.insert(i * 2. / domHalf);
  }
  std::vector<NumberSetFinite<double>> domain{NumberSetFinite<double>(domSet), NumberSetFinite<double>(domSet)};

  // min of function
  Eigen::VectorXd min(2);
  min << 1, 1;
  unsigned numEvidence = 10;      // number of samples allowed to make
  unsigned lhsNumSamples = 1000;  // number of sample to find min of acquisition function
  AcquisitionFunctionOption af = AcquisitionFunctionOption::lowerConfidenceBound;

  GaussianProcess<Eigen::VectorXd> gp(2, 0.001, rng);

  // add first evidence
  Eigen::VectorXd first(2);
  first << 0, 0;
  gp.addEvidence(first, functor(0, 0));
  std::cout << "Origin: " << functor(0, 0) << std::endl;

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

    // print acquisition map and mean map
    int xChunks = 20;
    int yChunks = 20;
    double xSpace = 4. / (xChunks - 1);
    double ySpace = 4. / (yChunks - 1);
    for (int y = yChunks - 1; y >= 0; --y) {
      // acqusition row
      for (int x = 0; x < xChunks; ++x) {
        Eigen::VectorXd sample(2);
        sample << (x * xSpace - 2), (y * ySpace - 2);
        double val = gp.calcAcquisition(af, sample);
        int color = static_cast<int>(val / 100 + 232);
        color = std::clamp(color, 232, 255);

        std::cout << "\033[48;5;" << color << "m  ";
      }
      std::cout << "\033[0m ";

      // mean row
      for (int x = 0; x < xChunks; ++x) {
        Eigen::VectorXd sample(2);
        sample << (x * xSpace - 2), (y * ySpace - 2);
        double val = gp.predictMean(sample);
        int color = static_cast<int>(val / 100 + 232);
        color = std::clamp(color, 232, 255);

        std::cout << "\033[48;5;" << color << "m  ";
      }
      std::cout << "\033[0m" << std::endl;
    }
    std::cout << "Acq min: " << std::endl << am << std::endl;
    std::cout << "Got: " << amOut << std::endl;

    gp.addEvidence(am, amOut);
  }

  // get min
  Eigen::VectorXd am = gp.getEvidenceMin();

  // check if min is near real min
  double predMin = functor(am[0], am[1]);
  double realMin = functor(min[0], min[1]);
  EXPECT_NEAR(predMin, realMin, epsilon);
}

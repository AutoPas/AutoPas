/**
 * @file GaussianClusterTest.cpp
 * @author Jan Nguyen
 * @date 21.04.20
 */

#include "GaussianClusterTest.h"

using namespace autopas;

TEST_F(GaussianClusterTest, wrongDimension) {
  Random rng(32);
  GaussianCluster cluster({2, 2}, 2, 0.001, rng);

  Eigen::VectorXi fd = Eigen::VectorXi::Zero(2);
  Eigen::VectorXd f1 = Eigen::VectorXd::Ones(1);
  Eigen::VectorXd f2 = Eigen::VectorXd::Zero(3);

  Eigen::VectorXi g1 = Eigen::VectorXi::Zero(1);
  Eigen::VectorXi g2 = Eigen::VectorXi::Zero(3);
  Eigen::VectorXd gc = Eigen::VectorXd::Ones(2);

  EXPECT_THROW(cluster.addEvidence(fd, f1, 0), utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(cluster.addEvidence(fd, f2, 1), utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(cluster.addEvidence(g1, gc, 1), utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(cluster.addEvidence(g2, gc, 1), utils::ExceptionHandler::AutoPasException);
}

TEST_F(GaussianClusterTest, maxNoNeighbours) {
  auto functor1 = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 0, 2) - 2; };
  auto functor2 = [](double i1, double i2) { return -std::pow(i1 + 0, 2) - std::pow(i2 - 1, 2) - 2; };
  auto functor3 = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 1, 2); };

  auto neighboursFun = [](Eigen::VectorXi) -> std::vector<Eigen::VectorXi> { return {}; };

  std::pair domain{NumberInterval<double>(-2, 2), NumberInterval<double>(-2, 2)};
  constexpr double maxError = 1.;

  // max of function
  int maxDiscrete = 2;
  auto maxContinuous = autopas::utils::Math::makeVectorXd({-1, 1});

  test2DFunctions({functor1, functor2, functor3}, neighboursFun, maxDiscrete, maxContinuous, maxError, domain,
                  AcquisitionFunctionOption::upperConfidenceBound, false);
}

TEST_F(GaussianClusterTest, maxAllNeighbours) {
  auto functor1 = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 1, 2) - 2; };
  auto functor2 = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 1, 2) - 1; };
  auto functor3 = [](double i1, double i2) { return -std::pow(i1 + 1, 2) - std::pow(i2 - 1, 2); };
  std::vector<std::function<double(double, double)>> functors = {functor1, functor2, functor3};

  std::vector<Eigen::VectorXi> allDiscrete;
  for (int i = 0; i < functors.size(); ++i) {
    Eigen::VectorXi vec(1);
    vec << i;
    allDiscrete.emplace_back(std::move(vec));
  }

  // return all but target as neighbours
  auto neighboursFun = [&allDiscrete](Eigen::VectorXi target) -> std::vector<Eigen::VectorXi> {
    std::vector<Eigen::VectorXi> result;
    for (const auto &vec : allDiscrete) {
      if (vec != target) {
        result.emplace_back(vec);
      }
    }
    return result;
  };

  std::pair domain{NumberInterval<double>(-2, 2), NumberInterval<double>(-2, 2)};
  constexpr double maxError = 1.;

  // max of function
  int maxDiscrete = 2;
  auto maxContinuous = autopas::utils::Math::makeVectorXd({-1, 1});

  test2DFunctions(functors, neighboursFun, maxDiscrete, maxContinuous, maxError, domain,
                  AcquisitionFunctionOption::upperConfidenceBound, false);
}

void GaussianClusterTest::printMaps(int xChunks, int yChunks, const autopas::NumberSet<double> &domainX,
                                    const autopas::NumberSet<double> &domainY, const autopas::GaussianCluster &gc,
                                    std::function<std::vector<Eigen::VectorXi>(Eigen::VectorXi)> neighboursFun,
                                    autopas::AcquisitionFunctionOption af) {
  // get distance between chunks
  double xSpace = (domainX.getMax() - domainX.getMin()) / (xChunks - 1);
  double ySpace = (domainY.getMax() - domainY.getMin()) / (yChunks - 1);

  size_t numFunctors = gc.getDimensions()[0];

  // precalculate acqMaps
  std::vector<std::vector<std::vector<double>>> acqMaps;
  acqMaps.reserve(numFunctors);
  double acqMin = std::numeric_limits<double>::max();
  double acqMax = std::numeric_limits<double>::lowest();
  for (size_t i = 0; i < numFunctors; ++i) {
    std::vector<std::vector<double>> acqMap;
    Eigen::VectorXi sampleDiscrete(1);
    sampleDiscrete << i;

    for (int y = 0; y < xChunks; ++y) {
      // calculate a row
      std::vector<double> row;
      for (int x = 0; x < xChunks; ++x) {
        // calculate value of chunk
        Eigen::Vector2d sampleContinuous(x * xSpace + domainX.getMin(), (y * ySpace + domainY.getMin()));
        double val = gc.calcAcquisition(af, sampleDiscrete, sampleContinuous, neighboursFun(sampleDiscrete));

        row.push_back(val);

        // keep track min and max value
        acqMin = std::min(val, acqMin);
        acqMax = std::max(val, acqMax);
      }
      acqMap.push_back(std::move(row));
    }
    acqMaps.push_back(std::move(acqMap));
  }

  // get scaling such that acqMax=1 and acqMin=0
  double acqScale = 1 / (acqMax - acqMin);

  // print map
  for (int y = yChunks - 1; y >= 0; --y) {
    // acquisiton row
    for (size_t i = 0; i < numFunctors; ++i) {
      for (int x = 0; x < xChunks; ++x) {
        // map value between 0 to 1 and square for clearer differences of high values
        double val = std::pow((acqMaps[i][y][x] - acqMin) * acqScale, 2);

        // map value to color
        int color = static_cast<int>(232 + 23 * val);
        color = std::clamp(color, 232, 255);

        // print two spaces of that color
        std::cout << "\033[48;5;" << color << "m  ";
      }
      // reset color, print a space
      std::cout << "\033[0m ";
    }
    std::cout << std::endl;
  }
}

void GaussianClusterTest::printEvidence(Eigen::VectorXi vecDiscrete, Eigen::VectorXd vecContinuous, double out,
                                        size_t evidenceNum) {
  std::cout << "Evidence " << evidenceNum << ": " << vecDiscrete[0] << "," << vecContinuous[0] << ","
            << vecContinuous[1] << " -> " << out << std::endl;
}

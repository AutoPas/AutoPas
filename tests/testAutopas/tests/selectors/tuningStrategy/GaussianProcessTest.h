/**
 * @file GaussianProcessTest.h
 * @author Jan Nguyen
 * @date 12.06.19
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "Eigen/Dense"
#include "autopas/selectors/FeatureVector.h"
#include "autopas/selectors/tuningStrategy/GaussianProcess.h"
#include "autopas/utils/NumberSet.h"

class GaussianProcessTest : public AutoPasTestBase {};

/**
 * Print a xChunks x yChunks mean and aquisition map.
 * @param xChunks width of maps
 * @param yChunks height of maps
 * @param domainX x domain of GaussianProcess
 * @param domainY y domain of GaussianProcess
 * @param gp GaussianProcess to calculate mean and acquisiton
 * @param af acquisition function
 * @param colorFactor Multiplies the mean by this factor to gain the color
 */
void printMap(int xChunks, int yChunks, const autopas::NumberSet<double> &domainX,
              const autopas::NumberSet<double> &domainY, const autopas::GaussianProcess<Eigen::VectorXd> &gp,
              autopas::AcquisitionFunctionOption af, double colorFactor) {
  // get distance between chunks
  double xSpace = (domainX.getMax() - domainX.getMin()) / (xChunks - 1);
  double ySpace = (domainY.getMax() - domainY.getMin()) / (yChunks - 1);

  // precalculate acqMap
  std::vector<std::vector<double>> acqMap;
  double acqMin = std::numeric_limits<double>::max();
  double acqMax = std::numeric_limits<double>::min();
  for (int y = 0; y < xChunks; ++y) {
    // calculate a row
    std::vector<double> row;
    for (int x = 0; x < xChunks; ++x) {
      // calculate value of chunk
      Eigen::Vector2d sample(x * xSpace + domainX.getMin(), (y * ySpace + domainY.getMin()));
      double val = gp.calcAcquisition(af, sample);

      // negate for special case lcb
      if (af == autopas::AcquisitionFunctionOption::lowerConfidenceBound) val = -val;

      row.push_back(val);

      // keep track min and max value
      acqMin = std::min(val, acqMin);
      acqMax = std::max(val, acqMax);
    }

    acqMap.push_back(row);
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
      int color = static_cast<int>(255 - val);
      color = std::clamp(color, 232, 255);

      // print two spaces of that color
      std::cout << "\033[48;5;" << color << "m  ";
    }
    // reset color, print a space
    std::cout << "\033[0m" << std::endl;
  }
}

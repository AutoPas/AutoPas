/**
 * @file CubeGauss.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Objects.h"

class CubeGauss : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * @param numParticles
   * @param boxLength
   * @param distributionMean
   * @param distributionStdDev
   * @param bottomLeftCorner
   */
  CubeGauss(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass,
            size_t numParticles, const std::array<double, 3> &boxLength, double distributionMean,
            double distributionStdDev, const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        numParticles(numParticles),
        boxLength(boxLength),
        distributionMean(distributionMean),
        distributionStdDev(distributionStdDev),
        bottomLeftCorner(bottomLeftCorner){}

            /**
             * Getter total number of Particles of Object
             * @return numParticles
             */
            [[nodiscard]] size_t getParticlesTotal() const override {
    return numParticles;
  }
  /**
   * Getter for distribution mean
   * @return distributionMean
   */
  [[nodiscard]] double getDistributionMean() const { return distributionMean; }
      /**
       * Getter for distributionStdDev
       * @return distributionStdDev
       */
      [[nodiscard]] double getDistributionStdDev() const {
    return distributionStdDev;
  }

  /**
   * Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   */
  const std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }
  /**
   * Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   */
  const std::array<double, 3> getBoxMax() const override {
    return autopas::ArrayMath::add(bottomLeftCorner, boxLength);
  }

  /**
   * Prints the Configuration of the current Object
   */
  std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "Distribution-Mean"
           << ":  " << distributionMean << std::endl;
    output << std::setw(_valueOffset) << std::left << "Distribution-StdDev"
           << ":  " << distributionStdDev << std::endl;
    output << std::setw(_valueOffset) << std::left << "NumberOfParticles"
           << ":  " << numParticles << std::endl;
    output << std::setw(_valueOffset) << std::left << "BoxLength"
           << ":  " << autopas::ArrayUtils::to_string(boxLength) << std::endl;
    output << Object::to_string();
    return output.str();
  }

 private:
  size_t numParticles;
  std::array<double, 3> boxLength;
  double distributionMean;
  double distributionStdDev;
  std::array<double, 3> bottomLeftCorner;
};
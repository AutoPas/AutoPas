/**
 * @file CubeGauss.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Objects.h"
#include "autopasTools/generators/GaussianGenerator.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Class describing an cuboid object filled with gaussian randomly distributed particles.
 */
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
            size_t numParticles, const std::array<double, 3> &boxLength, const std::array<double, 3> &distributionMean,
            const std::array<double, 3> &distributionStdDev, const std::array<double, 3> &bottomLeftCorner)
      : Object(velocity, typeId, epsilon, sigma, mass),
        numParticles(numParticles),
        boxLength(boxLength),
        distributionMean(distributionMean),
        distributionStdDev(distributionStdDev),
        bottomLeftCorner(bottomLeftCorner) {}

  /**
   * Getter for distribution mean
   * @return distributionMean
   */
  [[nodiscard]] const std::array<double, 3> &getDistributionMean() const { return distributionMean; }

  /**
   * Getter for distributionStdDev
   * @return distributionStdDev
   */
  [[nodiscard]] const std::array<double, 3> &getDistributionStdDev() const { return distributionStdDev; }

  [[nodiscard]] size_t getParticlesTotal() const override { return numParticles; }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength);
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "distribution-mean"
           << ":  " << autopas::utils::ArrayUtils::to_string(distributionMean) << std::endl;
    output << std::setw(_valueOffset) << std::left << "distribution-stddeviation"
           << ":  " << autopas::utils::ArrayUtils::to_string(distributionStdDev) << std::endl;
    output << std::setw(_valueOffset) << std::left << "numberOfParticles"
           << ":  " << numParticles << std::endl;
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(boxLength) << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

  void generate(autopas::AutoPas<ParticleType> &autopas) const override {
    ParticleType dummyParticle = getDummyParticle(autopas);
    autopasTools::generators::GaussianGenerator::fillWithParticles(autopas, getBoxMin(), getBoxMax(), numParticles,
                                                                   dummyParticle, distributionMean, distributionStdDev);
  }

 private:
  size_t numParticles;
  std::array<double, 3> boxLength;
  std::array<double, 3> distributionMean;
  std::array<double, 3> distributionStdDev;
  std::array<double, 3> bottomLeftCorner;
};
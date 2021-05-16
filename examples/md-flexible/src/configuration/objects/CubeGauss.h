/**
 * @file CubeGauss.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "src/ParticleAttributes.h"
#include "autopas/utils/ArrayMath.h"
#include "Object.h"

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
        _numParticles(numParticles),
        boxLength(boxLength),
        _distributionMean(distributionMean),
        _distributionStdDev(distributionStdDev),
        bottomLeftCorner(bottomLeftCorner) {}

  /**
   * Getter for distribution mean
   * @return distributionMean
   */
  [[nodiscard]] const std::array<double, 3> &getDistributionMean() const { return _distributionMean; }

  /**
   * Getter for distributionStdDev
   * @return distributionStdDev
   */
  [[nodiscard]] const std::array<double, 3> &getDistributionStdDev() const { return _distributionStdDev; }

  [[nodiscard]] size_t getParticlesTotal() const override { return _numParticles; }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return bottomLeftCorner; }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return autopas::utils::ArrayMath::add(bottomLeftCorner, boxLength);
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "distribution-mean"
           << ":  " << autopas::utils::ArrayUtils::to_string(_distributionMean) << std::endl;
    output << std::setw(_valueOffset) << std::left << "distribution-stddeviation"
           << ":  " << autopas::utils::ArrayUtils::to_string(_distributionStdDev) << std::endl;
    output << std::setw(_valueOffset) << std::left << "numberOfParticles"
           << ":  " << _numParticles << std::endl;
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(boxLength) << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

  void generate(std::vector<ParticleAttributes> &particles) const override {
    ParticleAttributes particle = getDummyParticle(particles.size());
  	std::default_random_engine generator(42);
  	std::array<std::normal_distribution<double>, 3> distributions = {
      std::normal_distribution<double>{_distributionMean[0], _distributionStdDev[0]},
      std::normal_distribution<double>{_distributionMean[1], _distributionStdDev[1]},
      std::normal_distribution<double>{_distributionMean[2], _distributionStdDev[2]}
		};
  	for (int i = 0; i < _numParticles; ++i) {
			particle.id++;
			particle.positionX = distributions[0](generator);
			particle.positionY = distributions[1](generator);
			particle.positionZ = distributions[2](generator);
			particles.push_back(particle);
  	}
  }

 private:
  size_t _numParticles;
  std::array<double, 3> boxLength;
  std::array<double, 3> _distributionMean;
  std::array<double, 3> _distributionStdDev;
  std::array<double, 3> bottomLeftCorner;
};

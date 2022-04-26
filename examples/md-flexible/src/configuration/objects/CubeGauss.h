/**
 * @file CubeGauss.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Object.h"
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
        _numParticles(numParticles),
        _boxLength(boxLength),
        _distributionMean(distributionMean),
        _distributionStdDev(distributionStdDev),
        _bottomLeftCorner(bottomLeftCorner) {}

  /**
   * Getter for distribution mean.
   * @return distributionMean
   */
  [[nodiscard]] const std::array<double, 3> &getDistributionMean() const { return _distributionMean; }

  /**
   * Getter for distributionStdDev.
   * @return distributionStdDev
   */
  [[nodiscard]] const std::array<double, 3> &getDistributionStdDev() const { return _distributionStdDev; }

  /**
   * Returns the number of particles generated in this CubeGauss object.
   * @return number of particles
   */
  [[nodiscard]] size_t getParticlesTotal() const override { return _numParticles; }

  /**
   * Returns the bottom left corner of the cube.
   * @return bottom left corner.
   */
  [[nodiscard]] std::array<double, 3> getBoxMin() const override { return _bottomLeftCorner; }

  /**
   * Returns the top right corner of the cube.
   * @return top right corner.
   */
  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    return autopas::utils::ArrayMath::add(_bottomLeftCorner, _boxLength);
  }

  /**
   * Converts the cube gauss object to a human readable string.
   * @return human readable cube gauss object.
   */
  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "distribution-mean"
           << ":  " << autopas::utils::ArrayUtils::to_string(_distributionMean) << std::endl;
    output << std::setw(_valueOffset) << std::left << "distribution-stddeviation"
           << ":  " << autopas::utils::ArrayUtils::to_string(_distributionStdDev) << std::endl;
    output << std::setw(_valueOffset) << std::left << "numberOfParticles"
           << ":  " << _numParticles << std::endl;
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(_boxLength) << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << std::endl;
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates the particles based on the configuration of the cube gauss object provided in the yaml file.
   * @param particles The container where the new particles will be stored.
   */
  void generate(std::vector<MulticenteredParticleType> &particles) const override {
    MulticenteredParticleType particle = getDummyParticle(particles.size());

    std::default_random_engine generator(42);
    std::array<std::normal_distribution<double>, 3> distributions = {
        std::normal_distribution<double>{_distributionMean[0], _distributionStdDev[0]},
        std::normal_distribution<double>{_distributionMean[1], _distributionStdDev[1]},
        std::normal_distribution<double>{_distributionMean[2], _distributionStdDev[2]}};

    for (int i = 0; i < _numParticles; ++i) {
      particle.setR({
        _bottomLeftCorner[0] + distributions[0](generator),
        _bottomLeftCorner[1] + distributions[1](generator),
        _bottomLeftCorner[2] + distributions[2](generator)
      });

      particles.push_back(particle);
      particle.setID(particle.getID() + 1);
    }
  }

 private:
  /**
   * The number of particles which will be generated.
   */
  size_t _numParticles;

  /**
   * The length of the box in each dimension.
   */
  std::array<double, 3> _boxLength;

  /**
   * The mean value for the gaussian distribution of the particles.
   */
  std::array<double, 3> _distributionMean;

  /**
   * The standard deviation value for the gaussian distribution of the particles.
   */
  std::array<double, 3> _distributionStdDev;

  /**
   * The coordinates of the bottom left corner of the cube object.
   */
  std::array<double, 3> _bottomLeftCorner;
};

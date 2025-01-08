/**
 * @file CubeGauss.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

#include "Object.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/PseudoContainer.h"
#include "autopasTools/generators/GaussianGenerator.h"

/**
 * Class describing an cuboid object filled with gaussian randomly distributed particles.
 */
class CubeGauss : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param numParticles
   * @param boxLength
   * @param distributionMean
   * @param distributionStdDev
   * @param bottomLeftCorner
   */
  CubeGauss(const std::array<CalcPrecision, 3> &velocity, unsigned long typeId, size_t numParticles,
            const std::array<CalcPrecision, 3> &boxLength, const std::array<CalcPrecision, 3> &distributionMean,
            const std::array<CalcPrecision, 3> &distributionStdDev,
            const std::array<CalcPrecision, 3> &bottomLeftCorner)
      : Object(velocity, typeId),
        _numParticles(numParticles),
        _boxLength(boxLength),
        _distributionMean(distributionMean),
        _distributionStdDev(distributionStdDev),
        _bottomLeftCorner(bottomLeftCorner) {}

  /**
   * Getter for distribution mean.
   * @return distributionMean
   */
  [[nodiscard]] const std::array<CalcPrecision, 3> &getDistributionMean() const { return _distributionMean; }

  /**
   * Getter for distributionStdDev.
   * @return distributionStdDev
   */
  [[nodiscard]] const std::array<CalcPrecision, 3> &getDistributionStdDev() const { return _distributionStdDev; }

  /**
   * Returns the number of particles generated in this CubeGauss object.
   * @return number of particles
   */
  [[nodiscard]] size_t getParticlesTotal() const override { return _numParticles; }

  /**
   * Returns the bottom left front corner of the cube.
   * @return bottom left front corner.
   */
  [[nodiscard]] std::array<CalcPrecision, 3> getBoxMin() const override { return _bottomLeftCorner; }

  /**
   * Returns the top right back corner of the cube.
   * @return top right back corner.
   */
  [[nodiscard]] std::array<CalcPrecision, 3> getBoxMax() const override {
    using namespace autopas::utils::ArrayMath::literals;
    return _bottomLeftCorner + _boxLength;
  }

  /**
   * Converts the cube gauss object to a human readable string.
   * @return human readable cube gauss object.
   */
  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;

    output << std::setw(_valueOffset) << std::left << "distribution-mean"
           << ":  " << autopas::utils::ArrayUtils::to_string(_distributionMean) << "\n";
    output << std::setw(_valueOffset) << std::left << "distribution-stddeviation"
           << ":  " << autopas::utils::ArrayUtils::to_string(_distributionStdDev) << "\n";
    output << std::setw(_valueOffset) << std::left << "numberOfParticles"
           << ":  " << _numParticles << "\n";
    output << std::setw(_valueOffset) << std::left << "box-length"
           << ":  " << autopas::utils::ArrayUtils::to_string(_boxLength) << "\n";
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(_bottomLeftCorner) << "\n";
    output << Object::to_string();
    return output.str();
  }

  /**
   * Generates the particles based on the configuration of the cube gauss object provided in the yaml file.
   * @param particles The container where the new particles will be stored.
   */
  void generate(std::vector<ParticleType> &particles) const override {
    // Wrapper so that std::vector can be used as an AutoPas::ParticleContainer
    auto particlesWrapper = autopasTools::PseudoContainer(particles);

    using namespace autopas::utils::ArrayMath::literals;
    const auto boxMax = _bottomLeftCorner + _boxLength;

    // dummy particle used as a template with id of the first newly generated one
    const ParticleType dummyParticle = getDummyParticle(particles.size());

    autopasTools::generators::GaussianGenerator::fillWithParticles(particlesWrapper, _bottomLeftCorner, boxMax,
                                                                   _numParticles, dummyParticle, _distributionMean,
                                                                   _distributionStdDev);
  }

 private:
  /**
   * The number of particles which will be generated.
   */
  size_t _numParticles;

  /**
   * The length of the box in each dimension.
   */
  std::array<CalcPrecision, 3> _boxLength;

  /**
   * The mean value for the gaussian distribution of the particles.
   */
  std::array<CalcPrecision, 3> _distributionMean;

  /**
   * The standard deviation value for the gaussian distribution of the particles.
   */
  std::array<CalcPrecision, 3> _distributionStdDev;

  /**
   * The coordinates of the bottom left front corner of the cube object.
   */
  std::array<CalcPrecision, 3> _bottomLeftCorner;
};

/**
 * @file CubeUniform.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once

class CubeUniform : public Object {
 public:
  /**
   * constructor for CubeUniform that is created in YamlParser and then included into the Simulation via the Generator
   * class
   * @param numParticles
   * @param boxLength
   * @param velocity_arg
   * @param bottomLeftCorner
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   */
  CubeUniform(size_t numParticles, const std::array<double, 3> &boxLength,
              const std::array<double, 3> &bottomLeftCorner, const std::array<double, 3> &velocity_arg,
              const unsigned long &typeId_arg, const double &epsilon_arg, const double &sigma_arg,
              const double &mass_arg)
      : numParticles(numParticles), boxLength(boxLength), bottomLeftCorner(bottomLeftCorner) {
    velocity = velocity_arg;
    typeId = typeId_arg;
    epsilon = epsilon_arg;
    sigma = sigma_arg;
    mass = mass_arg;
  }

  /**
   * Getter for total number of Particles in Object
   * @return numParticles
   */
  [[nodiscard]] size_t getParticlesTotal() const override { return numParticles; }

  /**
   * Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   */
  const std::array<double, 3> getBoxMin() const override {
    return bottomLeftCorner;
  }
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
  void printConfig() override {
    using namespace std;

    cout << std::setw(valueOffset) << left << "Center"
         << ":  " << autopas::ArrayUtils::to_string(bottomLeftCorner) << endl;
    cout << std::setw(valueOffset) << left << "NumberOfParticles"
         << ":  " << numParticles << endl;
    cout << std::setw(valueOffset) << left << "BoxLength"
         << ":  " << autopas::ArrayUtils::to_string(boxLength) << endl;
    cout << std::setw(valueOffset) << left << "Initial velocities"
         << ":  " << autopas::ArrayUtils::to_string(velocity) << endl;
    cout << std::setw(valueOffset) << left << "Particle Properties in Object:" << endl;
    cout << setw(valueOffset) << left << "Particle TypeId"
         << ":  " << typeId << endl;
    cout << setw(valueOffset) << left << "Particles Epsilon"
         << ":  " << epsilon << endl;
    cout << setw(valueOffset) << left << "Particles Sigma"
         << ":  " << sigma << endl;
    cout << setw(valueOffset) << left << "Particles Mass"
         << ":  " << mass << endl
         << endl;
  }

 private:
  static constexpr size_t valueOffset = 32;
  size_t numParticles;
  std::array<double, 3> boxLength;
  std::array<double, 3> bottomLeftCorner;
};
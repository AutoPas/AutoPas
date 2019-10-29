/**
 * @file CubeGrid.h
 * @author N. Fottner
 * @date 29/10/19
 */
#pragma once
class CubeGrid : public Object {
 public:
  /**
   * constructor for CubeGrid that is created in YamlParser and then included into the Simulation via the Generator
   * class
   * @param particlesPerDim
   * @param particleSpacing
   * @param velocity
   * @param bottomLeftCorner
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   */
  CubeGrid(const std::array<size_t, 3> &particlesPerDim, double particleSpacing,
           const std::array<double, 3> &bottomLeftCorner, const std::array<double, 3> &velocity_arg,
           const unsigned long &typeId_arg, const double &epsilon_arg, const double &sigma_arg, const double &mass_arg)
      : particlesPerDim(particlesPerDim),
        particleSpacing(particleSpacing),
        particlesTotal(particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]),
        bottomLeftCorner(bottomLeftCorner) {
    velocity = velocity_arg;
    typeId = typeId_arg;
    epsilon = epsilon_arg;
    sigma = sigma_arg;
    mass = mass_arg;
  }

  /**
   * Getter for ParticlesPerDim
   * @return particlePerDim
   */
  [[nodiscard]] const std::array<size_t, 3> &getParticlesPerDim() const { return particlesPerDim; }

      /**
       * Getter for ParticleSpacing
       * @return particleSpacing
       */
      [[nodiscard]] double getParticleSpacing() const {
    return particleSpacing;
  }

  /**
   * Getter for total number of Particles for object
   * @return particlesTotal
   */
  [[nodiscard]] size_t getParticlesTotal() const override { return particlesTotal; }

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
    std::array<double, 3> dppD;
    // copy for type conversion
    std::copy(std::begin(particlesPerDim), std::end(particlesPerDim), std::begin(dppD));
    return autopas::ArrayMath::add(bottomLeftCorner, (autopas::ArrayMath::mulScalar(dppD, particleSpacing)));
  }
  /**
   * Prints the Configuration of the current Object
   */
  void printConfig() override {
    using namespace std;

    cout << std::setw(valueOffset) << left << "Particles per dimension"
         << ":  " << autopas::ArrayUtils::to_string(particlesPerDim) << endl;
    cout << std::setw(valueOffset) << left << "Particle spacing"
         << ":  " << particleSpacing << endl;
    cout << std::setw(valueOffset) << left << "Number of Particles"
         << ":  " << (particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]) << endl;
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
  std::array<size_t, 3> particlesPerDim;
  double particleSpacing;
  size_t particlesTotal;
  std::array<double, 3> bottomLeftCorner;
};

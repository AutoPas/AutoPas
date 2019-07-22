#pragma once
#include <array>

class CubeGrid {
public:
    CubeGrid(const std::array<size_t, 3> &particlesPerDim, double particleSpacing,
             const std::array<double, 3> &velocity) : particlesPerDim(particlesPerDim),
                                                      particleSpacing(particleSpacing), velocity(velocity) {}

    CubeGrid() : particlesPerDim({20,20,20}),
                 particleSpacing(0.4), velocity({0.,0.,0.}) {}

    const array<size_t, 3> &getParticlesPerDim() const {
        return particlesPerDim;
    }

    double getParticleSpacing() const {
        return particleSpacing;
    }

    const array<double, 3> &getVelocity() const {
        return velocity;
    }

    std::string printConfig(){
        cout << "Grid generator" << endl;
      cout << setw(valueOffset) << left << "Particle spacing"
           << ":  " << particleSpacing << endl;

      cout << "Particles" << endl;
      cout << setw(valueOffset) << left << "  per dimension"
           << ":  " << particlesPerDim << endl;
      cout << setw(valueOffset) << left << "  total"
           << ":  " << (particlesPerDim * particlesPerDim * particlesPerDim) << endl;
    }


private:
    static constexpr size_t valueOffset = 32;
    std::array<size_t, 3> particlesPerDim;
    double particleSpacing;
    std::array<double, 3> velocity;

};

class CubeGauss{
public:
    CubeGauss(const std::array<double, 3> &boxLength, size_t numParticles, double distributionMean,
              double distributionStdDev, const std::array<double, 3> &velocity) : boxLength(boxLength),
                                                                                  numParticles(numParticles),
                                                                                  distributionMean(distributionMean),
                                                                                  distributionStdDev(
                                                                                          distributionStdDev),
                                                                                  velocity(velocity) {}

    const array<double, 3> &getBoxLength() const {
        return boxLength;
    }

    size_t getNumParticles() const {
        return numParticles;
    }

    double getDistributionMean() const {
        return distributionMean;
    }

    double getDistributionStdDev() const {
        return distributionStdDev;
    }

    const array<double, 3> &getVelocity() const {
        return velocity;
    }

private:
    static constexpr size_t valueOffset = 32;
    std::array<double, 3> boxLength;
    size_t numParticles;
    double distributionMean;
    double distributionStdDev;
    std::array<double, 3> velocity;


};

class CubeUniform  {

public:
    CubeUniform(const std::array<double, 3> &boxLength, size_t numParticles, const std::array<double, 3> &velocity)
            : boxLength(boxLength), numParticles(numParticles), velocity(velocity) {}

    const array<double, 3> &getBoxLength() const {
        return boxLength;
    }

    size_t getNumParticles() const {
        return numParticles;
    }

    const array<double, 3> &getVelocity() const {
        return velocity;
    }

private:
    static constexpr size_t valueOffset = 32;
    std::array<double, 3> boxLength;
    size_t numParticles;
    std::array<double, 3> velocity;

};
class Sphere {


public:
    Sphere(const std::array<double, 3> &center, int radius, double particleSpacing, unsigned long id,
           const std::array<double, 3> &velocity) : center(center), radius(radius), particleSpacing(particleSpacing),
                                                    id(id), velocity(velocity) {}

    const std::array<double, 3> &getCenter() const {
        return center;
    }

    int getRadius() const {
        return radius;
    }

    double getParticleSpacing() const {
        return particleSpacing;
    }

    unsigned long getId() const {
        return id;
    }

    const std::array<double, 3> &getVelocity() const {
        return velocity;
    }

private:
    static constexpr size_t valueOffset = 32;
    std::array<double, 3> center;
    int radius;
    double particleSpacing;
    unsigned long id;
    std::array<double, 3> velocity;

};
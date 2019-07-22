#pragma once
#include <array>
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ArrayMath.h"
#include "Generator.h"
class CubeGrid {
public:
    CubeGrid(const std::array<size_t, 3> &particlesPerDim, double particleSpacing,
             const std::array<double, 3> &velocity) : particlesPerDim(particlesPerDim),
                                                      particleSpacing(particleSpacing), velocity(velocity),particlesTotal(particlesPerDim[0]*particlesPerDim[1]*particlesPerDim[2]) {}

    const array<size_t, 3> &getParticlesPerDim() const {
        return particlesPerDim;
    }

    double getParticleSpacing() const {
        return particleSpacing;
    }

    const array<double, 3> &getVelocity() const {
        return velocity;
    }

    int getParticlesTotal() const {
        return particlesTotal;
    }

    void printConfig(){
        cout << setw(valueOffset) << left << "Particles per dimension"
             << ":  " << ArrayUtils::to_string(particlesPerDim) << endl;
      cout << setw(valueOffset) << left << "Particle spacing"
           << ":  " << particleSpacing << endl;
      cout << setw(valueOffset) << left << "Initial particle velocity"
           << ":  " << ArrayUtils::to_string(velocity) << endl;
        cout << setw(valueOffset) << left << "Number of Particles"
             << ":  " << (particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]) << endl;
    }


private:
    static constexpr size_t valueOffset = 32;
    std::array<size_t, 3> particlesPerDim;
    double particleSpacing;
    std::array<double, 3> velocity;
    int particlesTotal;

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
    void printConfig(){
        cout << setw(valueOffset) << left << "Box Length"
             << ":  " << ArrayUtils::to_string(boxLength) << endl;
        cout << setw(valueOffset) << left << "Distribution-Mean"
             << ":  " << distributionMean << endl;
        cout << setw(valueOffset) << left << "Distribution-StdDev"
             << ":  " << distributionStdDev << endl;
        cout << setw(valueOffset) << left << "NumberOfParticles"
             << ":  " <<numParticles << endl;
        cout << setw(valueOffset) << left << ""
             << ":  " << ArrayUtils::to_string(velocity) << endl;
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
    void printConfig(){
        cout << setw(valueOffset) << left << "Box Length"
             << ":  " << ArrayUtils::to_string(boxLength) << endl;
        cout << setw(valueOffset) << left << "NumberOfParticles"
             << ":  " <<numParticles << endl;
        cout << setw(valueOffset) << left << ""
             << ":  " << ArrayUtils::to_string(velocity) << endl;
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
    //@todo besser implementieren: (anderen Sphere Generator?)
    int particlesTotal(){
        int counter=0;
        for (int z = 0; z <= radius; ++z) {
            for (int y = 0; y <= radius; ++y) {
                for (int x = 0; x <= radius; ++x) {
                    std::array<double, 3> posDelta ={(double) x, (double) y,(double) z};
                    for (int i = -1; i <= 1; i += 2) {
                        for (int k = -1; k <= 1; k += 2) {
                            for (int l = -1; l <= 1; l += 2) {
                                std::array<double,3> multipliers = {(double) i, (double) k, (double) l};
                                std::array<double, 3> posVector =ArrayMath::add(center,ArrayMath::mulScalar(ArrayMath::mul(posDelta,multipliers),particleSpacing));
                                double disCheck=Generator::L2Norm(ArrayMath::sub(posVector,center));
                                if(disCheck<=(double)(radius+1)*particleSpacing) {
                                    counter ++;
                                }
                                if (z == 0)
                                    break;
                            }
                            if (y == 0)
                                break;
                        }
                        if (x == 0)
                            break;
                    }
                }
            }
        }
        return counter;
    }

    void printConfig(){
        cout << setw(valueOffset) << left << "Center of Sphere"
             << ":  " << ArrayUtils::to_string(center) << endl;
        cout << setw(valueOffset) << left << "radius"
             << ":  " <<radius << endl;
        cout << setw(valueOffset) << left << "particleSpacing"
             << ":  " << particleSpacing << endl;
//        cout << setw(valueOffset) << left << "first Particle in Sphere"
//             << ":  " << id << endl;
        cout << setw(valueOffset) << left << "NumberOfParticles"
             << ":  " << this->particlesTotal() << endl;
        cout << setw(valueOffset) << left << ""
             << ":  " << ArrayUtils::to_string(velocity) << endl;
    }



private:
    static constexpr size_t valueOffset = 32;
    std::array<double, 3> center;
    int radius;
    double particleSpacing;
    unsigned long id;
    std::array<double, 3> velocity;

};
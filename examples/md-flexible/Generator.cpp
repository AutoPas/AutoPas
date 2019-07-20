#include "Generator.h"



template <class Particle, class ParticleCell>
void Generator<Particle, ParticleCell>::initContainerGrid(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                                           size_t particlesPerDim, double particelSpacing) {
    double ppDxpS = (particlesPerDim)*particelSpacing;
    std::array<double, 3> boxMin({0., 0., 0.});

    std::array<double, 3> boxMax{ppDxpS, ppDxpS, ppDxpS};

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    Particle dummyParticle;
    GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                     {particelSpacing, particelSpacing, particelSpacing},
                                     {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

template <class Particle, class ParticleCell>
void Generator<Particle, ParticleCell>::initContainerGauss(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                                            double boxLength, size_t numParticles,
                                                            double distributionMean, double distributionStdDev) {
    std::array<double, 3> boxMin({0., 0., 0.});

    std::array<double, 3> boxMax{boxLength, boxLength, boxLength};

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    Particle dummyParticle;
    GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

template <class Particle, class ParticleCell>
void Generator<Particle, ParticleCell>::initContainerUniform(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                                              double boxLength, size_t numParticles) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    Particle dummyParticle;
    RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}

template <class Particle,class ParticleCell>
void Generator<Particle,ParticleCell>::Sphere(autopas::AutoPas<Particle, ParticleCell> &autopas,const std::array<double,3> &center,int radius, double particleSpacing,unsigned long id){
    for (int z = 0; z <= radius; ++z) { // generate circles along the z-axis; uses symmetry of sphere
        for (int y = 0; y <= radius; ++y) { // generate lines among the y-axis
            for (int x = 0; x <= radius; ++x) { // generate particles among the x-axis
                std::array<double, 3> posDelta ={(double) x, (double) y,
                                                 (double) z}; // offset of center as array
                double modifiedRadius = (double) radius + (double) radius / ((double) radius + 20.0); // use slightly increased radius for better looking circles
                if (L2Norm(posDelta) > modifiedRadius) // (x^2 + y^2 + z^2)^.5 = r
                    break;
                for (int i = -1; i <= 1; i += 2) { // mirror x-coordinate
                    for (int k = -1; k <= 1; k += 2) { // mirror y-coordinate
                        for (int l = -1; l <= 1; l += 2) { // mirror z-coordinate
                            double multipliers[3] = {(double) i, (double) k, (double) l}; // multipliers for mirroring

                            std::array<double, 3> posVector = center + (arrayMul(arrayMul(posDelta,multipliers),particleSpacing)); // actual coordinates of new particle

                            Particle p(posVector, {0.,0.,0.},id);
                            autopas.addParticle(p);

                            if (z == 0) // prevent duplicates
                                break;
                        }
                        if (y == 0) // prevent duplicates
                            break;
                    }
                    if (x == 0) // prevent duplicates
                        break;
                }
            }
        }
    }
}

template <class Particle,class ParticleCell>
void Generator<Particle,ParticleCell>::SphereV(autopas::AutoPas<Particle, ParticleCell> &autopas,const std::array<double,3> &velocity,const std::array<double,3> &center,int radius, double particleSpacing,unsigned long id){
    for (int z = 0; z <= radius; ++z) { // generate circles along the z-axis; uses symmetry of sphere
        for (int y = 0; y <= radius; ++y) { // generate lines among the y-axis
            for (int x = 0; x <= radius; ++x) { // generate particles among the x-axis
                std::array<double, 3> posDelta ={(double) x, (double) y,
                                                 (double) z}; // offset of center as array
                double modifiedRadius = (double) radius + (double) radius / ((double) radius + 20.0); // use slightly increased radius for better looking circles
                if (L2Norm(posDelta) > modifiedRadius) // (x^2 + y^2 + z^2)^.5 = r
                    break;
                for (int i = -1; i <= 1; i += 2) { // mirror x-coordinate
                    for (int k = -1; k <= 1; k += 2) { // mirror y-coordinate
                        for (int l = -1; l <= 1; l += 2) { // mirror z-coordinate
                            double multipliers[3] = {(double) i, (double) k, (double) l}; // multipliers for mirroring

                            std::array<double, 3> posVector = center + (arrayMul(arrayMul(posDelta,multipliers),particleSpacing)); // actual coordinates of new particle

                            Particle p(posVector, velocity,id);
                            autopas.addParticle(p);

                            if (z == 0) // prevent duplicates
                                break;
                        }
                        if (y == 0) // prevent duplicates
                            break;
                    }
                    if (x == 0) // prevent duplicates
                        break;
                }
            }
        }
    }
}

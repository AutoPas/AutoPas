#pragma once

#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "PrintableMolecule.h"
#include "autopas/AutoPas.h"
#include <vector>

/**Class for contructing a container and generating Objects and Shapes filled with Particles
 * */
class Generator {
public:

    static double L2Norm(std::array<double, 3> array) {
        double square_sum = 0;
        for (unsigned int i = 0; i < array.size(); i++) {
            square_sum += (array[i] * array[i]);
        }
        return sqrt(square_sum);
    }

    /**
    * @brief Constructs a container and fills it with particles.
    *
 * According to the options passed, a %DirectSum or %'LinkedCells' container is
 * built. It consists of %`FullContainers` and is filled with
 * `PrintableMolecules`. The particles are aligned on a cuboid grid.
 *
 * @param autopas AutoPas object that should be initialized
 * @param particlesPerDim Number of desired particles per dimension.
 * @param particelSpacing Space between two particles along each axis of space.
 */
    template <class Particle,class ParticleCell>
    static void CubeGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, std::array<size_t,3> particlesPerDim,
                  double particelSpacing);
    template <class Particle,class ParticleCell>
    static void Gauss(autopas::AutoPas<Particle, ParticleCell> &autopas, std::array<double, 3> boxLength, size_t numParticles,
               double distributionMean, double distributionStdDev);
    template <class Particle,class ParticleCell>
    static void Random(autopas::AutoPas<Particle, ParticleCell> &autopas, double boxLength, size_t numParticles);

    /**Generates a Sphere with @param radius number of Particles with initial velocity 0
     * @param autopas
     * @param center
     * @param radius
     * @param particleSpacing
     * @param id
     * */
    template <class Particle,class ParticleCell>
    static void Sphere(autopas::AutoPas<Particle, ParticleCell> &autopas,const std::array<double, 3> &center,int radius, double particleSpacing,unsigned long id);

    /**Generates a Sphere with @param radius number of Particles with initial @param velocity
    * @param Autopas
    * @param center
    * @param radius
    * @param velocity
    * @param particleSpacing
    * @param id
    * */
    template <class Particle,class ParticleCell>
    static void SphereV(autopas::AutoPas<Particle, ParticleCell> &autopas,const std::array<double, 3> &center,const std::array<double,3> &velocity,int radius, double particleSpacing,unsigned long id);
};

template <class Particle, class ParticleCell>
void Generator::CubeGrid(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                                 std::array<size_t,3> particlesPerDim, double particelSpacing) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax{(particlesPerDim[0])*particelSpacing, (particlesPerDim[1])*particelSpacing, (particlesPerDim[2])*particelSpacing};

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    Particle dummyParticle;
    GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                     {particelSpacing, particelSpacing, particelSpacing},
                                     {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

template <class Particle, class ParticleCell>
void Generator::Gauss(autopas::AutoPas<Particle, ParticleCell> &autopas,
                      std::array<double, 3> boxLength, size_t numParticles,
                                              double distributionMean, double distributionStdDev) {
    std::array<double, 3> boxMin({0., 0., 0.});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxLength);

    autopas.init();

    Particle dummyParticle;
    GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

template <class Particle, class ParticleCell>
void Generator::Random(autopas::AutoPas<Particle, ParticleCell> &autopas,
                                               double boxLength, size_t numParticles) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    Particle dummyParticle;
    RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}
//@todo add Type ID
template <class Particle,class ParticleCell>
void Generator::Sphere(autopas::AutoPas<Particle, ParticleCell> &autopas,const std::array<double,3> &center,int radius, double particleSpacing,unsigned long id){
    //@todo autopas sets und init eher in simulation.initialization machen sonst muss man boxlenght anpassen mit particlespacing und radius

    //@todo function schreiben die die passende initialisierung von boxMin und boxMax macht. wenn man zb 2 object initialisieren will
    auto boxLength = (double)(radius+1)*particleSpacing;
            //2* (double)radius * particleSpacing +1. ;
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});
    std::array<double, 3> boxMin({-1.0*boxLength, -1.0*boxLength, -1.*boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    for (int z = 0; z <= radius; ++z) { // generate circles along the z-axis; uses symmetry of sphere
        for (int y = 0; y <= radius; ++y) { // generate lines among the y-axis
            for (int x = 0; x <= radius; ++x) { // generate particles among the x-axis
                std::array<double, 3> posDelta ={(double) x, (double) y,
                                                 (double) z}; // offset of center as array
                for (int i = -1; i <= 1; i += 2) { // mirror x-coordinate
                    for (int k = -1; k <= 1; k += 2) { // mirror y-coordinate
                        for (int l = -1; l <= 1; l += 2) { // mirror z-coordinate
                            std::array<double,3> multipliers = {(double) i, (double) k, (double) l}; // multipliers for mirroring
                            std::array<double, 3> posVector =ArrayMath::add(center,ArrayMath::mulScalar(ArrayMath::mul(posDelta,multipliers),particleSpacing)); // actual coordinates of new particle
                            double disCheck=L2Norm(ArrayMath::sub(posVector,center));
                            if(disCheck<=(double)(radius+1)*particleSpacing) {
                                Particle p(posVector, {0., 0., 0.}, id);
                                autopas.addParticle(p);
                                id ++;
                                }
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
void Generator::SphereV(autopas::AutoPas<Particle, ParticleCell> &autopas,const std::array<double,3> &velocity,const std::array<double,3> &center,int radius, double particleSpacing,unsigned long id){
    auto boxLength = (double)(radius+1)*particleSpacing;
    //2* (double)radius * particleSpacing +1. ;
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});
    std::array<double, 3> boxMin({-1.0*boxLength, -1.0*boxLength, -1.*boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    for (int z = 0; z <= radius; ++z) { // generate circles along the z-axis; uses symmetry of sphere
        for (int y = 0; y <= radius; ++y) { // generate lines among the y-axis
            for (int x = 0; x <= radius; ++x) { // generate particles among the x-axis
                std::array<double, 3> posDelta ={(double) x, (double) y,
                                                 (double) z}; // offset of center as array
                for (int i = -1; i <= 1; i += 2) { // mirror x-coordinate
                    for (int k = -1; k <= 1; k += 2) { // mirror y-coordinate
                        for (int l = -1; l <= 1; l += 2) { // mirror z-coordinate
                            std::array<double,3> multipliers = {(double) i, (double) k, (double) l}; // multipliers for mirroring
                            std::array<double, 3> posVector =ArrayMath::add(center,ArrayMath::mulScalar(ArrayMath::mul(posDelta,multipliers),particleSpacing)); // actual coordinates of new particle
                            double disCheck=L2Norm(ArrayMath::sub(posVector,center));
                            if(disCheck<=(double)(radius+1)*particleSpacing) {
                                Particle p(posVector, velocity, id);
                                autopas.addParticle(p);
                                id ++;
                            }
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

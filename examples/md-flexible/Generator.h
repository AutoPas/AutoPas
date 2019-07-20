#pragma once

#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "PrintableMolecule.h"
#include "autopas/AutoPas.h"
#include <vector>

/**Class for contructing a container and generating Objects and Shapes filled with Particles
 * */
template <class Particle,class ParticleCell>
class Generator {

    Generator()=default;
    ~Generator()=default;


    double L2Norm(std::array<double, 3> array) {
        double square_sum = 0;
        for (unsigned int i = 0; i < array.size(); i++) {
            square_sum += (array[i] * array[i]);
        }
        return sqrt(square_sum);
    }

    std::array<double,3> arrayMul(std::array<double,3> &array,double scalar){
        for(auto e: array){
            e*=scalar;
        }
        return array;
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

    void initContainerGrid(autopas::AutoPas<Particle, ParticleCell> &autopas, size_t particlesPerDim,
                           double particelSpacing);

    void initContainerGauss(autopas::AutoPas<Particle, ParticleCell> &autopas, double boxLength, size_t numParticles,
                            double distributionMean, double distributionStdDev);

    void initContainerUniform(autopas::AutoPas<Particle, ParticleCell> &autopas, double boxLength, size_t numParticles);

    /**Generates a Sphere with @param radius number of Particles with initial velocity 0
     * @param autopas
     * @param center
     * @param radius
     * @param particleSpacing
     * @param id
     * */
    void Sphere(autopas::AutoPas<Particle, ParticleCell> &autopas,const std::array<double, 3> &center,int radius, double particleSpacing,unsigned long id);

    /**Generates a Sphere with @param radius number of Particles with initial @param velocity
    * @param Autopas
    * @param center
    * @param radius
    * @param velocity
    * @param particleSpacing
    * @param id
    * */
    void SphereV(autopas::AutoPas<Particle, ParticleCell> &autopas,const std::array<double, 3> &center,const std::array<double,3> &velocity,int radius, double particleSpacing,unsigned long id);
};

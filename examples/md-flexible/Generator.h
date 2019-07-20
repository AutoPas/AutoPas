#pragma once

#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "PrintableMolecule.h"
#include "autopas/AutoPas.h"

/**Class for contructing a container and generating Objects and Shapes filled with Particles
 * */
template <class Particle,class ParticleCell>
class Generator {

    /**
 * @brief Constructs a container and fills it with particles.
 *
 * According to the options passed, a %DirectSum or %'LinkedCells' container is
 * built. It consists of %`FullParticleCells` and is filled with
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




};

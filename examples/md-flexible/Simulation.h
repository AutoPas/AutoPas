//
// Created by nicola on 12.05.19.
//

#ifndef AUTOPAS_SIMULATION_H
#define AUTOPAS_SIMULATION_H

#include <autopas/utils/MemoryProfiler.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "MDFlexParser.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"


using namespace autopas;
using namespace std;

template<class Particle,class ParticleCell>
class  Simulation {

private:
    AutoPas<Particle,ParticleCell>* _autopas;



public:
    Simulation();

    virtual ~Simulation();  //@todo

    Simulation(const AutoPas<Particle, ParticleCell> &autopas);

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
    void initContainerGrid(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                           size_t particlesPerDim, double particelSpacing);

    void initContainerGauss(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                            double boxLength, size_t numParticles, double distributionMean, double distributionStdDev);

    void initContainerUniform(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                              double boxLength, size_t numParticles);

    /** @brief This function is needed to create functors with the actual type through templates.
     * This function is used for the Force calculation with an desired Functor and "iterate Pairwise"
     * @tparam FunctorChoice
     * @tparam AutoPasTemplate
     * @param autopas
     * @param cutoff
     * @param numIterations
     * @return Time for all calculation iterations in microseconds.
     */
    template <class FunctorChoice, class AutoPasTemplate>
    long calculate(AutoPasTemplate &autopas, double cutoff, size_t numIterations);

    /** @brief This function processes the main simulation loop
     * -calls the time discretization class(calculate fores, etc ...)
     * -do the output each timestep
     */
    void simulate(){




        setautopas(time.intergration(autopas));

        classe: timeintegration -


        #pragma omp parallel
        for(auto iter = _autopas->getContainer().begin(); iter.isValid(); ++iter) {
            // user code:
            auto position = iter->getR();
        }
    }

    void calculateV();

    void calculateX();






    /** @brief This function
     * -sets/initializes the simulation domain with the particles generators
     * @todo -initialized Velocities and Positions (and forces?)
     */
    void initialize(MDFlexParser parser){
        auto generatorChoice=parser.getGeneratorOption();
        auto particlesPerDim=parser.getParticlesPerDim();
        auto particleSpacing=parser.getParticleSpacing();
        auto particlesTotal=parser.getParticlesTotal();
        auto distributionMean=parser.getDistributionMean();
        auto distributionStdDev=parser.getDistributionStdDev();
        auto boxLength=parser.getBoxLength();

        switch (generatorChoice) {
            case MDFlexParser::GeneratorOption::grid: {
                initContainerGrid(_autopas, particlesPerDim, particleSpacing);
                particlesTotal = particlesPerDim * particlesPerDim * particlesPerDim;
                break;
            }
            case MDFlexParser::GeneratorOption::uniform: {
                initContainerUniform(_autopas, boxLength, particlesTotal);
                break;
            }
            case MDFlexParser::GeneratorOption::gaussian: {
                initContainerGauss(_autopas, boxLength, particlesTotal, distributionMean, distributionStdDev);
                break;
            }
            default:
                std::string errorMessage= std::string("Unknown generetor Choice");
                throw std::runtime_error(errorMessage);
        }

        //  @todo initialize Velocities and Positions
    }

    /**Getter for Autopas
     * @return autopasObject
     */
    const AutoPas<Particle, ParticleCell> &getAutopas() const;


};

template<class Particle, class ParticleCell>
const AutoPas<Particle, ParticleCell> &Simulation<Particle, ParticleCell>::getAutopas() const {
    return _autopas;
}

template<class Particle, class ParticleCell>
Simulation<Particle, ParticleCell>::Simulation(const AutoPas<Particle, ParticleCell> &autopas):_autopas(autopas) {}

template<class Particle,class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerGrid(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                       size_t particlesPerDim, double particelSpacing) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax(
            {(particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                     {particelSpacing, particelSpacing, particelSpacing},
                                     {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

template<class Particle,class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerGauss(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                        double boxLength, size_t numParticles, double distributionMean, double distributionStdDev) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

template<class Particle,class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerUniform(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                          double boxLength, size_t numParticles) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}


template <class Particle,class ParticleCell>
template <class FunctorChoice, class AutoPasTemplate>
long Simulation<Particle, ParticleCell>::calculate(AutoPasTemplate &autopas, double cutoff, size_t numIterations) {

    auto functor = FunctorChoice(cutoff, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);

    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;

    startCalc = std::chrono::high_resolution_clock::now();

    // actual Calculation
    for (unsigned int i = 0; i < numIterations; ++i) {
        if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
            cout << "Iteration " << i << endl;
            cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << endl;
        }
        autopas.iteratePairwise(&functor);
    }
    stopCalc = std::chrono::high_resolution_clock::now();

    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;
}

#endif //AUTOPAS_SIMULATION_H

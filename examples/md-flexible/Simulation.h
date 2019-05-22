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
#include "TimeDiscretization.h"


using namespace autopas;
using namespace std;

template<class Particle,class ParticleCell>
class  Simulation {

private:
    AutoPas<Particle, ParticleCell> *_autopas;

    //@todo nicht sicher mit der implementierung-> anders: functor initialize -> simulate(S.init)
    autopas::Functor Functor;

    long durationX;
    long durationF;
    long durationV;
public:
    virtual ~Simulation();  //@todo

    explicit Simulation(const AutoPas<Particle, ParticleCell> &autopas);

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
    /**Getter for Duration of Position Calculation
     * @return durationX
     */
    long getDurationX() const;
    /**Adder for durationX
     * @param durationX
     */
    void addDurationX(long durationX);
    /**Getter for Duration of Force Calculation
     * @return durationF
     */
    long getDurationF() const;
    /**Adder for durationF
     * @param durationF
     */
    void addDurationF(long durationF);
    /**Getter for Duration of Velocity Calculation
     * @return durationV
     */
    long getDurationV() const;
    /**Adder for durationV
     * @param durationV
     */
    void addDurationV(long durationV);

    void initContainerGauss(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                            double boxLength, size_t numParticles, double distributionMean, double distributionStdDev);

    void initContainerUniform(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                              double boxLength, size_t numParticles);

    /** @brief This function
    * -initializes the autopas Object
    * -sets/initializes the simulation domain with the particles generators
    * @todo -initialized Velocities and Positions (and forces?)
    */
    void initialize(MDFlexParser parser);

    /**
     * Does the ForceCalculation
     * @param Force Calculation Functor
     * @return Duration of Calculation
     * */
    void CalcF(autopas::Functor functor);

    /**
     * This function processes the main simulation loop
     * -calls the time discretization class(calculate fores, etc ...)
     * -do the output each timestep
     * -collects the duration of every Calculation(Position,Force,Velocity)
     * @param duration of Simulation
     */
    long simulate();

    /**Getter for Autopas Oject
     * @return Autopas Object
     */
    AutoPas<Particle, ParticleCell> *getAutopas() const;

    void setFunctor(auto functor);

};


template<class Particle, class ParticleCell>
AutoPas<Particle, ParticleCell> *Simulation<Particle, ParticleCell>::getAutopas() const {
    return _autopas;
}

template<class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(MDFlexParser parser){

    string logFileName(parser.getLogFileName());

    auto measureFlops(parser.getMeasureFlops());
    auto numIterations(parser.getIterations());
    auto particlesTotal(parser.getParticlesTotal());
    auto verletRebuildFrequency(parser.getVerletRebuildFrequency());
    auto vtkFilename(parser.getWriteVTK()); //
    auto boxLength(parser.getBoxLength());
    auto containerChoice(parser.getContainerOptions());
    auto selectorStrategy(parser.getSelectorStrategy());
    auto cutoff(parser.getCutoff());
    auto dataLayoutOptions(parser.getDataLayoutOptions());
    auto distributionMean(parser.getDistributionMean());
    auto distributionStdDev(parser.getDistributionStdDev());
    auto functorChoice(parser.getFunctorOption());
    auto generatorChoice(parser.getGeneratorOption());
    auto logLevel(parser.getLogLevel());
    auto newton3Options(parser.getNewton3Options());
    auto particleSpacing(parser.getParticleSpacing());
    auto particlesPerDim(parser.getParticlesPerDim());
    auto traversalOptions(parser.getTraversalOptions());
    auto tuningInterval(parser.getTuningInterval());
    auto tuningSamples(parser.getTuningSamples());
    auto verletSkinRadius(parser.getVerletSkinRadius());
    _autopas.setCutoff(cutoff);
    _autopas.setVerletSkin(verletSkinRadius);
    _autopas.setVerletRebuildFrequency(verletRebuildFrequency);
    _autopas.setTuningInterval(tuningInterval);
    _autopas.setNumSamples(tuningSamples);
    _autopas.setSelectorStrategy(selectorStrategy);
    _autopas.setAllowedContainers(containerChoice);
    _autopas.setAllowedTraversals(traversalOptions);
    _autopas.setAllowedDataLayouts(dataLayoutOptions);
    _autopas.setAllowedNewton3Options(newton3Options);

    switch (functorChoice) {
        case MDFlexParser::FunctorOption::lj12_6: {
            this->setFunctor(LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>);
        }
        case MDFlexParser::FunctorOption::lj12_6_AVX: {
            this->setFunctor(LJFunctorAVX<PrintableMolecule, FullParticleCell<PrintableMolecule>>);
        }
    }


        //  @todo initialize Velocities and Positions---JA oder NEIN?
    TimeDiscretization td;
    td.VSCalculateV(this->_autopas);
    td.VSCalculateX(this->_autopas);

    }

template<class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::setFunctor(auto functor) {
    Functor = functor;
}

template<class Particle, class ParticleCell>
Simulation<Particle, ParticleCell>::Simulation(const AutoPas<Particle, ParticleCell> &autopas):_autopas(autopas) {
    durationF=0; durationV=0; durationX=0;
}

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


template<class Particle,ParticleCell>
void Simulation<Particle,ParticleCell>::CalcF(autopas::Functor functor){
    td::chrono::high_resolution_clock::time_point startCalc, stopCalc;
    startCalc = std::chrono::high_resolution_clock::now();
    _autopas->iteratePairwise(functor)
    stopCalc = std::chrono::high_resolution_clock::now();
    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    this->addDurationF(durationCalc);
}


template<class Particle,class ParticleCell>
long Simulation<Particle, ParticleCell>::simulate(){



    auto functor = LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>(cutoff, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);
    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
    startCalc = std::chrono::high_resolution_clock::now();


    double particleDelta_T=0.01;    //@todo -get Value from Parser
    double time=0;
    double timeEnd=10;              //@todo -get TimeEnd from Parser
    TimeDiscretization td(particleDelta_T);
    //main simulation loop
    while(time<timeEnd){
        this->addDurationX(td.VSCalculateX(_autopas));
        // -> nicht sicher ob man das IF-Case braucht
        if (this->autopas::Logger::get()->level() <= this->autopas::Logger::LogLevel::debug) {
             cout << "Iteration " << time/particleDelta_T << endl;
            cout << "Current Memory usage: " << this->autopas::memoryProfiler::currentMemoryUsage() << " kB" << endl;
        }
        this->CalcF();
        this->addDurationV(td.VSCalculateV(_autopas));
        time+=particleDelta_T;
    }



    stopCalc = std::chrono::high_resolution_clock::now();
    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;



}

template<class Particle, class ParticleCell>
long Simulation<Particle, ParticleCell>::getDurationX() const {
    return durationX;
}

template<class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::addDurationX(long durationX) {
    Simulation::durationX += durationX;
}

template<class Particle, class ParticleCell>
long Simulation<Particle, ParticleCell>::getDurationF() const {
    return durationF;
}

template<class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::addDurationF(long durationF) {
    Simulation::durationF += durationF;
}

template<class Particle, class ParticleCell>
long Simulation<Particle, ParticleCell>::getDurationV() const {
    return durationV;
}

template<class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::addDurationV(long durationV) {
    Simulation::durationV += durationV;
}


#endif //AUTOPAS_SIMULATION_H

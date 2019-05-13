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
#include "MDFlexParser.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"


using namespace autopas;
using namespace std;

template<class Particle,class ParticleCell>
class  Simulation {

private:
    AutoPas<Particle,ParticleCell> _autopas;

public:
    Simulation();

    Simulation(const AutoPas<Particle, ParticleCell> &autopas);

    /** @brief This function is needed to create functors with the actual type through templates.
     * @tparam FunctorChoice
     * @tparam AutoPasTemplate
     * @param autopas
     * @param cutoff
     * @param numIterations
     * @return Time for all calculation iterations in microseconds.
     */
    template <class FunctorChoice, class AutoPasTemplate>
    long calculate(AutoPasTemplate &autopas, double cutoff, size_t numIterations){
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



    /** @brief This function processes the main simulation loop
     * -calls the time discretization class(calculate fores, etc ...)
     * -do the output each timestep
     */
    void simulate();
    /** @brief This function
     * -sets/initializes the simulation domain with the particles generators
     */
    void initialize(){

    }

    virtual ~Simulation();
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



#endif //AUTOPAS_SIMULATION_H

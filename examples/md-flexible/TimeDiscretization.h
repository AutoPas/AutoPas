//
// Created by nicola on 13.05.19.
//

#ifndef AUTOPAS_TIMEDISCRETIZATION_H
#define AUTOPAS_TIMEDISCRETIZATION_H

#include <autopas/utils/MemoryProfiler.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "autopas/AutoPas.h"
#include "Simulation.h"

//idee: Die TimeDiskretization Methode als Template weiterzugeben, und somit pro verschiede Diskretisierungsmethode eine Klasse/Functor zu schreiben


template <class Particle,class ParticleCell, class FunctorChoice>
class TimeDiscretization {

public:

    virtual ~TimeDiscretization() {
    //@todo
    }

    TimeDiscretization(FunctorChoice functor, AutoPas<Particle, ParticleCell> *autopas) : Functor(functor),
                                                                                          autopas(autopas) {}

    void VerletStörmerTime();


    void VSCalcuteX();

    void VSCalculateV();


    //@todo  rückgabewert__long???
    long FunctorCalculate(double cutoff, size_t numIterations);


private:

    /** Force Calculation Functors
     */
    FunctorChoice  Functor;

    AutoPas<Particle,ParticleCell>* autopas;


};

template <class Particle,class ParticleCell,class FunctorChoice>
long TimeDiscretization<Particle, ParticleCell,FunctorChoice>::FunctorCalculate(double cutoff, size_t numIterations) {

    auto functor = Functor(cutoff, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);

    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;

    startCalc = std::chrono::high_resolution_clock::now();

    // actual Calculation
    for (unsigned int i = 0; i < numIterations; ++i) {
        if (this->autopas::Logger::get()->level() <= this->autopas::Logger::LogLevel::debug) {
            cout << "Iteration " << i << endl;
            cout << "Current Memory usage: " << this->autopas::memoryProfiler::currentMemoryUsage() << " kB" << endl;
        }
        this->autopas.iteratePairwise(&functor);
    }
    stopCalc = std::chrono::high_resolution_clock::now();

    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;
}


#endif //AUTOPAS_TIMEDISCRETIZATION_H

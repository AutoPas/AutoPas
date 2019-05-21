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
#include "autopas/utils/ArrayMath.h"


//idee: Die TimeDiskretization Methode als Template weiterzugeben, und somit pro verschiede Diskretisierungsmethode eine Klasse/Functor zu schreiben


template <class Particle>
class TimeDiscretization {

public:
    TimeDiscretization(double particleDeltaT);

    virtual ~TimeDiscretization() {
    //@todo
    }

    /**Calculate the new Position for every Praticle using the Iterator and the Störmer-Verlet Algorithm
     */
    long VSCalculateX(autopas::AutoPas<Particle, autopas::ParticleCell<Particle>> &autopas);

    /**Calculate the new Velocity for every Praticle using the Iterator and the Störmer-Verlet Algorithm
     */
    long VSCalculateV(autopas::AutoPas<Particle, autopas::ParticleCell<Particle>> &autopas);

    /**Proceeds the Calculation of Forces with the given functor and "autopas.iteratePairwise"
     * @param cutoff
     */

    double getParticleDeltaT() const;


private:
    /** \brief Duration of a timestep
     * */
    double particle_delta_t;        // @todo: add particle_delta, time_end -> main and parser

};

template <class Particle>
long TimeDiscretization<Particle>::VSCalculateX(autopas::AutoPas<Particle, autopas::ParticleCell<Particle>> &autopas) {
    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
    startCalc = std::chrono::high_resolution_clock::now();
#pragma omp parallel
    for (auto iter = autopas->getContainer().begin(); iter.isValid(); ++iter) {
        auto r = iter->getR();
        auto v = iter->getV();
        auto m = autopas::MoleculeLJ::getMass();
        auto f = iter->getF();
        v = autopas::ArrayMath::mulScalar(v, this->particle_delta_t);
        f= autopas::ArrayMath::mulScalar(f,(particle_delta_t * particle_delta_t / (2 * m)));
        auto newR = autopas::ArrayMath::add(v,f);
        iter->addR(newR);
    }
    stopCalc = std::chrono::high_resolution_clock::now();
    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;
}


template <class Particle>
long TimeDiscretization<Particle>::VSCalculateV(autopas::AutoPas<Particle, autopas::ParticleCell<Particle>> &autopas) {
    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
    startCalc = std::chrono::high_resolution_clock::now();
#pragma omp parallel
    for (auto iter = autopas->getContainer().begin(); iter.isValid(); ++iter) {
        auto v =iter->getV();
        auto m = autopas::MoleculeLJ::getMass();
        auto force = iter->getF();
        auto old_force= autopas::MoleculeLJ::getOldf();  //@todo : implement right old_force interface (Particle -> inherent to MoleculeLJ -> PrintableMolecule)
        auto newV= ArrayMath((autopas::ArrayMath::add(force,old_force)),particle_delta_t/(2*m));
        iter->addV(newV);
    }
    stopCalc = std::chrono::high_resolution_clock::now();
    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;
}

template<class Particle>
TimeDiscretization<Particle>::TimeDiscretization(double particleDeltaT):particle_delta_t(particleDeltaT) {}

template<class Particle>
double TimeDiscretization<Particle>::getParticleDeltaT() const {
    return particle_delta_t;
}

#endif //AUTOPAS_TIMEDISCRETIZATION_H

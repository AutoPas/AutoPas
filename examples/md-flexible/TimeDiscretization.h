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


template <class Particle,class ParticleCell, class FunctorChoice>
class TimeDiscretization {

public:

    virtual ~TimeDiscretization() {
    //@todo
    }

    TimeDiscretization(FunctorChoice functor, AutoPas<Particle, ParticleCell> *autopas,double particle_delta_t) : Functor(functor),
                                                                                          autopas(autopas),particle_delta_t(particle_delta_t) {}

    /**Applies one timestep in the Simulation with the Stoermer-Verlet Algorithm
     * @tparam Particle
     * @tparam ParticleCell
     * @tparam FunctorChoice
     */
    long VerletStörmerTime();

    /**Calculate the new Position for every Praticle using the Iterator and the Störmer-Verlet Algorithm
     */
    void VSCalculateX();

    /**Calculate the new Velocity for every Praticle using the Iterator and the Störmer-Verlet Algorithm
     */
    void VSCalculateV();

    /**Proceeds the Calculation of Forces with the given functor and "autopas.iteratePairwise"
     * @param cutoff
     */
    void VSCalculateF(double cutoff);


private:

    /** \brief Force Calculation Functors
     */
    FunctorChoice  Functor;

    /** \brief Autopas Object
     */
    AutoPas<Particle,ParticleCell>* autopas;

    /** \brief Duration of a timestep
     * */
    double particle_delta_t;        // @todo: add particle_delta, time_end -> main and parser

};

template <class Particle,class ParticleCell,class FunctorChoice>
long TimeDiscretization<Particle, ParticleCell,FunctorChoice>::VerletStörmerTime() {
    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
    startCalc = std::chrono::high_resolution_clock::now();
    this->VSCalculateX();
    this->VSCalculateF();
    this->VSCalculateV();
    stopCalc = std::chrono::high_resolution_clock::now();
    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;
}


template <class Particle,class ParticleCell,class FunctorChoice>
void TimeDiscretization<Particle, ParticleCell,FunctorChoice>::VSCalculateF(double cutoff) {
    auto functor = Functor(cutoff, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);
    if (this->autopas::Logger::get()->level() <= this->autopas::Logger::LogLevel::debug) {
        // cout << "Iteration " << i << endl;  --> man kann einfach mit delta_t und end_t die current iteration ausgeben
        cout << "Current Memory usage: " << this->autopas::memoryProfiler::currentMemoryUsage() << " kB" << endl;
    }
    this->autopas.iteratePairwise(&functor);
}

template <class Particle,class ParticleCell,class FunctorChoice>
void TimeDiscretization<Particle, ParticleCell,FunctorChoice>::VSCalculateX() {
#pragma omp parallel
    for (auto iter = this->autopas->getContainer().begin(); iter.isValid(); ++iter) {
        auto r = iter->getR();
        auto v = iter->getV();
        auto m = MoleculeLJ::getMass();
        auto f = iter->getF();
        v = ArrayMath::mulScalar(v, this->particle_delta_t);
        f= ArrayMath::mulScalar(f,(particle_delta_t * particle_delta_t / (2 * m)));
        auto newR = ArrayMath::add(v,f);
        iter->addR(newR);
    }
}


template <class Particle,class ParticleCell,class FunctorChoice>
void TimeDiscretization<Particle, ParticleCell,FunctorChoice>::VSCalculateV() {
#pragma omp parallel
    for (auto iter = this->autopas->getContainer().begin(); iter.isValid(); ++iter) {
        auto v =iter->getV();
        auto m = MoleculeLJ::getMass();
        auto force = iter->getF();
        auto old_force= MoleculeLJ::getOldf();  //@todo : implement right old_force interface (Particle -> inherent to MoleculeLJ -> PrintableMolecule)
        auto newV= ArrayMath((ArrayMath::add(force,old_force)),particle_delta_t/(2*m));
        iter->addV(newV);
    }
}

#endif //AUTOPAS_TIMEDISCRETIZATION_H

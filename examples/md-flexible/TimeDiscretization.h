//
// Created by nicola on 13.05.19.
//

#pragma once
#include <autopas/utils/MemoryProfiler.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"


using namespace std;

template <class AutoPasTemplate>
class TimeDiscretization {

public:
    TimeDiscretization(double particleDeltaT);

    virtual ~TimeDiscretization() {
    }
    string arrayString(array<double,3> input){
        std::ostringstream os;
        for (double i: input){
            os << i;
            os << " _ ";
        }
        std::string str(os.str());
        return str;
    }
    /**Calculate the new Position for every Praticle using the Iterator and the Störmer-Verlet Algorithm
     */
    long VSCalculateX(AutoPasTemplate autopas);

    /**Calculate the new Velocity for every Praticle using the Iterator and the Störmer-Verlet Algorithm
     */
    long VSCalculateV(AutoPasTemplate autopas);

    double getParticleDeltaT() const;


private:
    /**  Duration of a timestep
     * */
    double particle_delta_t;

};

template <class AutoPasTemplate>
long TimeDiscretization<AutoPasTemplate>::VSCalculateX(AutoPasTemplate autopas) {
    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
    startCalc = std::chrono::high_resolution_clock::now();
#pragma omp parallel
    for (auto iter = autopas->getContainer()->begin(); iter.isValid(); ++iter) {
        auto v = iter->getV();
        auto m = iter->getMass();
        auto f = iter->getF();
        iter->setOldf(f);
        v = autopas::ArrayMath::mulScalar(v, particle_delta_t);
        f= autopas::ArrayMath::mulScalar(f,(particle_delta_t * particle_delta_t / (2 * m)));
        auto newR = autopas::ArrayMath::add(v,f);
        iter->addR(newR);
    }
    stopCalc = std::chrono::high_resolution_clock::now();
    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;
}


template <class AutoPasTemplate>
long TimeDiscretization<AutoPasTemplate>::VSCalculateV(AutoPasTemplate autopas) {
    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
    startCalc = std::chrono::high_resolution_clock::now();
#pragma omp parallel
    for (auto iter = autopas->getContainer()->begin(); iter.isValid(); ++iter) {
        auto m = iter->getMass();
        auto force = iter->getF();
        auto old_force= iter->getOldf();
        auto newV = autopas::ArrayMath::mulScalar((autopas::ArrayMath::add(force,old_force)) , particle_delta_t/(2*m));
        iter->addV(newV);
    }
    stopCalc = std::chrono::high_resolution_clock::now();
    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;
}

template<class AutoPasTemplate>
TimeDiscretization<AutoPasTemplate>::TimeDiscretization(double particleDeltaT):particle_delta_t(particleDeltaT) {}

template<class AutoPasTemplate>
double TimeDiscretization<AutoPasTemplate>::getParticleDeltaT() const {
    return particle_delta_t;
}
//
// Created by seckler on 19.01.18.
//

#ifndef AUTOPAS_SPHPARTICLE_H
#define AUTOPAS_SPHPARTICLE_H

#include "particles/Particle.h"

class SPHParticle : public autopas::Particle {
public:
    SPHParticle() : autopas::Particle() {}

    virtual ~SPHParticle() {}

    SPHParticle(std::array<double, 3> r, std::array<double, 3> v, unsigned long id) : autopas::Particle(r, v, id) {}




    double getDensity() const {
        return _density;
    }

    double addDensity(double density) {
        _density += density;
    }

    void setDensity(double density) {
        _density = density;
    }

    double getPressure() const {
        return _pressure;
    }

    void setPressure(double pressure) {
        _pressure = pressure;
    }

    double getMass() const {
        return _mass;
    }

    void setMass(double mass) {
        _mass = mass;
    }

    double getSmth() const {
        return _smth;
    }

    void setSmth(double smth) {
        _smth = smth;
    }

private:
    double _density;
    double _pressure;
    double _mass;
    double _smth;
};


#endif //AUTOPAS_SPHPARTICLE_H

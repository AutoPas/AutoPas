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

    double getSnds() const {
        return _snds;
    }

    void setSnds(double snds) {
        _snds = snds;
    }

    double getVSigMax() const {
        return _v_sig_max;
    }

    void checkAndSetVSigMax(double v_sig) {
        _v_sig_max = std::max(v_sig, _v_sig_max);
    }

    void setVSigMax(double v_sig_max) {
        _v_sig_max = v_sig_max;
    }

    const std::array<double, 3>& getAcc() const {
        return _acc;
    }

    void addAcc(const std::array<double, 3>& acc);

    void subAcc(const std::array<double, 3>& acc);

    void setAcc(const std::array<double, 3>& acc) {
        _acc = acc;
    }

    double getEngDot() const {
        return _eng_dot;
    }

    void addEngDot(double eng_dot) {
        _eng_dot += eng_dot;
    }

    void setEngDot(double eng_dot) {
        _eng_dot = eng_dot;
    }

private:
    double _density;
    double _pressure;
    double _mass;
    double _smth;
    double _snds;

    //temporaries / helpers
    double _v_sig_max;
    std::array<double, 3> _acc;
    double _eng_dot;
};


#endif //AUTOPAS_SPHPARTICLE_H

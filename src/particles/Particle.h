/*
 * Particle.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLE_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLE_H_

#include "utils/arrayMath.h"
#include <array>

namespace autopas {

    class Particle {
    public:
        Particle() : _r({0.0, 0.0, 0.0}), _v({0., 0., 0.}), _f({0.0, 0.0, 0.0}), _id(0) {}

        Particle(std::array<double, 3> r, std::array<double, 3> v, unsigned long id) : _r(r), _v(v),
                                                                                       _f({0.0, 0.0, 0.0}), _id(id) {}

        virtual ~Particle() {}

        const std::array<double, 3> &getF() const { return _f; }

        void setF(const std::array<double, 3> &f) { _f = f; }

        void addF(const std::array<double, 3> &f) { _f = arrayMath::add(_f, f); }

        void subF(const std::array<double, 3> &f) { _f = arrayMath::sub(_f, f); }

        unsigned long getID() const { return _id; }

        void setID(unsigned long id) { _id = id; }

        const std::array<double, 3> &getR() const { return _r; }

        void setR(const std::array<double, 3> &r) { _r = r; }

        const std::array<double, 3> &getV() const { return _v; }

        void setV(const std::array<double, 3> &v) { _v = v; }

    private:
        std::array<double, 3> _r;
        std::array<double, 3> _v;
        std::array<double, 3> _f;
        unsigned long _id;
    };

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLE_H_ */

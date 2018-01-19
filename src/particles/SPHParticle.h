//
// Created by seckler on 19.01.18.
//

#ifndef AUTOPAS_SPHPARTICLE_H
#define AUTOPAS_SPHPARTICLE_H

#include "Particle.h"

namespace autopas {

    class SPHParticle : public Particle {
    public:
        SPHParticle() : Particle() {}

        virtual ~SPHParticle() {}

        SPHParticle(std::array<double, 3> r, std::array<double, 3> v, unsigned long id) : Particle(r, v, id) {}


    private:
        

        double density;
        double pressure;
        double surf_norm;
    };
}  // namespace autopas

#endif //AUTOPAS_SPHPARTICLE_H

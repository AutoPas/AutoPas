//
// Created by seckler on 19.01.18.
//

#include "autopas.h"
#include "sph/autopassph.h"
#include <array>
#include <iostream>
#include <cassert>
void
SetupIC(autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> &sph_system,
        double *end_time, const std::array<double, 3> &bBoxMax) {
    // Place SPH particles

    const double dx = 1.0 / 128.0;
    unsigned int i = 0;
    for (double x = 0; x < bBoxMax[0] * 0.5; x += dx) {
        for (double y = 0; y < bBoxMax[1]; y += dx) {
            for (double z = 0; z < bBoxMax[2]; z += dx) {
                //std::array<double, 3> r, std::array<double, 3> v, unsigned long id, double mass, double smth, double snds
                autopas::sph::SPHParticle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012, 0.);
                ith.setDensity(1.0);
                //ith.pos.x = x;
                //ith.pos.y = y;
                //ith.pos.z = z;
                //ith.dens = 1.0;
                //ith.mass = 0.75;
                //ith.eng = 2.5;
                //ith.id = i++;
                //ith.smth = 0.012;
                sph_system.addParticle(ith);
            }
        }
    }
    for (double x = bBoxMax[0] * 0.5; x < bBoxMax[0] * 1.; x += dx * 2.0) {
        for (double y = 0; y < bBoxMax[1]; y += dx) {
            for (double z = 0; z < bBoxMax[2]; z += dx) {
                //std::array<double, 3> r, std::array<double, 3> v, unsigned long id, double mass, double smth, double snds
                autopas::sph::SPHParticle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012, 0.);
                ith.setDensity(0.5);
                //ith.pos.x = x;
                //ith.pos.y = y;
                //ith.pos.z = z;
                //ith.dens = 0.5;
                //ith.mass = 0.75;
                //ith.eng = 2.5;
                //ith.id = i++;
                //ith.smth = 0.012;
                sph_system.addParticle(ith);
            }
        }
    }
    for (auto part = sph_system.begin(); part.isValid(); ++part) {
        part->setMass(part->getMass() * bBoxMax[0] * bBoxMax[1] * bBoxMax[2] / (double) (i));
    }
    std::cout << "# of ptcls is... " << i << std::endl;

    // Set the end time
    *end_time = 0.12;
    // Fin.
    std::cout << "setup..." << std::endl;
}



int main() {
    std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
    boxMax[0] = 1.0;
    boxMax[1] = boxMax[2] = boxMax[0] / 8.0;
    double cutoff = 3.;
    autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> sphSystem(
            boxMin, boxMax, cutoff);
    double t_end;
    SetupIC(sphSystem, &t_end, boxMax);

}
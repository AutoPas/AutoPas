//
// Created by newc on 21.11.24.
//

#pragma once

#include "MultisiteMoleculeLJ.h"

namespace mdLib {
   /**
    * Variant of the MultisiteMoleculeLJ which stores its exact site positions instead of using a quaternion in
    * combination
    */
   class MultisiteMolLJExactSitePosVar : public MultisiteMoleculeLJ {
   MultisiteMolLJExactSitePosVar() = default;

   /**
    * Constructor of the MultisiteMolLJExactSitePosVar Class
    * @param r Position of the particle.
    * @param v Velocity of the particle.
    * @param exactSitePositions positions of sites of the particles
    * @param angularVel Rotational velocity of the particle.
    * @param moleculeId Id of the particle.
    * @param typeId Id of the type of the particle. Used in conjunction with ParticlePropertiesLibrary to access
    * molecular information such as site types and relative site positions.
    */
   MultisiteMolLJExactSitePosVar(std::array<double, 3> r, std::array<double, 3> v, std::vector<std::array<double, 3>> exactSitePositions,
                       std::array<double, 3> angularVel, unsigned long moleculeId, unsigned long typeId = 0);

   /**
    * Destructor of the MultisiteMoleculeLJ class.
    */
   ~MultisiteMolLJExactSitePosVar() override = default;
};
}




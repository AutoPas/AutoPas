/**
 * @file ParticleAttributes.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include <array>

#include "autopas/particles/OwnershipState.h"


/*
* A struct containing all properties of autopas::MoleculeLJ<double>.
* This can be used to align the attributes of a particle in memory to make serialization and deserialization easier.
*/
struct ParticleAttributes {
	// ParticleBase attributes
  std::array<double, 3> position;
  std::array<double, 3> velocity;
  std::array<double, 3> force;
	unsigned long id; 		
	autopas::OwnershipState ownershipState;

	// MoleculeLJ attributes
	size_t typeId;
  std::array<double, 3> oldForce;
};

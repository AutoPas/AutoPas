/**
 * @file ParticleAttributes.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

/*
* A struct containing all properties of autopas::MoleculeLJ<double>.
* This can be used to align the attributes of a particle in memory to make serialization and deserialization easier.
*/
struct ParticleAttributes {
	// ParticleBase attributes
	double positionX; 		
	double positionY;
	double positionZ;
	double velocityX; 		
	double velocityY;
	double velocityZ;
	double forceX;
	double forceY;
	double forceZ;
	unsigned long id; 		

	// MoleculeLJ attributes
	size_t typeId;
	double oldForceX;
	double oldForceY;
	double oldForceZ;
};

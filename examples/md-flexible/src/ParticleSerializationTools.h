/**
 * @file ParticleSerializationTools.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include "TypeDefinitions.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"

#include <vector>

/* 
* Provides tools to de-/serialize particles of type autopas::MoleculeLJ<double>.
*/
namespace ParticleSerializationTools {
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

	/* 
	* Serializes a particle and appends it to the serializedParticles container.
	* @param particle The particle which will be serialized.
	*	@param serializedParticles The container to wich the serialized particle will be appended.
	*/
	void serializeParticle(ParticleType &particle, std::vector<char> &serializedParticles);

	/*
	* Deserializes a serialized particle.
	* @param particleData A pointer to the serialized particle data.
	* @param particle The particle to which the desierialized attributes will be applied. 
	*/
	void deserializeParticleData(char* particleData, ParticleType &particle);

	/*
	* Deserializes a container of serialized particles.
	* @param particlesData A pointer to the serialized particle data.
	* @param particles The particle container to which to append the deserialized particles to. 
	*/
	void deserializeParticleData(std::vector<char> &particlesData, std::vector<ParticleType> &particles);
}



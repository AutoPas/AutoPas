/**
 * @file ParticleSerializationTools.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "ParticleAttributes.h"
#include "TypeDefinitions.h"

#include <vector>

/* 
 * Provides tools to de-/serialize particles of type autopas::MoleculeLJ<double>.
 */
namespace ParticleSerializationTools {
	/** 
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

	/**
	 * Deserializes a container of serialized particles.
	 * @param particlesData A pointer to the serialized particle data.
	 * @param particles The particle container to which to append the deserialized particles to. 
	 */
	void deserializeParticleData(std::vector<char> &particlesData, std::vector<ParticleType> &particles);

	/**
	 * Converts ParticleAttributes to a particle of type MoleculeLJ<double>.
	 * @param particleAttributes The attributes which will be applied to the returne particle.
	 * @return returns a particle of type MoleculeLJ<double>
	 */
	ParticleType convertParticleAttributesToParticle(ParticleAttributes &attributes);
}



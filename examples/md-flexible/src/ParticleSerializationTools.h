/**
 * @file ParticleSerializationTools.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "TypeDefinitions.h"

/**
 * Provides tools to de-/serialize particles of type autopas::MoleculeLJ<double>.
 */
namespace ParticleSerializationTools {
/**
 * Serializes a particle and appends it to the serializedParticles container.
 * @param particle The particle which will be serialized.
 * @param serializedParticles The container to wich the serialized particle will be appended.
 */
template <class ParticleClass>
void serializeParticle(const ParticleClass &particle, std::vector<char> &serializedParticles);

/**
 * Deserializes a serialized particle.
 * @param particleData A pointer to the serialized particle data.
 * @param particle The particle to which the desierialized attributes will be applied.
 */
template <class ParticleClass>
void deserializeParticle(char *particleData, ParticleClass &particle);

/**
 * Deserializes a container of serialized particles.
 * @param particlesData A pointer to the serialized particle data.
 * @param particles The particle container to which to append the deserialized particles to.
 */
template <class ParticleClass>
void deserializeParticles(std::vector<char> &particlesData, std::vector<ParticleClass> &particles);
}  // namespace ParticleSerializationTools

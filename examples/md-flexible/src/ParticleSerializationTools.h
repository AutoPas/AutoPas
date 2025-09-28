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
 * The combined size in byte of the simple attributes which need to be communicated using MPI.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
constexpr size_t AttributesSize = 200;
#else
constexpr size_t AttributesSize = 120;
#endif

/**
 * Serializes a particle and appends it to the serializedParticles container.
 * @param particle The particle which will be serialized.
 * @param serializedParticles The container to wich the serialized particle will be appended.
 */
void serializeParticle(const ParticleType &particle, std::vector<char> &serializedParticles);

/**
 * Deserializes a serialized particle.
 * @param particleData A pointer to the serialized particle data.
 * @param particle The particle to which the deserialized attributes will be applied.
 */
void deserializeParticle(char *particleData, ParticleType &particle);

/**
 * Deserializes a container of serialized particles.
 * @param particlesData A pointer to the serialized particle data.
 * @param particles The particle container to which to append the deserialized particles to.
 */
void deserializeParticles(std::vector<char> &particlesData, std::vector<ParticleType> &particles);
}  // namespace ParticleSerializationTools

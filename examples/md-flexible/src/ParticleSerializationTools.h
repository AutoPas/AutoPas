/**
 * @file ParticleSerializationTools.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "ParticleAttributes.h"
#include "TypeDefinitions.h"

/*
 * Provides tools to de-/serialize particles of type autopas::MoleculeLJ<double>.
 */
namespace ParticleSerializationTools {
//using AttributeSequence = std::integer_sequence<size_t, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15>;
//enum AttributeNames { id = 1, posX = 2, posY = 3, posZ = 4, velocityX = 5, velocityY = 6, velocityZ = 7, forceX = 8, forceY = 9, forceZ = 10, oldForceX = 11, oldForceY = 12, oldForceZ = 13, typeId = 14, ownershipState = 15};
//constexpr static std::array<AttributeNames, 0> getNeededAttr();

//constexpr static auto Attributes = std::make_tuple(ParticleType::AttributeNames::id, ParticleType::AttributeNames::posX, ParticleType::AttributeNames::posY, ParticleType::AttributeNames::posZ, ParticleType::AttributeNames::velocityX, ParticleType::AttributeNames::velocityY, ParticleType::AttributeNames::velocityZ, ParticleType::AttributeNames::forceX, ParticleType::AttributeNames::forceY, ParticleType::AttributeNames::forceZ, ParticleType::AttributeNames::oldForceX, ParticleType::AttributeNames::oldForceY, ParticleType::AttributeNames::oldForceZ, ParticleType::AttributeNames::typeId, ParticleType::AttributeNames::ownershipState);

/**
 * Serializes a particle and appends it to the serializedParticles container.
 * @param particle The particle which will be serialized.
 * @param serializedParticles The container to wich the serialized particle will be appended.
 */
void serializeParticle(const ParticleType &particle, std::vector<char> &serializedParticles);

/*
 * Deserializes a serialized particle.
 * @param particleData A pointer to the serialized particle data.
 * @param particle The particle to which the desierialized attributes will be applied.
 */
void deserializeParticle(char *particleData, ParticleType &particle);

/**
 * Deserializes a container of serialized particles.
 * @param particlesData A pointer to the serialized particle data.
 * @param particles The particle container to which to append the deserialized particles to.
 */
void deserializeParticles(std::vector<char> &particlesData, std::vector<ParticleType> &particles);
}  // namespace ParticleSerializationTools

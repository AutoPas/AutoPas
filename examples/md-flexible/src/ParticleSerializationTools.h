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

// todo implement this for rotational molecules

namespace {

/**
 * Stores the AttributeNames of the attributes of ParticleType which have to be communicated using MPI.
 */
constexpr std::array<typename autopas::MoleculeLJ::AttributeNames, 15> SimpleAttributes = {
    autopas::MoleculeLJ::AttributeNames::id,
    autopas::MoleculeLJ::AttributeNames::posX,
    autopas::MoleculeLJ::AttributeNames::posY,
    autopas::MoleculeLJ::AttributeNames::posZ,
    autopas::MoleculeLJ::AttributeNames::velocityX,
    autopas::MoleculeLJ::AttributeNames::velocityY,
    autopas::MoleculeLJ::AttributeNames::velocityZ,
    autopas::MoleculeLJ::AttributeNames::forceX,
    autopas::MoleculeLJ::AttributeNames::forceY,
    autopas::MoleculeLJ::AttributeNames::forceZ,
    autopas::MoleculeLJ::AttributeNames::oldForceX,
    autopas::MoleculeLJ::AttributeNames::oldForceY,
    autopas::MoleculeLJ::AttributeNames::oldForceZ,
    autopas::MoleculeLJ::AttributeNames::typeId,
    autopas::MoleculeLJ::AttributeNames::ownershipState};

/**
 * Stores the AttributeNames of the attributes of ParticleType which have to be communicated using MPI.
 */
constexpr std::array<typename autopas::MulticenteredMoleculeLJ::AttributeNames, 25> RotationalAttributes = {
    autopas::MulticenteredMoleculeLJ::AttributeNames::id,
    autopas::MulticenteredMoleculeLJ::AttributeNames::posX,
    autopas::MulticenteredMoleculeLJ::AttributeNames::posY,
    autopas::MulticenteredMoleculeLJ::AttributeNames::posZ,
    autopas::MulticenteredMoleculeLJ::AttributeNames::velocityX,
    autopas::MulticenteredMoleculeLJ::AttributeNames::velocityY,
    autopas::MulticenteredMoleculeLJ::AttributeNames::velocityZ,
    autopas::MulticenteredMoleculeLJ::AttributeNames::forceX,
    autopas::MulticenteredMoleculeLJ::AttributeNames::forceY,
    autopas::MulticenteredMoleculeLJ::AttributeNames::forceZ,
    autopas::MulticenteredMoleculeLJ::AttributeNames::oldForceX,
    autopas::MulticenteredMoleculeLJ::AttributeNames::oldForceY,
    autopas::MulticenteredMoleculeLJ::AttributeNames::oldForceZ,
    autopas::MulticenteredMoleculeLJ::AttributeNames::quaternion0,
    autopas::MulticenteredMoleculeLJ::AttributeNames::quaternion1,
    autopas::MulticenteredMoleculeLJ::AttributeNames::quaternion2,
    autopas::MulticenteredMoleculeLJ::AttributeNames::quaternion3,
    autopas::MulticenteredMoleculeLJ::AttributeNames::angularVelX,
    autopas::MulticenteredMoleculeLJ::AttributeNames::angularVelY,
    autopas::MulticenteredMoleculeLJ::AttributeNames::angularVelZ,
    autopas::MulticenteredMoleculeLJ::AttributeNames::torqueX,
    autopas::MulticenteredMoleculeLJ::AttributeNames::torqueY,
    autopas::MulticenteredMoleculeLJ::AttributeNames::torqueZ,
    autopas::MulticenteredMoleculeLJ::AttributeNames::typeId,
    autopas::MulticenteredMoleculeLJ::AttributeNames::ownershipState};

/**
 * The combined size in byte of the simple attributes which need to be communicated using MPI.
 */
constexpr size_t simpleAttributesSize = 120;

/**
 * The combined size in byte of the rotational attributes which need to be communicated using MPI.
 */
constexpr size_t rotationalAttributesSize = 200;

/**
 * Serializes the attribute of simple molecule defined by I.
 * @param particle: The particle who's attribute needs to be serialized.
 * @param attributeVector: The container in which the serialized attribute will be stored.
 * @param startIndex: The startindex in the container where to store the serialized attribute.
 */
template <size_t I>
void serializeSimpleAttribute(const autopas::MoleculeLJ &particle, std::vector<char> &attributeVector, size_t &startIndex) {
  const auto attribute = particle.get<SimpleAttributes[I]>();
  const auto sizeOfValue = sizeof(attribute);
  std::memcpy(&attributeVector[startIndex], &attribute, sizeOfValue);
  startIndex += sizeOfValue;
}

/**
 * Deserializes the attribute of simple molecule defined by I.
 * @param attributeVector: The vector containing the data which needs to be deserialized.
 * @param particle: The particle to which the serialized data will be applied.
 * @param startIndex: The start index in the attributeVector of the attribute which needs to be deserialized.
 */
template <size_t I>
void deserializeSimpleAttribute(char *&attributeVector, autopas::MoleculeLJ &particle, size_t &startIndex) {
  auto attribute = particle.get<SimpleAttributes[I]>();
  const auto sizeOfValue = sizeof(attribute);
  std::memcpy(&attribute, &attributeVector[startIndex], sizeOfValue);
  particle.set<SimpleAttributes[I]>(attribute);
  startIndex += sizeOfValue;
}

/**
 * The implementation of serializeParticle using the expansion operator.
 * @param particle: The particle which will be serialized.
 * @param serializedParticle: The char array of the particles serialized attributes.
 */
template <size_t... I>
void serializeSimpleParticleImpl(const autopas::MoleculeLJ &particle, std::vector<char> &serializedParticle,
                                 std::index_sequence<I...>) {
  // Serialize particle attributes
  size_t startIndex = 0;
  std::vector<char> attributesVector(simpleAttributesSize);
  (serializeSimpleAttribute<I>(particle, attributesVector, startIndex), ...);

  // Add serialized attributes to serialized particle
  serializedParticle.insert(serializedParticle.end(), attributesVector.begin(), attributesVector.end());
}

/**
 * The implementation fo deserializeParticle using the expansion operator.
 * @param particleData: The particle data which will be deserialized.
 * @param particle: The particle to which the deserialized attributes will be applied.
 */
template <size_t... I>
void deserializeSimpleParticleImpl(char *particleData, autopas::MoleculeLJ &particle, std::index_sequence<I...>) {
  size_t startIndex = 0;
  (deserializeSimpleAttribute<I>(particleData, particle, startIndex), ...);
}
}  // namespace

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
void serializeParticle(const ParticleClass &particle, std::vector<char> &serializedParticles) {
  autopas::utils::ExceptionHandler::exception("ParticleSerializationTools not implemented for particle");
}

/**
 * Deserializes a serialized particle.
 * @param particleData A pointer to the serialized particle data.
 * @param particle The particle to which the desierialized attributes will be applied.
 */
template <class ParticleClass>
void deserializeParticle(char *particleData, ParticleClass &particle) {
  autopas::utils::ExceptionHandler::exception("ParticleSerializationTools not implemented for particle");
}

/**
 * Deserializes a container of serialized particles.
 * @param particlesData A pointer to the serialized particle data.
 * @param particles The particle container to which to append the deserialized particles to.
 */
template <class ParticleClass>
void deserializeParticles(std::vector<char> &particlesData, std::vector<ParticleClass> &particles) {
  autopas::utils::ExceptionHandler::exception("ParticleSerializationTools not implemented for particle");
}
}  // namespace ParticleSerializationTools

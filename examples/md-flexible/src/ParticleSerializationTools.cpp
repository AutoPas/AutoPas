/**
 * @file ParticleSerializationTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "ParticleSerializationTools.h"

#include <tuple>

namespace {

/**
 * Stores the AttributeNames of the attributes of ParticleType which have to be communicated using MPI.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
constexpr std::array<typename ParticleType::AttributeNames, 25> Attributes = {
    mdLib::MultisiteMoleculeLJ::AttributeNames::id,
    mdLib::MultisiteMoleculeLJ::AttributeNames::posX,
    mdLib::MultisiteMoleculeLJ::AttributeNames::posY,
    mdLib::MultisiteMoleculeLJ::AttributeNames::posZ,
    mdLib::MultisiteMoleculeLJ::AttributeNames::velocityX,
    mdLib::MultisiteMoleculeLJ::AttributeNames::velocityY,
    mdLib::MultisiteMoleculeLJ::AttributeNames::velocityZ,
    mdLib::MultisiteMoleculeLJ::AttributeNames::forceX,
    mdLib::MultisiteMoleculeLJ::AttributeNames::forceY,
    mdLib::MultisiteMoleculeLJ::AttributeNames::forceZ,
    mdLib::MultisiteMoleculeLJ::AttributeNames::oldForceX,
    mdLib::MultisiteMoleculeLJ::AttributeNames::oldForceY,
    mdLib::MultisiteMoleculeLJ::AttributeNames::oldForceZ,
    mdLib::MultisiteMoleculeLJ::AttributeNames::quaternion0,
    mdLib::MultisiteMoleculeLJ::AttributeNames::quaternion1,
    mdLib::MultisiteMoleculeLJ::AttributeNames::quaternion2,
    mdLib::MultisiteMoleculeLJ::AttributeNames::quaternion3,
    mdLib::MultisiteMoleculeLJ::AttributeNames::angularVelX,
    mdLib::MultisiteMoleculeLJ::AttributeNames::angularVelY,
    mdLib::MultisiteMoleculeLJ::AttributeNames::angularVelZ,
    mdLib::MultisiteMoleculeLJ::AttributeNames::torqueX,
    mdLib::MultisiteMoleculeLJ::AttributeNames::torqueY,
    mdLib::MultisiteMoleculeLJ::AttributeNames::torqueZ,
    mdLib::MultisiteMoleculeLJ::AttributeNames::typeId,
    mdLib::MultisiteMoleculeLJ::AttributeNames::ownershipState};
#elif defined(MD_FLEXIBLE_FUNCTOR_DEM)
constexpr std::array<typename ParticleType::AttributeNames, 21> Attributes = {
    demLib::GranularDEMParticle::AttributeNames::id,
    demLib::GranularDEMParticle::AttributeNames::posX,
    demLib::GranularDEMParticle::AttributeNames::posY,
    demLib::GranularDEMParticle::AttributeNames::posZ,
    demLib::GranularDEMParticle::AttributeNames::velocityX,
    demLib::GranularDEMParticle::AttributeNames::velocityY,
    demLib::GranularDEMParticle::AttributeNames::velocityZ,
    demLib::GranularDEMParticle::AttributeNames::forceX,
    demLib::GranularDEMParticle::AttributeNames::forceY,
    demLib::GranularDEMParticle::AttributeNames::forceZ,
    demLib::GranularDEMParticle::AttributeNames::oldForceX,
    demLib::GranularDEMParticle::AttributeNames::oldForceY,
    demLib::GranularDEMParticle::AttributeNames::oldForceZ,
    demLib::GranularDEMParticle::AttributeNames::angularVelX,
    demLib::GranularDEMParticle::AttributeNames::angularVelY,
    demLib::GranularDEMParticle::AttributeNames::angularVelZ,
    demLib::GranularDEMParticle::AttributeNames::torqueX,
    demLib::GranularDEMParticle::AttributeNames::torqueY,
    demLib::GranularDEMParticle::AttributeNames::torqueZ,
    demLib::GranularDEMParticle::AttributeNames::typeId,
    demLib::GranularDEMParticle::AttributeNames::ownershipState};
#else
constexpr std::array<typename ParticleType::AttributeNames, 15> Attributes = {
    mdLib::MoleculeLJ::AttributeNames::id,
    mdLib::MoleculeLJ::AttributeNames::posX,
    mdLib::MoleculeLJ::AttributeNames::posY,
    mdLib::MoleculeLJ::AttributeNames::posZ,
    mdLib::MoleculeLJ::AttributeNames::velocityX,
    mdLib::MoleculeLJ::AttributeNames::velocityY,
    mdLib::MoleculeLJ::AttributeNames::velocityZ,
    mdLib::MoleculeLJ::AttributeNames::forceX,
    mdLib::MoleculeLJ::AttributeNames::forceY,
    mdLib::MoleculeLJ::AttributeNames::forceZ,
    mdLib::MoleculeLJ::AttributeNames::oldForceX,
    mdLib::MoleculeLJ::AttributeNames::oldForceY,
    mdLib::MoleculeLJ::AttributeNames::oldForceZ,
    mdLib::MoleculeLJ::AttributeNames::typeId,
    mdLib::MoleculeLJ::AttributeNames::ownershipState};
#endif

/**
 * The combined size in byte of the simple attributes which need to be communicated using MPI.
 */
#if MD_FLEXIBLE_MODE == MULTISITE
constexpr size_t AttributesSize = 200;
#else
constexpr size_t AttributesSize = 120;
#endif

/**
 * Serializes the attribute of a molecule defined by I.
 * @param particle: The particle who's attribute needs to be serialized.
 * @param attributeVector: The container in which the serialized attribute will be stored.
 * @param startIndex: The startindex in the container where to store the serialized attribute.
 */
template <size_t I>
void serializeAttribute(const ParticleType &particle, std::vector<char> &attributeVector, size_t &startIndex) {
  const auto attribute = particle.get<Attributes[I]>();
  const auto sizeOfValue = sizeof(attribute);
  std::memcpy(&attributeVector[startIndex], &attribute, sizeOfValue);
  startIndex += sizeOfValue;
}

/**
 * Deserializes the attribute defined by I.
 * @param attributeVector: The vector containing the data which needs to be deserialized.
 * @param particle: The particle to which the serialized data will be applied.
 * @param startIndex: The start index in the attributeVector of the attribute which needs to be deserialized.
 */
template <size_t I>
void deserializeAttribute(char *&attributeVector, ParticleType &particle, size_t &startIndex) {
  auto attribute = particle.get<Attributes[I]>();
  const auto sizeOfValue = sizeof(attribute);
  std::memcpy(&attribute, &attributeVector[startIndex], sizeOfValue);
  particle.set<Attributes[I]>(attribute);
  startIndex += sizeOfValue;
}

/**
 * The implementation of serializeParticle using the expansion operator.
 * @param particle: The particle which will be serialized.
 * @param serializedParticle: The char array of the particles serialized attributes.
 */
template <size_t... I>
void serializeParticleImpl(const ParticleType &particle, std::vector<char> &serializedParticle,
                           std::index_sequence<I...>) {
  // Serialize particle attributes
  size_t startIndex = 0;
  std::vector<char> attributesVector(AttributesSize);
  (serializeAttribute<I>(particle, attributesVector, startIndex), ...);

  // Add serialized attributes to serialized particle
  serializedParticle.insert(serializedParticle.end(), attributesVector.begin(), attributesVector.end());
}

/**
 * The implementation of deserializeParticle using the expansion operator.
 * @param particleData: The particle data which will be deserialized.
 * @param particle: The particle to which the deserialized attributes will be applied.
 */
template <size_t... I>
void deserializeParticleImpl(char *particleData, ParticleType &particle, std::index_sequence<I...>) {
  size_t startIndex = 0;
  (deserializeAttribute<I>(particleData, particle, startIndex), ...);
}
}  // namespace

namespace ParticleSerializationTools {
void serializeParticle(const ParticleType &particle, std::vector<char> &serializedParticles) {
  serializeParticleImpl(particle, serializedParticles, std::make_index_sequence<Attributes.size()>{});
}

void deserializeParticle(char *particleData, ParticleType &particle) {
  deserializeParticleImpl(particleData, particle, std::make_index_sequence<Attributes.size()>{});
}

void deserializeParticles(std::vector<char> &particlesData, std::vector<ParticleType> &particles) {
  ParticleType particle;
  // Reserve space for particles before deserializing them
  const auto numParticles = particlesData.size() / AttributesSize;
  particles.reserve(particles.size() + numParticles);
  for (size_t i = 0; i < particlesData.size(); i += AttributesSize) {
    deserializeParticle(&particlesData[i], particle);
    particles.push_back(particle);
  }
}

}  // namespace ParticleSerializationTools

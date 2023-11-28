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
#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
constexpr std::array<typename ParticleType::AttributeNames, 24> Attributes = {
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::id,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::posX,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::posY,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::posZ,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::velocityX,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::velocityY,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::velocityZ,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::forceX,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::forceY,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::forceZ,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::oldForceX,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::oldForceY,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::oldForceZ,
    //mdLib::PositionStoringMultiSiteMolecule::AttributeNames::quaternion0,
    //mdLib::PositionStoringMultiSiteMolecule::AttributeNames::quaternion1,
    //mdLib::PositionStoringMultiSiteMolecule::AttributeNames::quaternion2,
    //mdLib::PositionStoringMultiSiteMolecule::AttributeNames::quaternion3,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::relativeSitePositionsX,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::relativeSitePositionsY,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::relativeSitePositionsZ,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::angularVelX,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::angularVelY,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::angularVelZ,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::torqueX,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::torqueY,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::torqueZ,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::typeId,
    mdLib::PositionStoringMultiSiteMolecule::AttributeNames::ownershipState};
#elif defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
constexpr std::array<typename ParticleType::AttributeNames, 10> Attributes ={
                                                                mdLib::Site::AttributeNames::id,
                                                                mdLib::Site::AttributeNames::owningMoleculeId,
                                                                mdLib::Site::AttributeNames::posX,
                                                                mdLib::Site::AttributeNames::posY,
                                                                mdLib::Site::AttributeNames::posZ,
                                                                mdLib::Site::AttributeNames::forceX,
                                                                mdLib::Site::AttributeNames::forceY,
                                                                mdLib::Site::AttributeNames::forceZ,
                                                                mdLib::Site::AttributeNames::typeId,
                                                                mdLib::Site::AttributeNames::ownershipState};
#else
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
#endif
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
 //@TODO: (Johnny) I am almost certain i need to adapt this value in some cases! Not entirely sure how though
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
  for (size_t i = 0; i < particlesData.size(); i += AttributesSize) {
    deserializeParticle(&particlesData[i], particle);
    particles.push_back(particle);
  }
}

}  // namespace ParticleSerializationTools

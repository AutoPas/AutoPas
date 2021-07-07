/**
 * @file ParticleSerializationTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "ParticleSerializationTools.h"

namespace {

constexpr static std::array<typename ParticleType::AttributeNames, 15> Attributes = { ParticleType::AttributeNames::id, ParticleType::AttributeNames::posX, ParticleType::AttributeNames::posY, ParticleType::AttributeNames::posZ, ParticleType::AttributeNames::velocityX, ParticleType::AttributeNames::velocityY, ParticleType::AttributeNames::velocityZ, ParticleType::AttributeNames::forceX, ParticleType::AttributeNames::forceY, ParticleType::AttributeNames::forceZ, ParticleType::AttributeNames::oldForceX, ParticleType::AttributeNames::oldForceY, ParticleType::AttributeNames::oldForceZ, ParticleType::AttributeNames::typeId, ParticleType::AttributeNames::ownershipState };

template<size_t I>
void serializeAttribute(const ParticleType &particle, std::vector<char> &attributesVector) {
  const auto startIndex = attributesVector.size();
  constexpr auto attributeType = Attributes[I];
  const auto value = particle.get<attributeType>();
  const auto sizeOfValue = sizeof(value);
  attributesVector.resize(startIndex + sizeOfValue);
  std::memcpy(&attributesVector[startIndex], &value, sizeOfValue);
}

template<size_t... I>
void serializeParticleImpl(const ParticleType &particle, std::vector<char> &serializedParticles, std::index_sequence<I...>){
  std::vector<char> attributesVector;

  (serializeAttribute<I>(particle, attributesVector), ...);

  serializedParticles.insert(serializedParticles.end(), attributesVector.begin(), attributesVector.end());
}
}

namespace ParticleSerializationTools {

//void serializeParticle(const ParticleType &particle, std::vector<char> &serializedParticles) {
//  ParticleAttributes attributes{.position{particle.getR()},
//                                .velocity{particle.getV()},
//                                .force{particle.getF()},
//                                .id = particle.getID(),
//                                .ownershipState = particle.getOwnershipState(),
//                                .typeId = particle.getTypeId(),
//                                .oldForce{particle.getOldF()}};
//
//  std::vector<char> attributesVector;
//  attributesVector.resize(sizeof(ParticleAttributes));
//
//  std::memcpy(&attributesVector[0], &attributes, sizeof(ParticleAttributes));
//
//  serializedParticles.insert(serializedParticles.end(), attributesVector.begin(), attributesVector.end());
//}

void serializeParticle(const ParticleType &particle, std::vector<char> &serializedParticles) {
  serializeParticleImpl(particle, serializedParticles, std::make_index_sequence<Attributes.size()>{});
}

void deserializeParticle(char *particleData, ParticleType &particle) {
  ParticleAttributes attributes;
  std::memcpy(&attributes, particleData, sizeof(ParticleAttributes));

  particle.setR(attributes.position);
  particle.setV(attributes.velocity);
  particle.setF(attributes.force);
  particle.setOldF(attributes.oldForce);
  particle.setID(attributes.id);
  particle.setOwnershipState(attributes.ownershipState);
  particle.setTypeId(attributes.typeId);
}

void deserializeParticles(std::vector<char> &particlesData, std::vector<ParticleType> &particles) {
  ParticleType particle;
  size_t sizeOfParticleAttributes = sizeof(ParticleAttributes);
  for (size_t i = 0; i < particlesData.size(); i += sizeOfParticleAttributes) {
    deserializeParticle(&particlesData[i], particle);
    particles.push_back(particle);
  }
}

}  // namespace ParticleSerializationTools

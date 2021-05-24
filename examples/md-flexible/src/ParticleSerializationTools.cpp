/**
 * @file ParticleSerializationTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "ParticleSerializationTools.h"

namespace ParticleSerializationTools {
void serializeParticle(ParticleType &particle, std::vector<char> &serializedParticles) {
  ParticleAttributes attributes;

  attributes.position = particle.getR();
  attributes.velocity = particle.getV();
  attributes.force = particle.getF();

  if (particle.isDummy()) {
    attributes.ownershipState = autopas::OwnershipState::dummy;
  } else if (particle.isOwned()) {
    attributes.ownershipState = autopas::OwnershipState::owned;
  } else if (particle.isHalo()) {
    attributes.ownershipState = autopas::OwnershipState::halo;
  }

  attributes.id = particle.getID();
  attributes.typeId = particle.getTypeId();
  attributes.oldForce = particle.getOldF();

  std::vector<char> attributesVector;
  attributesVector.resize(sizeof(ParticleAttributes));

  std::memcpy(&attributesVector[0], &attributes, sizeof(ParticleAttributes));

  ParticleType deserializedParticle;
  deserializeParticleData(attributesVector.data(), deserializedParticle);

  serializedParticles.insert(serializedParticles.end(), attributesVector.begin(), attributesVector.end());
}

void deserializeParticleData(char *particleData, ParticleType &particle) {
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

void deserializeParticleData(std::vector<char> &particlesData, std::vector<ParticleType> &particles) {
  ParticleType particle;
  size_t sizeOfParticleAttributes = sizeof(ParticleAttributes);
  for (size_t i = 0; i < particlesData.size(); i += sizeOfParticleAttributes) {
    deserializeParticleData(&particlesData[i], particle);
    particles.push_back(particle);
  }
}

ParticleType convertParticleAttributesToParticle(ParticleAttributes &attributes) {
  ParticleType particle;

  particle.setR(attributes.position);
  particle.setV(attributes.velocity);
  particle.setF(attributes.force);
  particle.setOldF(attributes.oldForce);
  particle.setID(attributes.id);
  particle.setTypeId(attributes.typeId);

  return particle;
}
}  // namespace ParticleSerializationTools

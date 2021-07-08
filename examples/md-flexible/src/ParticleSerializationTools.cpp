/**
 * @file ParticleSerializationTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "ParticleSerializationTools.h"

#include <tuple>

namespace {

constexpr std::array<typename ParticleType::AttributeNames, 15> Attributes = {
    ParticleType::AttributeNames::id,
    ParticleType::AttributeNames::posX,
    ParticleType::AttributeNames::posY,
    ParticleType::AttributeNames::posZ,
    ParticleType::AttributeNames::velocityX,
    ParticleType::AttributeNames::velocityY,
    ParticleType::AttributeNames::velocityZ,
    ParticleType::AttributeNames::forceX,
    ParticleType::AttributeNames::forceY,
    ParticleType::AttributeNames::forceZ,
    ParticleType::AttributeNames::oldForceX,
    ParticleType::AttributeNames::oldForceY,
    ParticleType::AttributeNames::oldForceZ,
    ParticleType::AttributeNames::typeId,
    ParticleType::AttributeNames::ownershipState};

constexpr size_t AttributesSize = 120;

template <size_t I>
void serializeAttribute(ParticleType &particle, std::vector<char> &attributeVector, size_t &startIndex) {
  const auto attribute = particle.get<Attributes[I]>();
  const auto sizeOfValue = sizeof(attribute);
  std::memcpy(&attributeVector[startIndex], &attribute, sizeOfValue);
  startIndex += sizeOfValue;
}

template <size_t I>
void deserializeAttribute(char *&attributeVector, ParticleType &particle, size_t &startIndex) {
  auto attribute = particle.get<Attributes[I]>();
  const auto sizeOfValue = sizeof(attribute);
  std::memcpy(&attribute, &attributeVector[startIndex], sizeOfValue);
  particle.set<Attributes[I]>(attribute);
  startIndex += sizeOfValue;
}

template <size_t... I>
void serializeParticleImpl(ParticleType &particle, std::vector<char> &serializedParticle, std::index_sequence<I...>) {
  // Calculate size of serialized attributes
  // size_t sizeOfAttributes = 0;
  // auto increaseAttributesSize = [&](size_t addition) { sizeOfAttributes += addition; };
  //(increaseAttributesSize(sizeof(particle.get<Attributes[I]>())), ...);

  // Serialize particle attributes
  size_t startIndex = 0;
  std::vector<char> attributesVector(AttributesSize);
  (serializeAttribute<I>(particle, attributesVector, startIndex), ...);

  // Add serialized attributes to serialized particle
  serializedParticle.insert(serializedParticle.end(), attributesVector.begin(), attributesVector.end());
}

template <size_t... I>
void deserializeParticleImpl(char *particleData, ParticleType &particle, std::index_sequence<I...>) {
  size_t startIndex = 0;
  (deserializeAttribute<I>(particleData, particle, startIndex), ...);
}
}  // namespace

namespace ParticleSerializationTools {
void serializeParticle(ParticleType particle, std::vector<char> &serializedParticles) {
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

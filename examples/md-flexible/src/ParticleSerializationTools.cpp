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

void doubleVectorToCharVector(std::vector<double> ds, std::vector<char> &cs) {
  cs.resize(ds.size() * sizeof(double));
  memcpy(cs.data(), ds.data(), ds.size() * sizeof(double));
}

void charVectorToDoubleVector(std::vector<char> cs, std::vector<double> &ds) {
  ds.resize(cs.size() / sizeof(double));
  memcpy(ds.data(), cs.data(), cs.size());
}

void serializeParticlePositions(const std::vector<ParticleType> &particles, std::vector<char> &serializedParticles) {
  std::vector<double> positions;
  positions.reserve(serializedParticles.size() + particles.size() * 8);
  for (const auto &particle : particles) {
    // save id of particle
    unsigned long id = particle.getID();
    double *did = (double *)(&id);
    positions.push_back(*did);

    // save type id of particle
    unsigned long type = particle.getTypeId();
    double *dtype = (double *)(&type);
    positions.push_back(*dtype);

    // save position of particle
    auto pos = particle.getR();
    positions.push_back(pos.at(0));
    positions.push_back(pos.at(1));
    positions.push_back(pos.at(2));

    // save velocity of particle
    auto v = particle.getV();
    positions.push_back(v.at(0));
    positions.push_back(v.at(1));
    positions.push_back(v.at(2));
  }
  std::vector<char> positionsCharData(reinterpret_cast<char *>(positions.data()),
                                      reinterpret_cast<char *>(positions.data() + positions.size()));
  serializedParticles.insert(serializedParticles.end(), positionsCharData.begin(), positionsCharData.end());
}

void deserializeParticlePositions(std::vector<char> particleData, std::vector<ParticleType> &particles) {
  std::vector<double> particleDoubleData(reinterpret_cast<double *>(particleData.data()),
                                         reinterpret_cast<double *>(particleData.data() + particleData.size()));
  particles.reserve(particles.size() + particleDoubleData.size() / 8);
  for (size_t i = 0; i < particleDoubleData.size(); i += 8) {
    unsigned long *id = (unsigned long *)(&particleDoubleData[i]);
    unsigned long *type = (unsigned long *)(&particleDoubleData[i + 1]);
    ParticleType particle({particleDoubleData[i + 2], particleDoubleData[i + 3], particleDoubleData[i + 4]},
                          {particleDoubleData[i + 5], particleDoubleData[i + 6], particleDoubleData[i + 7]}, *id,
                          *type);
    particles.push_back(particle);
  }
}

void serializeParticleForces(const std::vector<ParticleType> &particles, std::vector<char> &serializedParticles) {
  std::vector<double> forces;
  forces.reserve(serializedParticles.size() + particles.size() * 4);
  for (const auto &particle : particles) {
    // save id of particle
    unsigned long id = particle.getID();
    double *did = (double *)(&id);
    forces.push_back(*did);

    // save position of particle
    auto f = particle.getF();
    forces.push_back(f.at(0));
    forces.push_back(f.at(1));
    forces.push_back(f.at(2));
  }
  std::vector<char> forcesCharData(reinterpret_cast<char *>(forces.data()),
                                   reinterpret_cast<char *>(forces.data() + forces.size()));
  serializedParticles.insert(serializedParticles.end(), forcesCharData.begin(), forcesCharData.end());
}

void deserializeParticleForces(std::vector<char> particleData, std::vector<ParticleType> &particles) {
  std::vector<double> particleDoubleData(reinterpret_cast<double *>(particleData.data()),
                                         reinterpret_cast<double *>(particleData.data() + particleData.size()));
  for (size_t i = 0; i < particleDoubleData.size(); i += 4) {
    unsigned long *id = (unsigned long *)(&particleDoubleData[i]);
    ParticleType particle({0, 0, 0}, {0, 0, 0}, *id);
    particle.setF({particleDoubleData[i + 1], particleDoubleData[i + 2], particleDoubleData[i + 3]});
    particles.push_back(particle);
  }
}

}  // namespace ParticleSerializationTools

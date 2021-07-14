/**
 * @file ParticleSerializationToolsTest.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "ParticleSerializationToolsTest.h"

#include <gtest/gtest.h>

#include "src/ParticleSerializationTools.h"
#include "testingHelpers/commonTypedefs.h"

ParticleSerializationToolsTest::ParticleSerializationToolsTest() : AutoPasTestBase() {
  _particle = {5.0, 6.0,  7.0,  10.0, 20.0, 30.0, 0.05, 0.06, 0.07, 12, autopas::OwnershipState::halo,
               5,   0.01, 0.02, 0.03};

  _molecule.setR(_particle.position);
  _molecule.setV(_particle.velocity);
  _molecule.setF(_particle.force);
  _molecule.setID(_particle.id);
  _molecule.setOwnershipState(_particle.ownershipState);
  _molecule.setTypeId(_particle.typeId);
  _molecule.setOldF(_particle.oldForce);
}

TEST_F(ParticleSerializationToolsTest, testSeralizeParticle) {
  std::vector<char> serializedParticle;
  ParticleSerializationTools::serializeParticle(_molecule, serializedParticle);

  auto *attributes = reinterpret_cast<ParticleAttributes *>(&serializedParticle[0]);

  EXPECT_EQ(attributes->position, _particle.position);
  EXPECT_EQ(attributes->velocity, _particle.velocity);
  EXPECT_EQ(attributes->force, _particle.force);
  EXPECT_EQ(attributes->id, _particle.id);
  EXPECT_EQ(attributes->ownershipState, _particle.ownershipState);
  EXPECT_EQ(attributes->oldForce, _particle.oldForce);
}

TEST_F(ParticleSerializationToolsTest, testDeserializeParticle) {
  std::vector<char> serializedParticle;
  ParticleSerializationTools::serializeParticle(_molecule, serializedParticle);

  autopas::MoleculeLJ<double> deserializedMolecule;
  ParticleSerializationTools::deserializeParticle(&serializedParticle[0], deserializedMolecule);

  EXPECT_EQ(_molecule.getR(), deserializedMolecule.getR());
  EXPECT_EQ(_molecule.getV(), deserializedMolecule.getV());
  EXPECT_EQ(_molecule.getF(), deserializedMolecule.getF());
  EXPECT_EQ(_molecule.getOldF(), deserializedMolecule.getOldF());
  EXPECT_EQ(_molecule.getID(), deserializedMolecule.getID());
  EXPECT_EQ(_molecule.isDummy(), deserializedMolecule.isDummy());
  EXPECT_EQ(_molecule.isOwned(), deserializedMolecule.isOwned());
  EXPECT_EQ(_molecule.isHalo(), deserializedMolecule.isHalo());
  EXPECT_EQ(_molecule.getTypeId(), deserializedMolecule.getTypeId());
}

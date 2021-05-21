/**
 * @file ParticleSerializationToolsTest.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "ParticleSerializationToolsTest.h"

#include "src/ParticleSerializationTools.h"
#include "testingHelpers/commonTypedefs.h"

#include <gtest/gtest.h>

ParticleSerializationToolsTest::ParticleSerializationToolsTest() : AutoPasTestBase() {
	_particle = {
		5.0, 6.0, 7.0,
		10.0, 20.0, 30.0,
		0.05, 0.06, 0.07,
		12,
		autopas::OwnershipState::halo,
		5,
		0.01, 0.02, 0.03
	};

	_molecule.setR({_particle.positionX, _particle.positionY, _particle.positionZ});
	_molecule.setV({_particle.velocityX, _particle.velocityY, _particle.velocityZ});
	_molecule.setF({_particle.forceX, _particle.forceY, _particle.forceZ});
	_molecule.setID(_particle.id);
	_molecule.setOwnershipState(_particle.ownershipState);
	_molecule.setTypeId(_particle.typeId);
	_molecule.setOldF({_particle.oldForceX, _particle.oldForceY, _particle.oldForceZ});
}

TEST_F(ParticleSerializationToolsTest, testSeralizeParticle){
	std::vector<char> serializedParticle;
	ParticleSerializationTools::serializeParticle(_molecule, serializedParticle);
	
	ParticleAttributes* attributes = reinterpret_cast<ParticleAttributes*>(&serializedParticle[0]);

	EXPECT_EQ(attributes->positionX, _particle.positionX);
	EXPECT_EQ(attributes->positionY, _particle.positionY);
	EXPECT_EQ(attributes->positionZ, _particle.positionZ);
	EXPECT_EQ(attributes->velocityX, _particle.positionX);
	EXPECT_EQ(attributes->velocityY, _particle.positionY);
	EXPECT_EQ(attributes->velocityZ, _particle.positionZ);
	EXPECT_EQ(attributes->forceX, _particle.forceX);
	EXPECT_EQ(attributes->forceY, _particle.forceY);
	EXPECT_EQ(attributes->forceZ, _particle.forceZ);
	EXPECT_EQ(attributes->id, _particle.id);
	EXPECT_EQ(attributes->ownershipState, _particle.ownershipState);
	EXPECT_EQ(attributes->oldForceX, _particle.oldForceX);
	EXPECT_EQ(attributes->oldForceY, _particle.oldForceY);
	EXPECT_EQ(attributes->oldForceZ, _particle.oldForceZ);
}

TEST_F(ParticleSerializationToolsTest, testDeserializeParticleData){
	std::vector<char> serializedParticle;
	ParticleSerializationTools::serializeParticle(_molecule, serializedParticle);

	autopas::MoleculeLJ<double> deserializedMolecule;
	ParticleSerializationTools::deserializeParticleData(&serializedParticle[0], deserializedMolecule);

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

TEST_F(ParticleSerializationToolsTest, testConvertParticleAttributesToParticle){
	autopas::MoleculeLJ<double> molecule = ParticleSerializationTools::convertParticleAttributesToParticle(_particle);

	EXPECT_EQ(_molecule.getR(), molecule.getR());	
	EXPECT_EQ(_molecule.getV(), molecule.getV());	
	EXPECT_EQ(_molecule.getF(), molecule.getF());	
	EXPECT_EQ(_molecule.getOldF(), molecule.getOldF());	
	EXPECT_EQ(_molecule.getID(), molecule.getID());	
	EXPECT_EQ(_molecule.isDummy(), molecule.isDummy());	
	EXPECT_EQ(_molecule.isOwned(), molecule.isOwned());	
	EXPECT_EQ(_molecule.isHalo(), molecule.isHalo());	
	EXPECT_EQ(_molecule.getTypeId(), molecule.getTypeId());	
}


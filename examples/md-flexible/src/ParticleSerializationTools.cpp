/**
 * @file ParticleSerializationTools.cpp
 * @author J. Körner
 * @date 13.05.2021
 */
#include "ParticleSerializationTools.h"

namespace ParticleSerializationTools {
	void serializeParticle(ParticleType &particle, std::vector<char> &serializedParticles) {
 		std::vector<char> attributesVector;
		attributesVector.resize(sizeof(ParticleAttributes));

		ParticleAttributes* attributes = reinterpret_cast<ParticleAttributes*>(&attributesVector[0]);

		attributes->positionX = particle.getR()[0];
		attributes->positionY = particle.getR()[1];
		attributes->positionZ = particle.getR()[2];

		attributes->velocityX = particle.getV()[0];
		attributes->velocityY = particle.getV()[1];
		attributes->velocityZ = particle.getV()[2];

		attributes->forceX = particle.getF()[0];
		attributes->forceY = particle.getF()[1];
		attributes->forceZ = particle.getF()[2];

		if (particle.isDummy()){
			attributes->ownershipState = autopas::OwnershipState::dummy;
		}
		else if (particle.isOwned()){
			attributes->ownershipState = autopas::OwnershipState::owned;
		}
		else if (particle.isHalo()){
			attributes->ownershipState = autopas::OwnershipState::halo;
		}

		attributes->id = particle.getID();

		attributes->typeId = particle.getTypeId();
		attributes->oldForceX = particle.getOldF()[0];
		attributes->oldForceY = particle.getOldF()[1];
		attributes->oldForceZ = particle.getOldF()[2];

		serializedParticles.insert(serializedParticles.end(), attributesVector.begin(), attributesVector.end());
	}
	
	void deserializeParticleData(char* particleData, ParticleType &particle){
		ParticleAttributes* attributes = reinterpret_cast<ParticleAttributes*>(particleData);

		particle.setR({attributes->positionX, attributes->positionY, attributes->positionZ});
		particle.setV({attributes->velocityX, attributes->velocityY, attributes->velocityZ});
		particle.setF({attributes->forceX, attributes->forceY, attributes->forceZ});
		particle.setOldF({attributes->oldForceX, attributes->oldForceY, attributes->oldForceZ});
		particle.setID(attributes->id);
		particle.setOwnershipState(attributes->ownershipState);
		particle.setTypeId(attributes->typeId);
	}

	void deserializeParticleData(std::vector<char> &particlesData, std::vector<ParticleType> &particles){
   	ParticleType particle;
		size_t sizeOfParticleAttributes = sizeof(ParticleAttributes);
 		for (int i = 0; i != particlesData.size(); i += sizeOfParticleAttributes){
			deserializeParticleData(&particlesData[i], particle);
   		particles.push_back(particle);
 		}
	}

	ParticleType convertParticleAttributesToParticle(ParticleAttributes &attributes) {
		ParticleType particle;

		particle.setR({attributes.positionX, attributes.positionY, attributes.positionZ});
		particle.setV({attributes.velocityX, attributes.velocityY, attributes.velocityZ});
		particle.setF({attributes.forceX, attributes.forceY, attributes.forceZ});
		particle.setID(attributes.id);
		particle.setTypeId(attributes.typeId);
		particle.setOldF({attributes.oldForceX, attributes.oldForceY, attributes.oldForceZ});

		return particle;
	}
}

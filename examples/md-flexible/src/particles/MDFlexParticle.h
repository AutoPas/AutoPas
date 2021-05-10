/**
 * @file MDFlexParticle.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "autopas/particles/ParticleBase.h"

/**
* A wrapper for autopas' ParticleBase for simple serialization of a particles attributes.
* This class can be used to serialize and deserialize autopas' ParticleBase for use in MPI communication.
* Be careful when using this class. If the layout of ParticleBase.h changes, serialization and deserialization might
* result in wrong data.
* This class needs to be updated after a change to the layout of ParticleBase.h>
*/
class MDFlexParticle : public autopas::ParticleBase<double, unsigned long> {
 public:
  MDFlexParticle();
  ~MDFlexParticle() = default;
	
	/**
	* Struct to simplify particle serialization for MPI messages.
	*/
	struct Attributes { 
		double positionX; 		
		double positionY;
		double positionZ;
		double velocityX; 		
		double velocityY;
		double velocityZ;
		double forceX;
		double forceY;
		double forceZ;
 		unsigned long id; 		
		int64_t ownershipState;
	};

	/**
	* Serializes the attributes of the particle to an array of chars.
	* The size of the array depends on the size of the Attributes struct.
	*/
  char* getSerializedAttributes() const { return reinterpret_cast<char*>(&(_attributes->positionX)); }
	
	/**
	* Deserializes the attributeData into the Attributes struct.
	* This function does not validate the the data. Make sure the data passed with attributeData has the correct layout.
	* @param attributeData the data which will be reinterpreted as Attributes.
	*/
	void setAttributes(char* &attributeData);

	private:
		/**
		* Points to the start of the particles attributes
		*/
		Attributes* _attributes;
};

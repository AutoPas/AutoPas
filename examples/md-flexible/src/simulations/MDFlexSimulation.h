/**
 * @file MDFlexSimulation.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "../Simulation.h"

/**
* Provides base functionality for MD-Flexible simulations. 
* Maybe add unique identifier for each class or an MDFlex factory of some kind
*/
class MDFlexSimulation {
	public:
		MDFlexSimulation();
		~MDFlexSimulation();

		/**
		* Used for parsing and for all other preparations required for the simulation
		*/
		virtual void initialize(int argc, char** argv);

		/**
		* Initializes AutoPas.
		*/
		virtual void initializeAutoPas() = 0;

		/**
		* A wrapper for running the simulation
		*/
		virtual void run();

		/**
		* This may be used to process runtime data / information for further use and / or to
		* do additional required finalization (see MDFlexSingleNode.cpp)
		*/
		virtual void finalize(int argc, char** argv);

	protected:
  		autopas::AutoPas<Simulation::ParticleType>* _autopas;
  		MDFlexConfig* _configuration;
		Simulation* _simulation;

};

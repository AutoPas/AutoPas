/**
 * @file MdFlex.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "Simulation.h"

/**
* A wrapper for MD-Flexible simulations.
* Maybe add unique identifier for each class or an MdFlex factory of some kind
*/
class MdFlex {
	public:
		/**
		* Used for parsing and for all other preparations required for the simulation
		*/
		virtual bool Initialize(int argc, char** argv);

		/**
		* A wrapper for running the simulation
		*/
		virtual bool Run();

		/**
		* This may be used to process runtime data / information for further use and / or to
		* do additional requred finalization (see MdFlexSingleNode.cpp)
		*/
		virtual bool Finalize();

	private:
		Simulation _simulation;
};

/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

class MDFlexMPI : MDFlexSimulation {
	public:
		MDFlexMPI(int argc, char** argv);
		~MDFlexMPI() = default;

		void run() override;

	private:
		int _rank;
		std::array<unsigned int, 3> _processorCoordinates;
};


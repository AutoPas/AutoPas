/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

class MDFlexMPI : MDFlexSimulation {
	public:
		MDFlexMPI() = default;
		~MDFlexMPI() = default;

		void run() override;
};


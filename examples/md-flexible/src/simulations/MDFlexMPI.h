/**
 * @file MDFlexMPI.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

class MDFlexMPI : MDFlexSimulation {
	public:
		void initialize(int argc, char** argv) override;
		void initializeAutoPas() override;
		void run() override;
		void finalize(int argc, char** argv) override;
};

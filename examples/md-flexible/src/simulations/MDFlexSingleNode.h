/**
 * @file MDFlexSingleNode.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

#include <fstream>

class MDFlexSingleNode : MDFlexSimulation {
	public:
		void initialize(int argc, char** argv) override;
		void initializeAutoPas() override;
		void run() override;
		void finalize(int argc, char** argv) override;

	private:
 		std::streambuf* _streamBuffer;
 		std::ofstream _logFile;
 		std::ostream* _outputStream;
};

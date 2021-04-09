/**
 * @file MdFlexSingleNode.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MdFlex.h"

class MdFlexSingleNode : MdFlex {
	public:
		bool Initialize(int argc, char** argv) override;
		bool Run() override;
		bool Finalize() override;
};

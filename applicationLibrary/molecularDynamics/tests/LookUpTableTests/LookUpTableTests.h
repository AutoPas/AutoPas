/**
* @file LookUpTableTest.h
* @author J. Hampe
* @date 4.4.2024
*/
#pragma once
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <vector>

#include "../../../../tests/testAutopas/AutoPasTestBase.h"
#include "molecularDynamicsLibrary/ATLookUpTable.h"
#include "molecularDynamicsLibrary/LJLookUpTable.h"
#include "molecularDynamicsLibrary/LookUpTableTypes.h"
#include "autopas/utils/Timer.h"

#ifndef AUTOPAS_LOOKUPTABLETESTS_H
#define AUTOPAS_LOOKUPTABLETESTS_H

class LookUpTableTest : public AutoPasTestBase {};

#endif  // AUTOPAS_LOOKUPTABLETESTS_H

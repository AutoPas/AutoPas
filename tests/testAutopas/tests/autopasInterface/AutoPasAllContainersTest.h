/**
 * @file AutoPasAllContainersTest.h
 * @author F. Gratl
 * @date 29.03.23
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "testingHelpers/GenerateValidConfigurations.h"

using ParamType = ContainerConfiguration;

class AutoPasAllContainersTest : public AutoPasTestBase, public ::testing::WithParamInterface<ParamType> {};

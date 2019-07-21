/**
 * @file GaussianProcessTest.h
 * @author Jan Nguyen
 * @date 12.06.19
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "Eigen/Dense"
#include "autopas/selectors/FeatureVector.h"
#include "autopas/selectors/tuningStrategy/GaussianProcess.h"
#include "autopas/utils/NumberSet.h"

class GaussianProcessTest : public AutoPasTestBase {};

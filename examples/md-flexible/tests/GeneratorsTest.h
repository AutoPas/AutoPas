/**
 * @file GeneratorsTest.h
 * @author N. Fottner
 * @date 02/08/19
 */
#pragma once
#include <AutoPasTestBase.h>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "Generator.h"
#include "Objects/Objects.h"
#include "PrintableMolecule.h"
#include "TimeDiscretization.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "parsing/YamlParser.h"
#include "testingHelpers/GaussianGenerator.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class GeneratorsTest : public AutoPasTestBase {
 public:
  GeneratorsTest() = default;

 protected:
  double epsilon{1.0};
  double sigma{1.0};
  double cutoff{1.};
  std::array<double, 3> boxmin{{0., 0., 0.}};
  std::array<double, 3> boxmax{{5., 5., 5.}};
};
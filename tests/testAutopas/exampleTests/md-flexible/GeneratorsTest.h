#pragma once
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "../../../../examples/md-flexible/Generator.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../examples/md-flexible/parsing/YamlParser.h"
#include "AutoPasTestBase.h"
#include "Objects/Objects.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "testingHelpers/GaussianGenerator.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"
class GeneratorsTest : public AutoPasTestBase {
 public:
  GeneratorsTest()
      : AutoPasTestBase(), epsilon{1.0}, sigma{1.0}, cutoff{1.}, boxmin{{0., 0., 0.}}, boxmax{{5., 5., 5.}} {}

 protected:
  double epsilon;
  double sigma;
  double cutoff;
  std::array<double, 3> boxmin;
  std::array<double, 3> boxmax;
};
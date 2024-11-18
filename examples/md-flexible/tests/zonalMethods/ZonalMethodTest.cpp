
#include "tests/zonalMethods/ZonalMethodTest.h"

#include "src/zonalMethods/HalfShell.h"
#include "src/zonalMethods/region/RectRegion.h"

/**
 * Define function to allow instantiation
 */
void ZonalMethodTest::collectParticles(AutoPasType &autoPasContainer) {}

/**
 * Define function to allow instantiation
 */
void ZonalMethodTest::SendAndReceiveExports(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm,
                                            std::array<int, 26> allNeighbourIndices) {}
/**
 * Define function to allow instantiation
 */
void ZonalMethodTest::SendAndReceiveResults(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm,
                                            std::array<int, 26> allNeighbourIndices) {}

/**
 * Tests the functionality of the getRectRegionConditional function, by
 * specifying a condition which should result in no regions.
 */
TEST_F(ZonalMethodTest, testRectRegionFalseCalculation) {
  RectRegion homeBoxRegion({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

  std::vector<RectRegion> regions;

  auto condition = [](const int d[3]) { return false; };

  getRectRegionsConditional(homeBoxRegion, 0.5, 0.1, regions, condition);

  EXPECT_EQ(regions.size(), 0);
}

/**
 * Tests the functionality of the getRectRegionConditional function, by
 * specifying a condition which should result in all 26 regions.
 */
TEST_F(ZonalMethodTest, testRectRegionTrueCalculation) {
  RectRegion homeBoxRegion({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

  std::vector<RectRegion> regions;

  auto condition = [](const int d[3]) { return true; };

  getRectRegionsConditional(homeBoxRegion, 0.5, 0.1, regions, condition);

  EXPECT_EQ(regions.size(), 26);
}

/**
 * Tests the functionality of the getRectRegionConditional function, by
 * specifying a condition which should result in the Halfshell import stencil.
 */
TEST_F(ZonalMethodTest, testRectRegionHSImportCalculation) {
  RectRegion homeBoxRegion({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

  std::vector<RectRegion> regions;

  auto condition = [](const int d[3]) {
    /**
     * Stencil:
     *  z > 0 +
     *  z == 0 and y > 0 +
     *  z == 0 and y == 0 and x > 0
     */
    return d[2] > 0 or (d[2] == 0 and (d[1] > 0 or (d[1] == 0 and d[0] > 0)));
  };

  getRectRegionsConditional(homeBoxRegion, 0.5, 0.1, regions, condition);

  EXPECT_EQ(regions.size(), 13);

  std::vector<RectRegion> expectedRegions{};

  // regions with z > 0:
  // 0, 0, 1
  expectedRegions.push_back(RectRegion{{0.0, 0.0, 1.0}, {1.0, 1.0, 0.6}});
  // 0, 1, 1
  expectedRegions.push_back(RectRegion{{0.0, 1.0, 1.0}, {1.0, 0.6, 0.6}});
  // 1, 0, 1
  expectedRegions.push_back(RectRegion{{1.0, 0.0, 1.0}, {0.6, 1.0, 0.6}});
  // 1, 1, 1
  expectedRegions.push_back(RectRegion{{1.0, 1.0, 1.0}, {0.6, 0.6, 0.6}});
  // 1, -1, 1
  expectedRegions.push_back(RectRegion{{1.0, 0.0, 1.0}, {0.6, -0.6, 0.6}});
  // 0, -1, 1
  expectedRegions.push_back(RectRegion{{0.0, 0.0, 1.0}, {1.0, -0.6, 0.6}});
  // -1, -1, 1
  expectedRegions.push_back(RectRegion{{0.0, 0.0, 1.0}, {-0.6, -0.6, 0.6}});
  // -1, 0, 1
  expectedRegions.push_back(RectRegion{{0.0, 0.0, 1.0}, {-0.6, 1.0, 0.6}});
  // -1, 1, 1
  expectedRegions.push_back(RectRegion{{0.0, 1.0, 1.0}, {-0.6, 0.6, 0.6}});

  // regions with z == 0 and y > 0
  // 0, 1, 0
  expectedRegions.push_back(RectRegion{{0.0, 1.0, 0.0}, {1.0, 0.6, 1.0}});
  // -1, 1, 0
  expectedRegions.push_back(RectRegion{{0.0, 1.0, 0.0}, {-0.6, 0.6, 1.0}});
  // 1, 1, 0
  expectedRegions.push_back(RectRegion{{1.0, 1.0, 0.0}, {0.6, 0.6, 1.0}});

  // regions with z == 0 and y == 0 and x > 0
  // 1, 0, 0
  expectedRegions.push_back(RectRegion{{1.0, 0.0, 0.0}, {0.6, 1.0, 1.0}});

  // match exact results
  for (auto r : regions) {
    bool found = false;
    for (size_t i = 0; i < expectedRegions.size(); i++) {
      if (r == expectedRegions[i]) {
        found = true;
        expectedRegions.erase(expectedRegions.begin() + i);
        break;
      }
    }
    EXPECT_TRUE(found);
  }
  EXPECT_EQ(expectedRegions.size(), 0);
}

/**
 * Tests the functionality of the getRectRegionConditional function, by
 * specifying a condition which should result in the Halfshell export stencil.
 */
TEST_F(ZonalMethodTest, testRectRegionHSExportCalculation) {
  RectRegion homeBoxRegion({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

  std::vector<RectRegion> regions;

  auto condition = [](const int d[3]) {
    /**
     * Stencil:
     *  z > 0 +
     *  z == 0 and y > 0 +
     *  z == 0 and y == 0 and x > 0
     */
    return d[2] > 0 or (d[2] == 0 and (d[1] > 0 or (d[1] == 0 and d[0] > 0)));
  };

  getRectRegionsConditional(homeBoxRegion, 0.5, 0.1, regions, condition, false);

  EXPECT_EQ(regions.size(), 13);

  std::vector<RectRegion> expectedRegions{};

  // regions with z > 0:
  // 0, 0, 1
  expectedRegions.push_back(RectRegion{{0.0, 0.0, 1.0}, {1.0, 1.0, -0.6}});
  // 0, 1, 1
  expectedRegions.push_back(RectRegion{{0.0, 1.0, 1.0}, {1.0, -0.6, -0.6}});
  // 1, 0, 1
  expectedRegions.push_back(RectRegion{{1.0, 0.0, 1.0}, {-0.6, 1.0, -0.6}});
  // 1, 1, 1
  expectedRegions.push_back(RectRegion{{1.0, 1.0, 1.0}, {-0.6, -0.6, -0.6}});
  // 1, -1, 1
  expectedRegions.push_back(RectRegion{{1.0, 0.0, 1.0}, {-0.6, 0.6, -0.6}});
  // 0, -1, 1
  expectedRegions.push_back(RectRegion{{0.0, 0.0, 1.0}, {1.0, 0.6, -0.6}});
  // -1, -1, 1
  expectedRegions.push_back(RectRegion{{0.0, 0.0, 1.0}, {0.6, 0.6, -0.6}});
  // -1, 0, 1
  expectedRegions.push_back(RectRegion{{0.0, 0.0, 1.0}, {0.6, 1.0, -0.6}});
  // -1, 1, 1
  expectedRegions.push_back(RectRegion{{0.0, 1.0, 1.0}, {0.6, -0.6, -0.6}});

  // regions with z == 0 and y > 0
  // 0, 1, 0
  expectedRegions.push_back(RectRegion{{0.0, 1.0, 0.0}, {1.0, -0.6, 1.0}});
  // -1, 1, 0
  expectedRegions.push_back(RectRegion{{0.0, 1.0, 0.0}, {0.6, -0.6, 1.0}});
  // 1, 1, 0
  expectedRegions.push_back(RectRegion{{1.0, 1.0, 0.0}, {-0.6, -0.6, 1.0}});

  // regions with z == 0 and y == 0 and x > 0
  // 1, 0, 0
  expectedRegions.push_back(RectRegion{{1.0, 0.0, 0.0}, {-0.6, 1.0, 1.0}});

  // match exact results
  for (auto r : regions) {
    bool found = false;
    for (size_t i = 0; i < expectedRegions.size(); i++) {
      if (r == expectedRegions[i]) {
        found = true;
        expectedRegions.erase(expectedRegions.begin() + i);
        break;
      }
    }
    EXPECT_TRUE(found);
  }
  EXPECT_EQ(expectedRegions.size(), 0);
}

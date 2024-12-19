
#include "tests/zonalMethods/MidpointTest.h"

TEST_F(MidpointTest, testMidpointInitialization) {
  ASSERT_EQ(_exportRegions.size(), 26);
  ASSERT_EQ(_importRegions.size(), 26);
  ASSERT_EQ(_interactionZones.size(), 26);

  auto sum = 0;
  for (auto schedule : _interactionSchedule) {
    sum += schedule.second.size();
  }

  /**
   * This is a weak check, but it suffices for now.
   * (see Spahl 2016)
   * Number of total interactions =
   *  13 interactions with the neighbour on the other side +
   *  6 times 8 interaction (see Figure 3.5) +
   *  12 times 2 interaction (see Figure 3.6) =
   *  13 + 48 + 24 = 85
   */
  ASSERT_EQ(sum, 85);

  // some individual checks
  auto schedule = _interactionSchedule.at(std::to_string(convRelNeighboursToIndex({-1, 0, 0})));
  ASSERT_EQ(schedule.size(), 9);
  std::vector<std::string> expected;
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      expected.push_back(std::to_string(convRelNeighboursToIndex({1, i, j})));
    }
  }
  std::sort(schedule.begin(), schedule.end());
  std::sort(expected.begin(), expected.end());
  ASSERT_EQ(schedule, expected);

  schedule = _interactionSchedule.at(std::to_string(convRelNeighboursToIndex({1, 1, 0})));
  ASSERT_EQ(schedule.size(), 2);
  expected.clear();
  expected.push_back(std::to_string(convRelNeighboursToIndex({-1, -1, 1})));
  expected.push_back(std::to_string(convRelNeighboursToIndex({-1, -1, -1})));

  schedule = _interactionSchedule.at(std::to_string(convRelNeighboursToIndex({-1, -1, -1})));
  ASSERT_EQ(schedule.size(), 1);

  schedule = _interactionSchedule.at(std::to_string(convRelNeighboursToIndex({1, 1, 1})));
  ASSERT_EQ(schedule.size(), 0);
}

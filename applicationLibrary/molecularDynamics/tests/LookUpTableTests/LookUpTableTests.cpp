/**
* @file LookUpTableTest.cpp
* @author J. Hampe
* @date 4.4.2024
*/

#include "LookUpTableTests.h"
#include "testingHelpers/commonTypedefs.h"
#include "autopas/utils/ArrayMath.h"


/**
 * Very simple test to check that rotate still works
 */
TEST_F(LookUpTableTest, MostBasicRotate) {
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  double cutoff = 2.5;
  double nu = 0.073;
  std::vector<autopas::utils::Timer> LUTtimers = {autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer()};
  ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                  ForceLookUpTable::nextNeighbor, double, int> ATLookUpTable = std::move(ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                                            ForceLookUpTable::nextNeighbor, double, int>(
      {cutoff * cutoff, nu, 8.0}, &LUTtimers));
  ATLookUpTable.retrieveValue({0., 0., 0.}, {std::sqrt(0.75), 0., 0.5}, {0., 0., 1.}, 1., 1., 1.); // Change so targetB ends up on 1,0,0 to check second rotation
}

TEST_F(LookUpTableTest, SehrFalsch) {
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  double cutoff = 2.5;
  double nu = 0.073;
  std::vector<autopas::utils::Timer> LUTtimers = {autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer()};
  ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                  ForceLookUpTable::nextNeighbor, double, int> ATLookUpTable = std::move(ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                                                ForceLookUpTable::nextNeighbor, double, int>(
          {cutoff * cutoff, nu, 100.0}, &LUTtimers));
  ATLookUpTable.retrieveValue({2.703957905911746, -0.0026572954522118756, 6.610864923489291}, {3.3706888660320056, 0.3937881695307143, 7.713014487295848}, {2.031792390399456, -1.60843618566604, 7.721119056770151}, 1.816432840887397, 5.801611825050684, 4.262996564967491); // Change so targetB ends up on 1,0,0 to check second rotation
}

/**
 * Find the rotate error from first rotation
 */
TEST_F(LookUpTableTest, RotateErrorFirstRot) {
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::critical);
  double cutoff = 2.5;
  double nu = 0.073;
  std::vector<autopas::utils::Timer> LUTtimers = {autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer()};
  ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                  ForceLookUpTable::nextNeighbor, double, int> ATLookUpTable = std::move(ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                                                ForceLookUpTable::nextNeighbor, double, int>(
          {cutoff * cutoff, nu, 8.0}, &LUTtimers));
  //ATLookUpTable.retrieveValue({0., 0., 0.}, {std::sqrt(0.75), 0., 0.5}, {0., 0., 1.}, 1., 1., 1.); // Change so targetB ends up on 1,0,0 to check second rotation
  ATLookUpTable.retrieveValue({-0.0038729557731159102, 0.00040231088697253504, 0.00012919291748858297}, {1.346814631105968, 0.001708787507448397, -0.0013156195867966409}, {-0.0006312167672606656, 1.558781017264059, 1.1007112931787657}, 1.8243607517135754, 5.454547558089624, 3.6398356607768276);

}

/**
 * Find the rotate error from second rotation
 */
TEST_F(LookUpTableTest, RotateErrorSecondRot) {
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::critical);
  double cutoff = 2.5;
  double nu = 0.073;
  std::vector<autopas::utils::Timer> LUTtimers = {autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer()};
  ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                  ForceLookUpTable::nextNeighbor, double, int> ATLookUpTable = std::move(ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                                                ForceLookUpTable::nextNeighbor, double, int>(
          {cutoff * cutoff, nu, 8.0}, &LUTtimers));
  ATLookUpTable.retrieveValue({0., 0., 0.}, {std::sqrt(0.75), 0., 0.5}, {0., 0., 1.}, 1., 1., 1.); // Change so targetB ends up on 1,0,0 to check second rotation
}

/**
 * Find the rotate error increasing distances
 */
TEST_F(LookUpTableTest, RotateErrorIncreasingDistances) {
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::critical);
  double cutoff = 2.5;
  double nu = 0.073;
  std::vector<autopas::utils::Timer> LUTtimers = {autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer()};
  ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                  ForceLookUpTable::nextNeighbor, double, int> ATLookUpTable = std::move(ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                                                ForceLookUpTable::nextNeighbor, double, int>(
          {cutoff * cutoff, nu, 8.0}, &LUTtimers));
  ATLookUpTable.retrieveValue({0., 0., 0.}, {std::sqrt(0.75), 0., 0.5}, {0., 0., 1.}, 1., 1., 1.); // Change so targetB ends up on 1,0,0 to check second rotation
}


/**
 * Find the rotate error over internal triangle angles
 */
TEST_F(LookUpTableTest, RotateErrorIncreasingAngle) {
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  double cutoff = 2.5;
  double nu = 0.073;
  std::vector<autopas::utils::Timer> LUTtimers = {autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer()};
  ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                  ForceLookUpTable::nextNeighbor, double, int> ATLookUpTable = std::move(ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                                                ForceLookUpTable::nextNeighbor, double, int>(
          {cutoff * cutoff, nu, 8.0}, &LUTtimers));
  ATLookUpTable.retrieveValue({0., 0., 0.}, {std::sqrt(0.75), 0., 0.5}, {0., 0., 1.}, 1., 1., 1.); // Change so targetB ends up on 1,0,0 to check second rotation
}


/**
 * Find the rotate error over LUT size
 */
TEST_F(LookUpTableTest, RotateErrorIncreasingLUTsize) {
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  double cutoff = 2.5;
  double nu = 0.073;
  std::vector<autopas::utils::Timer> LUTtimers = {autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer()};
  ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                  ForceLookUpTable::nextNeighbor, double, int> ATLookUpTable = std::move(ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                                                ForceLookUpTable::nextNeighbor, double, int>(
          {cutoff * cutoff, nu, 8.0}, &LUTtimers));
  ATLookUpTable.retrieveValue({0., 0., 0.}, {std::sqrt(0.75), 0., 0.5}, {0., 0., 1.}, 1., 1., 1.); // Change so targetB ends up on 1,0,0 to check second rotation
}

/**
 * Find the rotate error over triangle rotation angles
 */
TEST_F(LookUpTableTest, RotateErrorTriangleRotation) {
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  double cutoff = 2.5;
  double nu = 0.073;
  std::vector<autopas::utils::Timer> LUTtimers = {autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer(), autopas::utils::Timer()};
  ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                  ForceLookUpTable::nextNeighbor, double, int> ATLookUpTable = std::move(ForceLookUpTable::ATLookUpTable<ForceLookUpTable::relative, ForceLookUpTable::evenSpacing,
                                                                ForceLookUpTable::nextNeighbor, double, int>(
          {cutoff * cutoff, nu, 8.0}, &LUTtimers));
  ATLookUpTable.retrieveValue({0., 0., 0.}, {std::sqrt(0.75), 0., 0.5}, {0., 0., 1.}, 1., 1., 1.); // Change so targetB ends up on 1,0,0 to check second rotation
}
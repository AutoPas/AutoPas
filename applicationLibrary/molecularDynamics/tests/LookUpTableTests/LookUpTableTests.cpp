/**
* @file LookUpTableTest.cpp
* @author J. Hampe
* @date 4.4.2024
*/

#include "LookUpTableTests.h"
#include "testingHelpers/commonTypedefs.h"


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

/**
 * Find the rotate error increasing distances
 */
TEST_F(LookUpTableTest, RotateErrorIncreasingDistances) {
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
 * Find the rotate error over angles
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
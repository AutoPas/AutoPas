/**
 * @file SoATest.cpp
 * @author seckler
 * @date 08.06.18
 */

#include "SoATest.h"

namespace SoATest {

TEST_F(SoATest, testInitialization) { autopas::SoA<autopas::Particle::SoAArraysType> soa; }

TEST_F(SoATest, SoATypeTest) {
  static_assert(std::is_same<autopas::utils::SoAType<size_t, double, double, double>::Type,
                             std::tuple<std::vector<size_t, autopas::AlignedAllocator<size_t>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>>>::value,
                "type must be same");
}

using soatype = autopas::utils::SoAType<size_t, double, double, double>::Type;

TEST_F(SoATest, SoAStorageTestGet) {
  autopas::utils::SoAStorage<soatype> soAStorage;

  EXPECT_EQ(soAStorage.get<0>().size(), 0);
  EXPECT_EQ(soAStorage.get<1>().size(), 0);
  EXPECT_EQ(soAStorage.get<2>().size(), 0);
  EXPECT_EQ(soAStorage.get<3>().size(), 0);

  soAStorage.get<0>().resize(4);
  soAStorage.get<1>().resize(5);
  soAStorage.get<2>().resize(6);
  soAStorage.get<3>().resize(7);

  EXPECT_EQ(soAStorage.get<0>().size(), 4);
  EXPECT_EQ(soAStorage.get<1>().size(), 5);
  EXPECT_EQ(soAStorage.get<2>().size(), 6);
  EXPECT_EQ(soAStorage.get<3>().size(), 7);

  static_assert(std::is_same<decltype(soAStorage.get<0>().data()), size_t *>::value, "id type must be proper(size_t)");
  static_assert(std::is_same<decltype(soAStorage.get<1>().data()), double *>::value,
                "position type must be proper(double)");
  static_assert(std::is_same<decltype(soAStorage.get<2>().data()), double *>::value,
                "position type must be proper(double)");
  static_assert(std::is_same<decltype(soAStorage.get<3>().data()), double *>::value,
                "position type must be proper(double)");
}

TEST_F(SoATest, SoAStorageTestAlignment) {
  autopas::utils::SoAStorage<soatype> soAStorage;

  // check alignment to autopas::DEFAULT_CACHE_LINE_SIZE
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soAStorage.get<0>().data()) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soAStorage.get<1>().data()) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soAStorage.get<2>().data()) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soAStorage.get<3>().data()) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);

  soAStorage.get<0>().resize(4);
  soAStorage.get<1>().resize(5);
  soAStorage.get<2>().resize(6);
  soAStorage.get<3>().resize(7);

  // check alignment to autopas::DEFAULT_CACHE_LINE_SIZE
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soAStorage.get<0>().data()) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soAStorage.get<1>().data()) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soAStorage.get<2>().data()) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soAStorage.get<3>().data()) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);
}

TEST_F(SoATest, SoAStorageTestApply) {
  autopas::utils::SoAStorage<soatype> soAStorage;

  EXPECT_EQ(soAStorage.get<0>().size(), 0);
  EXPECT_EQ(soAStorage.get<1>().size(), 0);
  EXPECT_EQ(soAStorage.get<2>().size(), 0);
  EXPECT_EQ(soAStorage.get<3>().size(), 0);

  soAStorage.apply([](auto &list) { list.resize(2); });

  EXPECT_EQ(soAStorage.get<0>().size(), 2);
  EXPECT_EQ(soAStorage.get<1>().size(), 2);
  EXPECT_EQ(soAStorage.get<2>().size(), 2);
  EXPECT_EQ(soAStorage.get<3>().size(), 2);
}

TEST_F(SoATest, SoATestPush) {
  // default soa using autopas::Particle
  using autopas::Particle;
  autopas::SoA<Particle::SoAArraysType> soa;

  EXPECT_EQ(soa.getNumParticles(), 0);

  soa.push<Particle::AttributeNames::id>(2);
  soa.push<Particle::AttributeNames::posX>(0.3);
  soa.push<Particle::AttributeNames::posY>(0.1);
  soa.push<Particle::AttributeNames::posZ>(0.5);
  soa.push<Particle::AttributeNames::forceX>(-0.2);
  soa.push<Particle::AttributeNames::forceY>(0.7);
  soa.push<Particle::AttributeNames::forceZ>(0.07);

  EXPECT_EQ(soa.getNumParticles(), 1);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 0.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 0.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 0.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), -0.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 0.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 0.07);
}

TEST_F(SoATest, SoATestAppend) {
  // default soa using autopas::Particle
  using autopas::Particle;
  std::array<autopas::SoA<Particle::SoAArraysType>, 2> soaBuffer;

  soaBuffer[0].push<Particle::AttributeNames::id>(2);
  soaBuffer[0].push<Particle::AttributeNames::posX>(0.3);
  soaBuffer[0].push<Particle::AttributeNames::posY>(0.1);
  soaBuffer[0].push<Particle::AttributeNames::posZ>(0.5);
  soaBuffer[0].push<Particle::AttributeNames::forceX>(-0.2);
  soaBuffer[0].push<Particle::AttributeNames::forceY>(0.7);
  soaBuffer[0].push<Particle::AttributeNames::forceZ>(0.07);

  EXPECT_EQ(soaBuffer[0].getNumParticles(), 1);
  EXPECT_EQ(soaBuffer[1].getNumParticles(), 0);
  // Append to empty buffer
  soaBuffer[1].append(soaBuffer[0]);
  EXPECT_EQ(soaBuffer[0].getNumParticles(), 1);
  EXPECT_EQ(soaBuffer[1].getNumParticles(), 1);

  EXPECT_EQ(soaBuffer[1].read<Particle::AttributeNames::id>(0), 2);
  EXPECT_EQ(soaBuffer[1].read<Particle::AttributeNames::posX>(0), 0.3);
  EXPECT_EQ(soaBuffer[1].read<Particle::AttributeNames::posY>(0), 0.1);
  EXPECT_EQ(soaBuffer[1].read<Particle::AttributeNames::posZ>(0), 0.5);
  EXPECT_EQ(soaBuffer[1].read<Particle::AttributeNames::forceX>(0), -0.2);
  EXPECT_EQ(soaBuffer[1].read<Particle::AttributeNames::forceY>(0), 0.7);
  EXPECT_EQ(soaBuffer[1].read<Particle::AttributeNames::forceZ>(0), 0.07);
  // Append to filled buffer
  soaBuffer[0].append(soaBuffer[1]);
  EXPECT_EQ(soaBuffer[0].getNumParticles(), 2);
  EXPECT_EQ(soaBuffer[1].getNumParticles(), 1);
}

TEST_F(SoATest, SoATestSwap) {
  // default soa using autopas::Particle
  using autopas::Particle;
  autopas::SoA<Particle::SoAArraysType> soa;

  soa.resizeArrays(2);

  soa.begin<Particle::AttributeNames::id>()[0] = 3;
  soa.begin<Particle::AttributeNames::posX>()[0] = 1.3;
  soa.begin<Particle::AttributeNames::posY>()[0] = 1.1;
  soa.begin<Particle::AttributeNames::posZ>()[0] = 1.5;
  soa.begin<Particle::AttributeNames::forceX>()[0] = -1.2;
  soa.begin<Particle::AttributeNames::forceY>()[0] = 1.7;
  soa.begin<Particle::AttributeNames::forceZ>()[0] = 1.07;

  soa.begin<Particle::AttributeNames::id>()[1] = 1;
  soa.begin<Particle::AttributeNames::posX>()[1] = 10.3;
  soa.begin<Particle::AttributeNames::posY>()[1] = 10.1;
  soa.begin<Particle::AttributeNames::posZ>()[1] = 10.5;
  soa.begin<Particle::AttributeNames::forceX>()[1] = -10.2;
  soa.begin<Particle::AttributeNames::forceY>()[1] = 10.7;
  soa.begin<Particle::AttributeNames::forceZ>()[1] = 10.07;

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 1.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 1.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 1.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), -1.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 1.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 1.07);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(1), 1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(1), 10.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(1), 10.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(1), 10.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(1), -10.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(1), 10.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(1), 10.07);

  soa.swap(0, 1);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(1), 3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(1), 1.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(1), 1.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(1), 1.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(1), -1.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(1), 1.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(1), 1.07);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 10.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 10.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 10.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), -10.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 10.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 10.07);
}

TEST_F(SoATest, SoATestMultiWriteRead) {
  // default soa using autopas::Particle
  using autopas::Particle;
  autopas::SoA<Particle::SoAArraysType> soa;

  soa.resizeArrays(1);

  soa.begin<Particle::AttributeNames::id>()[0] = 1;

  EXPECT_EQ(soa.getNumParticles(), 1);

  soa.writeMultiple<Particle::AttributeNames::posX, Particle::AttributeNames::posY, Particle::AttributeNames::posZ>(
      0, {4., 5., 6.});

  std::array<double, 3> f = {7., 8., 9.};
  soa.writeMultiple<Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
                    Particle::AttributeNames::forceZ>(0, f);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 4.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 5.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 6.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), 7.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 8.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 9.);

  auto res = soa.readMultiple<Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
                              Particle::AttributeNames::forceZ>(0);

  static_assert(std::is_same<decltype(res), std::array<double, 3>>::value, "id type must be proper(size_t)");

  EXPECT_EQ(res[0], 7.);
  EXPECT_EQ(res[1], 8.);
  EXPECT_EQ(res[2], 9.);
}

// this test makes certain that methods don't destroy the inner state of the test
TEST_F(SoATest, SoATestComplicatedAccess) {
  // default soa using autopas::Particle
  using autopas::Particle;
  autopas::SoA<Particle::SoAArraysType> soa;

  EXPECT_EQ(soa.getNumParticles(), 0);

  soa.push<Particle::AttributeNames::id>(2);
  soa.push<Particle::AttributeNames::posX>(0.3);
  soa.push<Particle::AttributeNames::posY>(0.1);
  soa.push<Particle::AttributeNames::posZ>(0.5);
  soa.push<Particle::AttributeNames::forceX>(-0.2);
  soa.push<Particle::AttributeNames::forceY>(0.7);
  soa.push<Particle::AttributeNames::forceZ>(0.07);

  EXPECT_EQ(soa.getNumParticles(), 1);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 0.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 0.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 0.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), -0.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 0.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 0.07);

  soa.clear();

  EXPECT_EQ(soa.getNumParticles(), 0);

  soa.push<Particle::AttributeNames::id>(3);
  soa.push<Particle::AttributeNames::posX>(1.3);
  soa.push<Particle::AttributeNames::posY>(1.1);
  soa.push<Particle::AttributeNames::posZ>(1.5);
  soa.push<Particle::AttributeNames::forceX>(-1.2);
  soa.push<Particle::AttributeNames::forceY>(1.7);
  soa.push<Particle::AttributeNames::forceZ>(1.07);

  EXPECT_EQ(soa.getNumParticles(), 1);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 1.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 1.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 1.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), -1.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 1.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 1.07);

  soa.resizeArrays(2);

  soa.begin<Particle::AttributeNames::id>()[1] = 1;
  soa.begin<Particle::AttributeNames::posX>()[1] = 10.3;
  soa.begin<Particle::AttributeNames::posY>()[1] = 10.1;
  soa.begin<Particle::AttributeNames::posZ>()[1] = 10.5;
  soa.begin<Particle::AttributeNames::forceX>()[1] = -10.2;
  soa.begin<Particle::AttributeNames::forceY>()[1] = 10.7;
  soa.begin<Particle::AttributeNames::forceZ>()[1] = 10.07;

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 1.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 1.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 1.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), -1.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 1.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 1.07);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(1), 1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(1), 10.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(1), 10.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(1), 10.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(1), -10.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(1), 10.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(1), 10.07);

  soa.swap(0, 1);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(1), 3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(1), 1.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(1), 1.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(1), 1.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(1), -1.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(1), 1.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(1), 1.07);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 10.3);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 10.1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 10.5);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), -10.2);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 10.7);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 10.07);

  soa.pop_back();

  EXPECT_EQ(soa.getNumParticles(), 1);

  soa.writeMultiple<Particle::AttributeNames::posX, Particle::AttributeNames::posY, Particle::AttributeNames::posZ>(
      0, {4., 5., 6.});

  std::array<double, 3> f = {7., 8., 9.};
  soa.writeMultiple<Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
                    Particle::AttributeNames::forceZ>(0, f);

  EXPECT_EQ(soa.read<Particle::AttributeNames::id>(0), 1);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posX>(0), 4.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posY>(0), 5.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::posZ>(0), 6.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceX>(0), 7.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceY>(0), 8.);
  EXPECT_EQ(soa.read<Particle::AttributeNames::forceZ>(0), 9.);

  auto res = soa.readMultiple<Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
                              Particle::AttributeNames::forceZ>(0);

  static_assert(std::is_same<decltype(res), std::array<double, 3>>::value, "id type must be proper(size_t)");

  EXPECT_EQ(res[0], 7.);
  EXPECT_EQ(res[1], 8.);
  EXPECT_EQ(res[2], 9.);
}
} // end namespace SoATest

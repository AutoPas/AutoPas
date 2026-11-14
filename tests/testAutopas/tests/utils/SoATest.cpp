/**
 * @file SoATest.cpp
 * @author seckler
 * @date 08.06.18
 */

#include "SoATest.h"

#include <span>
#include <utility>

/**
 * Tests default construction of SoA<ParticleBaseFP64::SoAArraysType>.
 */
TEST_F(SoATest, testInitialization) { autopas::SoA<autopas::ParticleBaseFP64::SoAArraysType> soa; }

/**
 * Checks at compile time that utils::SoAType produces a tuple of AlignedAllocator-backed std::vectors.
 */
TEST_F(SoATest, SoATypeTest) {
  static_assert(std::is_same<autopas::utils::SoAType<size_t, double, double, double>::Type,
                             std::tuple<std::vector<size_t, autopas::AlignedAllocator<size_t>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>,
                                        std::vector<double, autopas::AlignedAllocator<double>>>>::value,
                "type must be same");
}

using soatype = autopas::utils::SoAType<size_t, double, double, double>::Type;

/**
 * Tests SoAStorage::get() default empty sizes, independent per-attribute resizing, and pointer element types.
 */
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

/**
 * Verifies that SoAStorage base vectors remain aligned to autopas::DEFAULT_CACHE_LINE_SIZE before and after resizing.
 */
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

/**
 * Tests SoAStorage::apply() applying an operation uniformly across all attribute vectors in the storage tuple.
 */
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

/**
 * Tests SoA::push() populating base attributes for a particle and verifies size() and read().
 */
TEST_F(SoATest, SoATestPush) {
  // default soa using autopas::ParticleBaseFP64
  using autopas::ParticleBaseFP64;
  autopas::SoA<ParticleBaseFP64::SoAArraysType> soa;

  EXPECT_EQ(soa.size(), 0);

  soa.push<ParticleBaseFP64::AttributeNames::id>(2);
  soa.push<ParticleBaseFP64::AttributeNames::posX>(0.3);
  soa.push<ParticleBaseFP64::AttributeNames::posY>(0.1);
  soa.push<ParticleBaseFP64::AttributeNames::posZ>(0.5);
  soa.push<ParticleBaseFP64::AttributeNames::forceX>(-0.2);
  soa.push<ParticleBaseFP64::AttributeNames::forceY>(0.7);
  soa.push<ParticleBaseFP64::AttributeNames::forceZ>(0.07);

  EXPECT_EQ(soa.size(), 1);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 0.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 0.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 0.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), -0.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 0.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 0.07);
}

/**
 * Tests SoA::append(const SoA &) into both empty and already-populated destination buffers.
 */
TEST_F(SoATest, SoATestAppend) {
  // default soa using autopas::ParticleBaseFP64
  using autopas::ParticleBaseFP64;
  std::array<autopas::SoA<ParticleBaseFP64::SoAArraysType>, 2> soaBuffer;

  soaBuffer[0].push<ParticleBaseFP64::AttributeNames::id>(2);
  soaBuffer[0].push<ParticleBaseFP64::AttributeNames::posX>(0.3);
  soaBuffer[0].push<ParticleBaseFP64::AttributeNames::posY>(0.1);
  soaBuffer[0].push<ParticleBaseFP64::AttributeNames::posZ>(0.5);
  soaBuffer[0].push<ParticleBaseFP64::AttributeNames::forceX>(-0.2);
  soaBuffer[0].push<ParticleBaseFP64::AttributeNames::forceY>(0.7);
  soaBuffer[0].push<ParticleBaseFP64::AttributeNames::forceZ>(0.07);

  EXPECT_EQ(soaBuffer[0].size(), 1);
  EXPECT_EQ(soaBuffer[1].size(), 0);
  // Append to empty buffer
  soaBuffer[1].append(soaBuffer[0]);
  EXPECT_EQ(soaBuffer[0].size(), 1);
  EXPECT_EQ(soaBuffer[1].size(), 1);

  EXPECT_EQ(soaBuffer[1].read<ParticleBaseFP64::AttributeNames::id>(0), 2);
  EXPECT_EQ(soaBuffer[1].read<ParticleBaseFP64::AttributeNames::posX>(0), 0.3);
  EXPECT_EQ(soaBuffer[1].read<ParticleBaseFP64::AttributeNames::posY>(0), 0.1);
  EXPECT_EQ(soaBuffer[1].read<ParticleBaseFP64::AttributeNames::posZ>(0), 0.5);
  EXPECT_EQ(soaBuffer[1].read<ParticleBaseFP64::AttributeNames::forceX>(0), -0.2);
  EXPECT_EQ(soaBuffer[1].read<ParticleBaseFP64::AttributeNames::forceY>(0), 0.7);
  EXPECT_EQ(soaBuffer[1].read<ParticleBaseFP64::AttributeNames::forceZ>(0), 0.07);
  // Append to filled buffer
  soaBuffer[0].append(soaBuffer[1]);
  EXPECT_EQ(soaBuffer[0].size(), 2);
  EXPECT_EQ(soaBuffer[1].size(), 1);
}

/**
 * Tests SoA::swap() exchanging all attribute values between two particle indices.
 */
TEST_F(SoATest, SoATestSwap) {
  // default soa using autopas::ParticleBaseFP64
  using autopas::ParticleBaseFP64;
  autopas::SoA<ParticleBaseFP64::SoAArraysType> soa;

  soa.resizeArrays(2);

  soa.begin<ParticleBaseFP64::AttributeNames::id>()[0] = 3;
  soa.begin<ParticleBaseFP64::AttributeNames::posX>()[0] = 1.3;
  soa.begin<ParticleBaseFP64::AttributeNames::posY>()[0] = 1.1;
  soa.begin<ParticleBaseFP64::AttributeNames::posZ>()[0] = 1.5;
  soa.begin<ParticleBaseFP64::AttributeNames::forceX>()[0] = -1.2;
  soa.begin<ParticleBaseFP64::AttributeNames::forceY>()[0] = 1.7;
  soa.begin<ParticleBaseFP64::AttributeNames::forceZ>()[0] = 1.07;

  soa.begin<ParticleBaseFP64::AttributeNames::id>()[1] = 1;
  soa.begin<ParticleBaseFP64::AttributeNames::posX>()[1] = 10.3;
  soa.begin<ParticleBaseFP64::AttributeNames::posY>()[1] = 10.1;
  soa.begin<ParticleBaseFP64::AttributeNames::posZ>()[1] = 10.5;
  soa.begin<ParticleBaseFP64::AttributeNames::forceX>()[1] = -10.2;
  soa.begin<ParticleBaseFP64::AttributeNames::forceY>()[1] = 10.7;
  soa.begin<ParticleBaseFP64::AttributeNames::forceZ>()[1] = 10.07;

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 1.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 1.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 1.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), -1.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 1.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 1.07);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(1), 1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(1), 10.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(1), 10.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(1), 10.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(1), -10.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(1), 10.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(1), 10.07);

  soa.swap(0, 1);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(1), 3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(1), 1.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(1), 1.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(1), 1.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(1), -1.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(1), 1.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(1), 1.07);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 10.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 10.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 10.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), -10.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 10.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 10.07);
}

/**
 * Tests SoA::writeMultiple() (with initializer-list and std::array inputs) and SoA::readMultiple() on base attributes.
 */
TEST_F(SoATest, SoATestMultiWriteRead) {
  // default soa using autopas::ParticleBaseFP64
  using autopas::ParticleBaseFP64;
  autopas::SoA<ParticleBaseFP64::SoAArraysType> soa;

  soa.resizeArrays(1);

  soa.begin<ParticleBaseFP64::AttributeNames::id>()[0] = 1;

  EXPECT_EQ(soa.size(), 1);

  soa.writeMultiple<ParticleBaseFP64::AttributeNames::posX, ParticleBaseFP64::AttributeNames::posY,
                    ParticleBaseFP64::AttributeNames::posZ>(0, {4., 5., 6.});

  std::array<double, 3> f = {7., 8., 9.};
  soa.writeMultiple<ParticleBaseFP64::AttributeNames::forceX, ParticleBaseFP64::AttributeNames::forceY,
                    ParticleBaseFP64::AttributeNames::forceZ>(0, f);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 4.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 5.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 6.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), 7.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 8.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 9.);

  auto res = soa.readMultiple<ParticleBaseFP64::AttributeNames::forceX, ParticleBaseFP64::AttributeNames::forceY,
                              ParticleBaseFP64::AttributeNames::forceZ>(0);

  static_assert(std::is_same<decltype(res), std::array<double, 3>>::value, "id type must be proper(size_t)");

  EXPECT_EQ(res[0], 7.);
  EXPECT_EQ(res[1], 8.);
  EXPECT_EQ(res[2], 9.);
}

/**
 * Exercises a mixed sequence of push, clear, resizeArrays, begin, swap, pop_back, writeMultiple, and readMultiple
 * to ensure operations do not corrupt internal SoA state.
 */
TEST_F(SoATest, SoATestComplicatedAccess) {
  // default soa using autopas::ParticleBaseFP64
  using autopas::ParticleBaseFP64;
  autopas::SoA<ParticleBaseFP64::SoAArraysType> soa;

  EXPECT_EQ(soa.size(), 0);

  soa.push<ParticleBaseFP64::AttributeNames::id>(2);
  soa.push<ParticleBaseFP64::AttributeNames::posX>(0.3);
  soa.push<ParticleBaseFP64::AttributeNames::posY>(0.1);
  soa.push<ParticleBaseFP64::AttributeNames::posZ>(0.5);
  soa.push<ParticleBaseFP64::AttributeNames::forceX>(-0.2);
  soa.push<ParticleBaseFP64::AttributeNames::forceY>(0.7);
  soa.push<ParticleBaseFP64::AttributeNames::forceZ>(0.07);

  EXPECT_EQ(soa.size(), 1);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 0.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 0.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 0.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), -0.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 0.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 0.07);

  soa.clear();

  EXPECT_EQ(soa.size(), 0);

  soa.push<ParticleBaseFP64::AttributeNames::id>(3);
  soa.push<ParticleBaseFP64::AttributeNames::posX>(1.3);
  soa.push<ParticleBaseFP64::AttributeNames::posY>(1.1);
  soa.push<ParticleBaseFP64::AttributeNames::posZ>(1.5);
  soa.push<ParticleBaseFP64::AttributeNames::forceX>(-1.2);
  soa.push<ParticleBaseFP64::AttributeNames::forceY>(1.7);
  soa.push<ParticleBaseFP64::AttributeNames::forceZ>(1.07);

  EXPECT_EQ(soa.size(), 1);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 1.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 1.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 1.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), -1.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 1.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 1.07);

  soa.resizeArrays(2);

  soa.begin<ParticleBaseFP64::AttributeNames::id>()[1] = 1;
  soa.begin<ParticleBaseFP64::AttributeNames::posX>()[1] = 10.3;
  soa.begin<ParticleBaseFP64::AttributeNames::posY>()[1] = 10.1;
  soa.begin<ParticleBaseFP64::AttributeNames::posZ>()[1] = 10.5;
  soa.begin<ParticleBaseFP64::AttributeNames::forceX>()[1] = -10.2;
  soa.begin<ParticleBaseFP64::AttributeNames::forceY>()[1] = 10.7;
  soa.begin<ParticleBaseFP64::AttributeNames::forceZ>()[1] = 10.07;

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 1.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 1.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 1.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), -1.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 1.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 1.07);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(1), 1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(1), 10.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(1), 10.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(1), 10.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(1), -10.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(1), 10.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(1), 10.07);

  soa.swap(0, 1);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(1), 3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(1), 1.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(1), 1.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(1), 1.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(1), -1.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(1), 1.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(1), 1.07);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 10.3);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 10.1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 10.5);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), -10.2);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 10.7);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 10.07);

  soa.pop_back();

  EXPECT_EQ(soa.size(), 1);

  soa.writeMultiple<ParticleBaseFP64::AttributeNames::posX, ParticleBaseFP64::AttributeNames::posY,
                    ParticleBaseFP64::AttributeNames::posZ>(0, {4., 5., 6.});

  std::array<double, 3> f = {7., 8., 9.};
  soa.writeMultiple<ParticleBaseFP64::AttributeNames::forceX, ParticleBaseFP64::AttributeNames::forceY,
                    ParticleBaseFP64::AttributeNames::forceZ>(0, f);

  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::id>(0), 1);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 4.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posY>(0), 5.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::posZ>(0), 6.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceX>(0), 7.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceY>(0), 8.);
  EXPECT_EQ(soa.read<ParticleBaseFP64::AttributeNames::forceZ>(0), 9.);

  auto res = soa.readMultiple<ParticleBaseFP64::AttributeNames::forceX, ParticleBaseFP64::AttributeNames::forceY,
                              ParticleBaseFP64::AttributeNames::forceZ>(0);

  static_assert(std::is_same<decltype(res), std::array<double, 3>>::value, "id type must be proper(size_t)");

  EXPECT_EQ(res[0], 7.);
  EXPECT_EQ(res[1], 8.);
  EXPECT_EQ(res[2], 9.);
}

/**
 * Tests dynamic extra columns across resizeArraysForAttributes, write, read, begin, swap, pop_back, push, and clear.
 */
TEST_F(SoATest, SoATestExtraColumns) {
  using autopas::ParticleBaseFP64;
  using SoAType = autopas::SoA<ParticleBaseFP64::SoAArraysType>;
  constexpr size_t NumBase = SoAType::NumBaseAttrs;
  constexpr size_t Extra0 = NumBase + 0;
  constexpr size_t Extra1 = NumBase + 1;

  SoAType soa;
  EXPECT_EQ(soa.size(), 0);

  // Resize for posX, posY, posZ and two extra columns
  soa.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, ParticleBaseFP64::AttributeNames::posY,
                                ParticleBaseFP64::AttributeNames::posZ, Extra0, Extra1>(3);

  EXPECT_EQ(soa.size(), 3);

  // Write values to base and extra columns
  soa.write<ParticleBaseFP64::AttributeNames::posX>(0, 1.0);
  soa.write<Extra0>(0, 10.0);
  soa.write<Extra1>(0, 20.0);

  soa.write<ParticleBaseFP64::AttributeNames::posX>(1, 2.0);
  soa.write<Extra0>(1, 30.0);
  soa.write<Extra1>(1, 40.0);

  soa.write<ParticleBaseFP64::AttributeNames::posX>(2, 3.0);
  soa.write<Extra0>(2, 50.0);
  soa.write<Extra1>(2, 60.0);

  // Read via read()
  EXPECT_DOUBLE_EQ(soa.read<Extra0>(0), 10.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra1>(0), 20.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra0>(1), 30.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra1>(1), 40.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra0>(2), 50.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra1>(2), 60.0);

  // Test pointer access via begin()
  auto *extra0Ptr = soa.begin<Extra0>();
  auto *extra1Ptr = soa.begin<Extra1>();
  EXPECT_DOUBLE_EQ(extra0Ptr[0], 10.0);
  EXPECT_DOUBLE_EQ(extra1Ptr[0], 20.0);
  EXPECT_DOUBLE_EQ(extra0Ptr[1], 30.0);
  EXPECT_DOUBLE_EQ(extra1Ptr[1], 40.0);

  // Test swap
  soa.swap(0, 1);
  EXPECT_DOUBLE_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(0), 2.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra0>(0), 30.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra1>(0), 40.0);
  EXPECT_DOUBLE_EQ(soa.read<ParticleBaseFP64::AttributeNames::posX>(1), 1.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra0>(1), 10.0);
  EXPECT_DOUBLE_EQ(soa.read<Extra1>(1), 20.0);

  // Test pop_back
  soa.pop_back();
  EXPECT_EQ(soa.size(), 2);

  // Test push to extra column
  soa.push<Extra0>(70.0);
  EXPECT_EQ(soa.size(), 3);
  EXPECT_DOUBLE_EQ(soa.read<Extra0>(2), 70.0);

  // Test clear
  soa.clear();
  EXPECT_EQ(soa.size(), 0);
}

/**
 * Tests selective attribute resizing via resizeArraysForAttributes() followed by SoA::append(const SoA &).
 */
TEST_F(SoATest, SoATestSelectiveResizingAndAppend) {
  using autopas::ParticleBaseFP64;
  using SoAType = autopas::SoA<ParticleBaseFP64::SoAArraysType>;
  constexpr size_t NumBase = SoAType::NumBaseAttrs;
  constexpr size_t Extra0 = NumBase + 0;

  SoAType soa1;
  soa1.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, Extra0>(2);
  EXPECT_EQ(soa1.size(), 2);

  soa1.write<ParticleBaseFP64::AttributeNames::posX>(0, 1.5);
  soa1.write<Extra0>(0, 3.5);
  soa1.write<ParticleBaseFP64::AttributeNames::posX>(1, 2.5);
  soa1.write<Extra0>(1, 4.5);

  SoAType soa2;
  soa2.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, Extra0>(1);
  soa2.write<ParticleBaseFP64::AttributeNames::posX>(0, 5.5);
  soa2.write<Extra0>(0, 6.5);

  soa1.append(soa2);
  EXPECT_EQ(soa1.size(), 3);
  EXPECT_DOUBLE_EQ(soa1.read<ParticleBaseFP64::AttributeNames::posX>(2), 5.5);
  EXPECT_DOUBLE_EQ(soa1.read<Extra0>(2), 6.5);
}

/**
 * Tests SoAView::size<attribute>() on partially resized columns and SoA::append(const SoAView &) without reading
 * unresized columns.
 */
TEST_F(SoATest, SoATestSelectiveResizingAndViewAppend) {
  using autopas::ParticleBaseFP64;
  using SoAType = autopas::SoA<ParticleBaseFP64::SoAArraysType>;
  constexpr size_t NumBase = SoAType::NumBaseAttrs;
  constexpr size_t Extra0 = NumBase + 0;

  SoAType soa1;
  soa1.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, Extra0>(2);
  EXPECT_EQ(soa1.size(), 2);

  soa1.write<ParticleBaseFP64::AttributeNames::posX>(0, 1.5);
  soa1.write<Extra0>(0, 3.5);
  soa1.write<ParticleBaseFP64::AttributeNames::posX>(1, 2.5);
  soa1.write<Extra0>(1, 4.5);

  SoAType soa2;
  soa2.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, Extra0>(2);
  soa2.write<ParticleBaseFP64::AttributeNames::posX>(0, 5.5);
  soa2.write<Extra0>(0, 6.5);
  soa2.write<ParticleBaseFP64::AttributeNames::posX>(1, 7.5);
  soa2.write<Extra0>(1, 8.5);

  auto view = soa2.constructView(0, 1);
  EXPECT_EQ(view.size(), 1);
  EXPECT_EQ(view.size<ParticleBaseFP64::AttributeNames::posX>(), 1);
  EXPECT_EQ(view.size<ParticleBaseFP64::AttributeNames::posY>(), 0);
  EXPECT_EQ(view.size<Extra0>(), 1);

  // Appending the view must not read from unresized columns (e.g. posY, ptr)
  soa1.append(view);
  EXPECT_EQ(soa1.size(), 3);
  EXPECT_DOUBLE_EQ(soa1.read<ParticleBaseFP64::AttributeNames::posX>(2), 5.5);
  EXPECT_DOUBLE_EQ(soa1.read<Extra0>(2), 6.5);
}

/**
 * Verifies that read(), write(), readMultiple(), and writeMultiple() throw std::out_of_range for unallocated extra
 * columns.
 */
TEST_F(SoATest, SoATestUnallocatedExtraColumnAccessThrows) {
  using autopas::ParticleBaseFP64;
  using SoAType = autopas::SoA<ParticleBaseFP64::SoAArraysType>;
  constexpr size_t NumBase = SoAType::NumBaseAttrs;
  constexpr size_t Extra0 = NumBase + 0;
  constexpr size_t Extra1 = NumBase + 1;

  SoAType soa;
  soa.resizeArrays(2);

  // Extra columns were not allocated, read and write must throw std::out_of_range
  EXPECT_THROW(soa.read<Extra0>(0), std::out_of_range);
  EXPECT_THROW(soa.write<Extra0>(0, 1.0), std::out_of_range);
  EXPECT_THROW((soa.readMultiple<Extra0>(0)), std::out_of_range);
  EXPECT_THROW((soa.writeMultiple<Extra0>(0, std::array<double, 1>{1.0})), std::out_of_range);

  // Even if Extra0 is allocated, Extra1 is not
  soa.resizeArraysForAttributes<Extra0>(2);
  EXPECT_NO_THROW(soa.write<Extra0>(0, 2.5));
  EXPECT_DOUBLE_EQ(soa.read<Extra0>(0), 2.5);

  EXPECT_THROW(soa.read<Extra1>(0), std::out_of_range);
  EXPECT_THROW(soa.write<Extra1>(0, 3.5), std::out_of_range);
  EXPECT_THROW((soa.readMultiple<Extra1>(0)), std::out_of_range);
  EXPECT_THROW((soa.writeMultiple<Extra1>(0, std::array<double, 1>{3.5})), std::out_of_range);
}

/**
 * Tests std::span, size, pointer, alignment, and const accessors on sparsely allocated extra columns.
 */
TEST_F(SoATest, SoATestExtraColumnAccessors) {
  using autopas::ParticleBaseFP64;
  using SoAType = autopas::SoA<ParticleBaseFP64::SoAArraysType>;
  constexpr size_t NumBase = SoAType::NumBaseAttrs;
  constexpr size_t Extra1 = NumBase + 1;

  SoAType soa;
  EXPECT_EQ(soa.numExtraColumns(), 0);
  EXPECT_TRUE(soa.extraColumn(0).empty());
  EXPECT_EQ(soa.extraColumnSize(0), 0);
  EXPECT_EQ(soa.extraColumnBegin(0), nullptr);

  // Allocating Extra1 leaves extra column 0 empty while growing numExtraColumns() to 2.
  soa.resizeArraysForAttributes<Extra1>(3);
  EXPECT_EQ(soa.numExtraColumns(), 2);

  EXPECT_TRUE(soa.extraColumn(0).empty());
  EXPECT_EQ(soa.extraColumnSize(0), 0);
  EXPECT_EQ(soa.extraColumnBegin(0), nullptr);

  EXPECT_EQ(soa.extraColumnSize(1), 3);
  ASSERT_NE(soa.extraColumnBegin(1), nullptr);
  EXPECT_EQ(reinterpret_cast<uintptr_t>(soa.extraColumnBegin(1)) % autopas::DEFAULT_CACHE_LINE_SIZE, 0);

  // Out-of-bounds extra-column index must safely return empty / 0 / nullptr.
  EXPECT_TRUE(soa.extraColumn(2).empty());
  EXPECT_EQ(soa.extraColumnSize(2), 0);
  EXPECT_EQ(soa.extraColumnBegin(2), nullptr);
  EXPECT_EQ(std::as_const(soa).extraColumnBegin(2), nullptr);

  double *mutCol = soa.extraColumnBegin(1);
  mutCol[0] = 10.5;
  mutCol[1] = 20.5;
  mutCol[2] = 30.5;

  const auto constSpan = std::as_const(soa).extraColumn(1);
  ASSERT_EQ(constSpan.size(), 3);
  EXPECT_DOUBLE_EQ(constSpan[0], 10.5);
  EXPECT_DOUBLE_EQ(constSpan[1], 20.5);
  EXPECT_DOUBLE_EQ(constSpan[2], 30.5);
  EXPECT_EQ(std::as_const(soa).extraColumnBegin(1), constSpan.data());
  EXPECT_EQ(std::as_const(soa).begin<Extra1>(), constSpan.data());

  auto view = soa.constructView();
  EXPECT_EQ(view.numExtraColumns(), 2);
  EXPECT_TRUE(view.extraColumn(0).empty());
  EXPECT_EQ(view.extraColumnBegin(0), nullptr);

  // Non-const extraColumnBegin on a view must mutate the underlying SoA storage.
  view.extraColumnBegin(1)[1] = 25.5;
  const auto viewSpan = std::as_const(view).extraColumn(1);
  ASSERT_EQ(viewSpan.size(), 3);
  EXPECT_DOUBLE_EQ(viewSpan[1], 25.5);
  EXPECT_DOUBLE_EQ(std::as_const(view).begin<Extra1>()[1], 25.5);
  EXPECT_DOUBLE_EQ(soa.read<Extra1>(1), 25.5);
}

/**
 * Tests SoAView extra-column slicing with default views, non-zero start indices, partial columns, and append(view).
 */
TEST_F(SoATest, SoATestViewExtraColumnSlicing) {
  using autopas::ParticleBaseFP64;
  using SoAType = autopas::SoA<ParticleBaseFP64::SoAArraysType>;
  constexpr size_t NumBase = SoAType::NumBaseAttrs;
  constexpr size_t Extra1 = NumBase + 1;
  constexpr size_t Extra2 = NumBase + 2;

  const autopas::SoAView<ParticleBaseFP64::SoAArraysType> defaultView;
  EXPECT_EQ(defaultView.numExtraColumns(), 0);
  EXPECT_TRUE(defaultView.extraColumn(0).empty());
  EXPECT_EQ(defaultView.extraColumnSize(0), 0);
  EXPECT_EQ(defaultView.extraColumnBegin(0), nullptr);

  SoAType src;
  src.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, Extra1>(4);
  src.resizeSingleAttribute<Extra2>(2);
  for (size_t i = 0; i < 4; ++i) {
    src.write<ParticleBaseFP64::AttributeNames::posX>(i, static_cast<double>(i + 1));
    src.write<Extra1>(i, static_cast<double>((i + 1) * 10));
  }
  src.write<Extra2>(0, 100.0);
  src.write<Extra2>(1, 200.0);

  // Subslice [1, 3) has non-zero startIndex and overlaps partial column Extra2 (length 2).
  auto view = src.constructView(1, 3);
  EXPECT_EQ(view.size(), 2);
  EXPECT_TRUE(view.extraColumn(0).empty());
  EXPECT_EQ(view.extraColumnSize(0), 0);
  EXPECT_EQ(view.extraColumnBegin(0), nullptr);

  const auto col1Span = view.extraColumn(1);
  ASSERT_EQ(col1Span.size(), 2);
  EXPECT_DOUBLE_EQ(col1Span[0], 20.0);
  EXPECT_DOUBLE_EQ(col1Span[1], 30.0);
  EXPECT_EQ(view.extraColumnSize(1), 2);
  EXPECT_DOUBLE_EQ(std::as_const(view).begin<Extra1>()[0], 20.0);
  EXPECT_DOUBLE_EQ(std::as_const(view).begin<Extra1>()[1], 30.0);

  const auto col2Span = view.extraColumn(2);
  ASSERT_EQ(col2Span.size(), 1);
  EXPECT_DOUBLE_EQ(col2Span[0], 200.0);
  EXPECT_EQ(view.size<Extra2>(), 1);

  // Slice starting past Extra2's allocated size must expose an empty span and nullptr begin.
  auto tailView = src.constructView(2, 4);
  EXPECT_EQ(tailView.size<Extra2>(), 0);
  EXPECT_TRUE(tailView.extraColumn(2).empty());
  EXPECT_EQ(tailView.extraColumnBegin(2), nullptr);
  ASSERT_EQ(tailView.extraColumn(1).size(), 2);
  EXPECT_DOUBLE_EQ(tailView.extraColumn(1)[0], 30.0);
  EXPECT_DOUBLE_EQ(tailView.extraColumn(1)[1], 40.0);

  // Zero-length slice at non-zero index must yield empty extra-column spans and nullptr begin.
  auto emptySlice = src.constructView(2, 2);
  EXPECT_EQ(emptySlice.size(), 0);
  EXPECT_TRUE(emptySlice.extraColumn(1).empty());
  EXPECT_EQ(emptySlice.extraColumnBegin(1), nullptr);

  SoAType dst;
  dst.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, Extra1>(1);
  dst.write<ParticleBaseFP64::AttributeNames::posX>(0, 9.0);
  dst.write<Extra1>(0, 5.0);

  dst.append(defaultView);
  dst.append(emptySlice);
  EXPECT_EQ(dst.size(), 1);
  EXPECT_EQ(dst.numExtraColumns(), 2);

  // Appending [1, 3) grows dst._extraColumns for Extra2 while skipping empty Extra0.
  dst.append(view);
  EXPECT_EQ(dst.size(), 3);
  EXPECT_EQ(dst.numExtraColumns(), 3);
  EXPECT_DOUBLE_EQ(dst.read<ParticleBaseFP64::AttributeNames::posX>(0), 9.0);
  EXPECT_DOUBLE_EQ(dst.read<ParticleBaseFP64::AttributeNames::posX>(1), 2.0);
  EXPECT_DOUBLE_EQ(dst.read<ParticleBaseFP64::AttributeNames::posX>(2), 3.0);
  EXPECT_TRUE(dst.extraColumn(0).empty());
  ASSERT_EQ(dst.extraColumnSize(1), 3);
  EXPECT_DOUBLE_EQ(dst.read<Extra1>(0), 5.0);
  EXPECT_DOUBLE_EQ(dst.read<Extra1>(1), 20.0);
  EXPECT_DOUBLE_EQ(dst.read<Extra1>(2), 30.0);
  ASSERT_EQ(dst.extraColumnSize(2), 1);
  EXPECT_DOUBLE_EQ(dst.read<Extra2>(0), 200.0);
}

/**
 * Tests templated append<attributes...>() on SoA and SoAView sources (exercising appendSingleArray), including
 * count == 0.
 */
TEST_F(SoATest, SoATestTemplatedAppendExtraColumns) {
  using autopas::ParticleBaseFP64;
  using SoAType = autopas::SoA<ParticleBaseFP64::SoAArraysType>;
  constexpr size_t NumBase = SoAType::NumBaseAttrs;
  constexpr size_t Extra0 = NumBase + 0;
  constexpr size_t Extra1 = NumBase + 1;

  SoAType src;
  src.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, Extra1>(3);
  src.write<ParticleBaseFP64::AttributeNames::posX>(0, 2.0);
  src.write<ParticleBaseFP64::AttributeNames::posX>(1, 3.0);
  src.write<ParticleBaseFP64::AttributeNames::posX>(2, 4.0);
  src.write<Extra1>(0, 22.0);
  src.write<Extra1>(1, 33.0);
  src.write<Extra1>(2, 44.0);

  SoAType dstFromSoA;
  dstFromSoA.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX, Extra1>(1);
  dstFromSoA.write<ParticleBaseFP64::AttributeNames::posX>(0, 1.0);
  dstFromSoA.write<Extra1>(0, 11.0);

  // Extra0 has size 0 in src, so append<..., Extra0, ...> must skip it cleanly.
  dstFromSoA.append<ParticleBaseFP64::AttributeNames::posX, Extra0, Extra1>(src);
  EXPECT_EQ(dstFromSoA.size(), 4);
  EXPECT_EQ(dstFromSoA.size<Extra0>(), 0);
  EXPECT_TRUE(dstFromSoA.extraColumn(0).empty());
  ASSERT_EQ(dstFromSoA.size<ParticleBaseFP64::AttributeNames::posX>(), 4);
  ASSERT_EQ(dstFromSoA.size<Extra1>(), 4);
  EXPECT_DOUBLE_EQ(dstFromSoA.read<ParticleBaseFP64::AttributeNames::posX>(0), 1.0);
  EXPECT_DOUBLE_EQ(dstFromSoA.read<ParticleBaseFP64::AttributeNames::posX>(3), 4.0);
  EXPECT_DOUBLE_EQ(dstFromSoA.read<Extra1>(0), 11.0);
  EXPECT_DOUBLE_EQ(dstFromSoA.read<Extra1>(3), 44.0);

  SoAType dstFromView;
  dstFromView.resizeArraysForAttributes<ParticleBaseFP64::AttributeNames::posX>(1);
  dstFromView.write<ParticleBaseFP64::AttributeNames::posX>(0, 5.0);

  auto srcView = src.constructView(1, 3);
  dstFromView.append<ParticleBaseFP64::AttributeNames::posX, Extra0, Extra1>(srcView);
  EXPECT_EQ(dstFromView.size(), 3);
  EXPECT_EQ(dstFromView.size<Extra0>(), 0);
  EXPECT_TRUE(dstFromView.extraColumn(0).empty());
  EXPECT_DOUBLE_EQ(dstFromView.read<ParticleBaseFP64::AttributeNames::posX>(0), 5.0);
  EXPECT_DOUBLE_EQ(dstFromView.read<ParticleBaseFP64::AttributeNames::posX>(1), 3.0);
  EXPECT_DOUBLE_EQ(dstFromView.read<ParticleBaseFP64::AttributeNames::posX>(2), 4.0);
  ASSERT_EQ(dstFromView.size<Extra1>(), 2);
  EXPECT_DOUBLE_EQ(dstFromView.read<Extra1>(0), 33.0);
  EXPECT_DOUBLE_EQ(dstFromView.read<Extra1>(1), 44.0);

  // Empty view must leave destination untouched.
  dstFromView.append<ParticleBaseFP64::AttributeNames::posX, Extra1>(src.constructView(1, 1));
  EXPECT_EQ(dstFromView.size(), 3);
  EXPECT_EQ(dstFromView.size<Extra1>(), 2);
}

/**
 * Tests push, resizeSingleAttribute, resizeArrays, writeMultiple/readMultiple, clear, and empty-pack resize on extra
 * columns.
 */
TEST_F(SoATest, SoATestExtraColumnLifecycle) {
  using autopas::ParticleBaseFP64;
  using SoAType = autopas::SoA<ParticleBaseFP64::SoAArraysType>;
  constexpr size_t NumBase = SoAType::NumBaseAttrs;
  constexpr size_t Extra0 = NumBase + 0;
  constexpr size_t Extra1 = NumBase + 1;
  constexpr size_t Extra2 = NumBase + 2;

  SoAType soa;
  soa.push<Extra1>(12.5);
  EXPECT_EQ(soa.size(), 1);
  EXPECT_EQ(soa.numExtraColumns(), 2);
  EXPECT_EQ(soa.size<Extra0>(), 0);
  EXPECT_EQ(soa.size<Extra1>(), 1);
  EXPECT_DOUBLE_EQ(soa.read<Extra1>(0), 12.5);

  soa.resizeSingleAttribute<Extra2>(2);
  EXPECT_EQ(soa.numExtraColumns(), 3);
  EXPECT_EQ(soa.size<Extra2>(), 2);

  // resizeArrays must resize both base storage and all already-allocated extra columns.
  soa.resizeArrays(3);
  EXPECT_EQ(soa.size(), 3);
  EXPECT_EQ(soa.size<Extra0>(), 3);
  EXPECT_EQ(soa.size<Extra1>(), 3);
  EXPECT_EQ(soa.size<Extra2>(), 3);
  EXPECT_DOUBLE_EQ(soa.read<Extra1>(0), 12.5);

  soa.writeMultiple<ParticleBaseFP64::AttributeNames::posX, Extra0, Extra2>(2, std::array<double, 3>{7.5, 8.5, 9.5});
  const auto values = soa.readMultiple<ParticleBaseFP64::AttributeNames::posX, Extra0, Extra2>(2);
  EXPECT_DOUBLE_EQ(values[0], 7.5);
  EXPECT_DOUBLE_EQ(values[1], 8.5);
  EXPECT_DOUBLE_EQ(values[2], 9.5);

  soa.clear();
  EXPECT_EQ(soa.size(), 0);
  EXPECT_EQ(soa.size<Extra0>(), 0);
  EXPECT_EQ(soa.size<Extra1>(), 0);
  EXPECT_EQ(soa.size<Extra2>(), 0);

  const size_t prevExtraCols = soa.numExtraColumns();
  soa.resizeArraysForAttributes<>(4);
  EXPECT_EQ(soa.size(), 4);
  EXPECT_EQ(soa.numExtraColumns(), prevExtraCols);
}

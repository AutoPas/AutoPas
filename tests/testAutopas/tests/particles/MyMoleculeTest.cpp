/**
 * @file myMoleculeTest.cpp
 * @author seckler
 * @date 18.01.18
 */

#include <gtest/gtest.h>

#include "autopas/particles/ParticleDefinitions.h"
#include "tests/autopasInterface/AutoPasTest.h"

using namespace autopas;

class MyMolecule : public ParticleFP64 {
 public:
  MyMolecule() : ParticleFP64(), _myvar(0) {}
  MyMolecule(const std::array<double, 3> &r, const std::array<double, 3> &v, unsigned long i, int myvar)
      : ParticleFP64(r, v, i), _myvar(myvar) {}
  void print() {
    std::cout << "Molecule with position: ";
    for (auto &r : getR()) {
      std::cout << r << ", ";
    }
    std::cout << "and force: ";

    for (auto &f : getF()) {
      std::cout << f << ", ";
    }
    std::cout << "ID: " << getID();
    std::cout << " myvar: " << _myvar << std::endl;
  }

  int getMyvar() const { return _myvar; }

 private:
  int _myvar;
};

class MyMoleculeTest : public AutoPasTest {};

TEST_F(MyMoleculeTest, testConstructorAndGetters) {
  std::array<double, 3> r({1.1, 2.2, 3.3});
  int myvar = 5;
  std::array<double, 3> vel({4.4, 5.5, 6.6});
  unsigned long id = 17ul;
  MyMolecule m(r, vel, id, myvar);

  for (int d = 0; d < 3; ++d) {
    ASSERT_DOUBLE_EQ(m.getR()[d], r[d]);
    ASSERT_DOUBLE_EQ(m.getV()[d], vel[d]);
  }
  ASSERT_EQ(id, m.getID());
  ASSERT_EQ(myvar, m.getMyvar());
}

TEST_F(MyMoleculeTest, testMovingParticles) {
  const std::array<double, 3> zero{{0.0, 0.0, 0.0}};
  const double maxAllowedDistSquared = 1;
  {
    MyMolecule m(zero, zero, 0, 0);
    EXPECT_TRUE(m.setRDistanceCheck({0.5, 0.0, 0.0}, maxAllowedDistSquared));
  }
  {
    MyMolecule m(zero, zero, 0, 0);
    EXPECT_FALSE(m.setRDistanceCheck({1.5, 0.0, 0.0}, maxAllowedDistSquared));
  }
  {
    MyMolecule m(zero, zero, 0, 0);
    EXPECT_TRUE(m.addRDistanceCheck({0.5, 0.0, 0.0}, maxAllowedDistSquared));
  }
  {
    MyMolecule m(zero, zero, 0, 0);
    EXPECT_FALSE(m.addRDistanceCheck({1.5, 0.0, 0.0}, maxAllowedDistSquared));
  }
}
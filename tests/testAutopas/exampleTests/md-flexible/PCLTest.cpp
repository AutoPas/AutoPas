
#include "PCLTest.h"
using namespace std;

double PCLTest::mixingE(double e1, double e2) { return std::sqrt(e1 * e2); }
double PCLTest::mixingS(double s1, double s2) { return ((s1 + s2) / 2); }
TEST_F(PCLTest, Functions) {
  // Testing PCL function with default Epsilon24 and Sigma values
  ASSERT_EQ(mass, PCL.getMass(dummyParticle.getTypeId()));
  PrintableMolecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
  PrintableMolecule p2({0., 0., 0.}, {0., 0., 0.}, 1);
  ASSERT_EQ(PCL.mixing24E(p1.getTypeId(), p2.getTypeId()), 24 * mixingE(epsilon, epsilon));
  ASSERT_EQ(PCL.mixingSS(p1.getTypeId(), p2.getTypeId()),
            mixingS(sigma,sigma)*mixingS(sigma,sigma));
}
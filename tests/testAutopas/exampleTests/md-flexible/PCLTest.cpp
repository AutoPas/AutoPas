
#include "PCLTest.h"
using namespace std;

double PCLTest::mixingE(double e1, double e2) { return std::sqrt(e1 * e2); }
double PCLTest::mixingS(double s1, double s2) { return ((s1 + s2) / 2); }
TEST_F(PCLTest, Functions) {
  // Testing PCL function with default Epsilon and Sigma values
  ASSERT_EQ(1.0, PCL.getMass(dummyParticle));
  PrintableMolecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
  PrintableMolecule p2({0., 0., 0.}, {0., 0., 0.}, 1);
  ASSERT_EQ(PCL.mixing24E(p1.getID(), p2.getID()), 24 * mixingE(p1.getEpsilon(), p2.getEpsilon()));
  ASSERT_EQ(PCL.mixingSS(p1.getID(), p2.getID()),
            mixingS(p1.getSigma(), p2.getSigma()) * mixingS(p1.getSigma(), p2.getSigma()));
}
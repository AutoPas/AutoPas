
#include "PCLTest.h"
using namespace std;

double PCLTest::mixingE(double e1, double e2) { return std::sqrt(e1 * e2); }
double PCLTest::mixingS(double s1, double s2) { return ((s1 + s2) / 2); }


TEST_F(PCLTest, Functions) {
   PCL.addType(0,epsilon,sigma,mass);PCL.addType(1,epsilon2,sigma2,mass);
  PrintableMolecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
  PrintableMolecule p2({0., 0., 0.}, {0., 0., 0.}, 1);
    p2.setTypeId(1);
ASSERT_EQ(mass, PCL.getMass(p1.getTypeId()));
    ASSERT_EQ(PCL.mixing24Epsilon(p1.getTypeId(), p2.getTypeId()), 24 * mixingE(epsilon, epsilon2));
    ASSERT_EQ(PCL.mixingSigmaSquare(p1.getTypeId(), p2.getTypeId()),
            mixingS(sigma,sigma2)*mixingS(sigma,sigma2));
}
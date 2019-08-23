
#include "ObjectsTest.h"

TEST_F(ObjectsTest, BoxMin_Max) {
  // as initilized in constructor of ObjectTests: BoxMin or BoxMax are eather equal boxlength or velocity array
  std::array<double, 3> SphereBoxMax = {5., 5., 5.};
  std::array<double, 3> SphereBoxMin = {-5., -5., -5.};
  EXPECT_EQ(_CGrid.getBoxMax(), boxlength);
  EXPECT_EQ(_CGauss.getBoxMax(), boxlength);
  EXPECT_EQ(_CUniform.getBoxMax(), boxlength);
  EXPECT_EQ(_Sphere.getBoxMax(), SphereBoxMax);
  EXPECT_EQ(_CGrid.getBoxMin(), velocity);
  EXPECT_EQ(_CGauss.getBoxMin(), velocity);
  EXPECT_EQ(_CUniform.getBoxMin(), velocity);
  EXPECT_EQ(_Sphere.getBoxMin(), SphereBoxMin);
}


#include "ObjectsTest.h"

TEST_F(ObjectsTest, BoxMin_Max) {
  // as initilized in constructor of ObjectTests: BoxMin or BoxMax are eather equal boxlength or velocity array
  EXPECT_EQ(_CGrid.getBoxMax(), boxlength);
  EXPECT_EQ(_CGauss.getBoxMax(), boxlength);
  EXPECT_EQ(_CUniform.getBoxMax(), boxlength);
  EXPECT_EQ(_Sphere.getBoxMax(), boxlength);
  EXPECT_EQ(_CGrid.getBoxMin(), velocity);
  EXPECT_EQ(_CGauss.getBoxMin(), velocity);
  EXPECT_EQ(_CUniform.getBoxMin(), velocity);
  EXPECT_EQ(_Sphere.getBoxMin(), velocity);
}

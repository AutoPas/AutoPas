/**
 * @file LJFunctorCudaTest.cu
 * @author jspahl
 * @date 14.12.18
 */

#include "LJFunctorCudaTest.h"
#include <autopas/particles/MoleculeLJ.h>
#include "testingHelpers/commonTypedefs.h"

void LJFunctorCudaTest::testAoSNoGlobals(bool newton3) {
  autopas::LJFunctorCuda<Molecule, FMCell> functor(cutoff, epsilon, sigma, shift);

  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1);

  double3 rs [2]= {{0., 0., 0.}, {0.1, 0.2, 0.3}};
  double3 fs [2]= {{0., 0., 0.}, {0.0, 0.0, 0.0}};

  double3* d_cell0_r;
  double3* d_cell0_f;

  cudaError_t e;
  //e = cudaMalloc((void **) &d_cell0_r, sizeof(double3) * 2);
  //e = cudaMalloc((void **) &d_cell0_f, sizeof(double3) * 2);

  //e = cudaMemcpy(d_cell0_r, rs, 2, cudaMemcpyHostToDevice);
  //e = cudaMemcpy(d_cell0_f, fs, 2, cudaMemcpyHostToDevice);

  functor.CudaFunctor(2, d_cell0_r, d_cell0_f, false);

  //e = cudaMemcpy(fs, d_cell0_f, 2, cudaMemcpyDeviceToHost);
  //cudaFree(d_cell0_r);
  //cudaFree(d_cell0_f);

  auto f1one = p1.getF();
  auto f2one = p2.getF();
  EXPECT_NEAR(f1one[0], expectedForce[0], absDelta);
  EXPECT_NEAR(f1one[1], expectedForce[1], absDelta);
  EXPECT_NEAR(f1one[2], expectedForce[2], absDelta);

  if (newton3) {
    EXPECT_NEAR(f2one[0], -expectedForce[0], absDelta);
    EXPECT_NEAR(f2one[1], -expectedForce[1], absDelta);
    EXPECT_NEAR(f2one[2], -expectedForce[2], absDelta);
  } else {
    EXPECT_DOUBLE_EQ(f2one[0], 0);
    EXPECT_DOUBLE_EQ(f2one[1], 0);
    EXPECT_DOUBLE_EQ(f2one[2], 0);
  }

  //functor.AoSFunctor(p2, p1, newton3);

  auto f1two = p1.getF();
  auto f2two = p2.getF();

  double factor = newton3 ? 2. : 1.;

  EXPECT_NEAR(f1two[0], factor * expectedForce[0], absDelta);
  EXPECT_NEAR(f1two[1], factor * expectedForce[1], absDelta);
  EXPECT_NEAR(f1two[2], factor * expectedForce[2], absDelta);

  EXPECT_NEAR(f2two[0], -factor * expectedForce[0], absDelta);
  EXPECT_NEAR(f2two[1], -factor * expectedForce[1], absDelta);
  EXPECT_NEAR(f2two[2], -factor * expectedForce[2], absDelta);
}

TEST_F(LJFunctorCudaTest, testAoSFunctorNoGlobalsNoN3) {
  bool newton3 = false;
  testAoSNoGlobals(newton3);
}


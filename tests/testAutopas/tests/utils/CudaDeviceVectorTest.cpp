/*
 * @file CudaDeviceVectorTest.cpp
 * @author: jspahl
 * @date 18.2.19
 */

#if defined(AUTOPAS_CUDA)

#include "tests/utils/CudaDeviceVectorTest.h"

#include <vector>

#include "autopas/utils/CudaDeviceVector.h"

namespace CudaDeviceVectorTest {

TEST_F(CudaDeviceVectorTest, CopyTest) {
  std::vector<int> source = {1, -2, 3, -4, 5, -6, 7};
  autopas::utils::CudaDeviceVector<int> test_vector;
  test_vector.copyHostToDevice(7, source.data());

  std::vector<int> target = {0, 0, 0, 0, 0, 0, 0};
  test_vector.copyDeviceToHost(7, target.data());

  for (size_t i = 0; i < source.size(); ++i) {
    EXPECT_EQ(source[i], target[i]) << "Vectors differ at index " << i;
  }
}

}  // end namespace CudaDeviceVectorTest
#endif

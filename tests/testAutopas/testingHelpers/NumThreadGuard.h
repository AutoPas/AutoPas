/**
 * @file NumThreadGuard.h
 * @author C. Menges
 * @date 14.04.2019
 */

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

class NumThreadGuard {
 public:
  NumThreadGuard(int newNum) {
#ifdef AUTOPAS_OPENMP
    numThreadsBefore = omp_get_max_threads();
    omp_set_num_threads(newNum);
#endif
  }

  ~NumThreadGuard() {
#ifdef AUTOPAS_OPENMP
    omp_set_num_threads(numThreadsBefore);
#endif
  }

 private:
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore;
#endif
};

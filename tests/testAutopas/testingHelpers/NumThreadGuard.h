/**
 * @file NumThreadGuard.h
 * @author C. Menges
 * @date 14.04.2019
 */

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

/**
 * NumThreadGuard sets current number of threads to newNum and resets number of threads during destruction.
 */
class NumThreadGuard final {
 public:
  /**
   * Construct a new NumThreadGuard object and sets current number of threads to newNum.
   * @param newNum new number of threads
   */
  NumThreadGuard(const int newNum) {
#ifdef AUTOPAS_OPENMP
    numThreadsBefore = omp_get_max_threads();
    omp_set_num_threads(newNum);
#endif
  }

  /**
   * Destroy the NumThreadGuard object and reset number of threads.
   */
  ~NumThreadGuard() {
#ifdef AUTOPAS_OPENMP
    omp_set_num_threads(numThreadsBefore);
#endif
  }

  /**
   * delete copy constructor.
   */
  NumThreadGuard(const NumThreadGuard&) = delete;

  /**
   * delete copy assignment constructor.
   * @return deleted, so not important
   */
  NumThreadGuard& operator=(const NumThreadGuard&) = delete;

 private:
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore;
#endif
};

#include <rocm_smi/rocm_smi.h>

#include "ROCM.h"

namespace pmt::rocm {
class ROCMState {
 public:
  operator State();
  Timestamp timestamp_;
  double watt_ = 0;
  double joules_ = 0;
};

class ROCMImpl : public ROCM {
 public:
  ROCMImpl(const unsigned device_number);
  ~ROCMImpl();

 private:
  virtual State GetState() override;

  virtual const char *GetDumpFilename() override { return "/tmp/pmt_rocm.out"; }

  unsigned int device_number_;

  ROCMState state_previous_;
  ROCMState GetROCMState();
};
}  // end namespace pmt::rocm
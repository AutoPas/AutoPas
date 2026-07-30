#include <amd_smi/amdsmi.h>

#include "AMDSMI.h"

namespace pmt::amdsmi {

class AMDSMIState {
 public:
  operator State();
  Timestamp timestamp_;
  double watt_ = 0;
  double joules_ = 0;
};

class AMDSMIImpl : public AMDSMI {
 public:
  AMDSMIImpl(const unsigned device_number);
  ~AMDSMIImpl();

 private:
  virtual State GetState() override;

  virtual const char *GetDumpFilename() override {
    return "/tmp/pmt_amdsmi.out";
  }

  amdsmi_processor_handle processor_;

  AMDSMIState state_previous_;
  AMDSMIState GetAMDSMIState();
};

}  // end namespace pmt::amdsmi
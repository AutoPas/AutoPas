#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "Tegra.h"
#include "common/PMT.h"

namespace pmt::tegra {

using TegraMeasurement = std::pair<std::string, int>;

class TegraState {
 public:
  operator State();
  Timestamp timestamp_;
  std::vector<TegraMeasurement> measurements;
  unsigned int watt_ = 0;
  unsigned int joules_ = 0;
};

class TegraImpl : public Tegra {
 public:
  TegraImpl();
  virtual ~TegraImpl();

  State GetState() override { return GetTegraState(); }

  virtual const char *GetDumpFilename() override {
    return "/tmp/pmt_tegra.out";
  }

 private:
  TegraState GetTegraState();
  std::vector<TegraMeasurement> GetMeasurements();

  std::string filename_ = "";
  bool started_tegrastats_ = false;
  const int measurement_interval_ = 10;  // milliseconds
  TegraState previous_state_;
};

}  // end namespace pmt::tegra
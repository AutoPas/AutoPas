#ifndef PMT_RAPLIMPL_H_
#define PMT_RAPLIMPL_H_

#include <cstddef>
#include <mutex>
#include <string>
#include <vector>
#include <array>

#include "Rapl.h"
#include "common/PMT.h"
#include "common/io.h"

namespace pmt::rapl {

const int kNumRaplDomains = 4;
const int kKeepAliveInterval = 10;  // call Measure() roughly every nth update

struct RaplMeasurement {
  std::string name;
  std::size_t joules;
};

class RaplImpl : public Rapl {
 public:
  RaplImpl();

  State GetState() override;

  virtual const char *GetDumpFilename() override { return "/tmp/pmt_rapl.out"; }

 private:
  void Init();
  std::vector<int> fd_energy_uj_;

  std::vector<RaplMeasurement> GetMeasurements();

  Timestamp previous_timestamp_;
  std::vector<RaplMeasurement> previous_measurements_;

  std::vector<std::string> packages_names_;
  std::vector<os::file_descriptor> energy_fds_;

  // The numbers in the rapl /energy_uj files range from zero up to a maximum
  // specified in /max_energy_range_uj. This class reports monotonically
  // increasing values starting with zero for the first measurement. Therefore,
  // for every counter, the result is computed as follows:
  //  now = <read value>
  //  offset += now < previous ? max : 0
  //  result += offset + now - first
  std::vector<std::size_t> uj_max_;
  std::vector<std::size_t> uj_first_;
  std::vector<std::size_t> uj_previous_;
  std::vector<std::size_t> uj_offset_;

  // Mutex used to guard GetMeasurements()
  std::mutex mutex_;
};

}  // end namespace pmt::rapl

#endif  // PMT_RAPLIMPL_H_

#ifndef PMT_COMMON_H_
#define PMT_COMMON_H_

#include <chrono>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace pmt {

const std::string kDumpFilenameVariable = "PMT_DUMP_FILE";
const std::string kDumpIntervalVariable = "PMT_DUMP_INTERVAL";

using Timestamp = std::chrono::high_resolution_clock::time_point;

class State {
 public:
  State &operator=(const State &state) {
    timestamp_ = state.timestamp_;
    nr_measurements_ = state.nr_measurements_;
    name_ = state.name_;
    joules_ = state.joules_;
    watt_ = state.watt_;
    return *this;
  }

  State(int nr_measurements = 1) : nr_measurements_(nr_measurements) {
    name_.resize(nr_measurements);
    joules_.resize(nr_measurements);
    watt_.resize(nr_measurements);
  }

  int NrMeasurements() { return nr_measurements_; }

  Timestamp timestamp() { return timestamp_; }
  std::string name(int i);
  float joules(int i);
  float watts(int i);

  Timestamp timestamp_;
  int nr_measurements_;
  std::vector<std::string> name_;
  std::vector<float> joules_;
  std::vector<float> watt_;
};

class PMT {
 public:
  virtual ~PMT();

  static double seconds(const Timestamp &timestamp);

  static double seconds(const Timestamp &first, const Timestamp &second);

  static double seconds(const State &first, const State &second);

  static double joules(const State &first, const State &second);

  static double watts(const State &first, const State &second);

  void StartDump(const char *filename = nullptr);
  void StopDump();

  virtual void Mark(const State &start, const State &current,
                    const std::string &message) const;

  void SetMeasurementInterval(unsigned int milliseconds = 0);
  unsigned int GetMeasurementInterval() const {
    return measurement_interval_;
  };                               // in milliseconds
  unsigned int GetDumpInterval();  // in milliseconds

  State Read();

 protected:
  virtual State GetState() { return state_latest_; };

  virtual const char *GetDumpFilename() = 0;

  void Dump(const State &state);

  Timestamp GetTime();

 private:
  unsigned int measurement_interval_ = 100;  // milliseconds

  // The last state set by the thread
  State state_latest_;

  // This thread continuously call GetState to update state_latest_. It is
  // started automatically upon the first Read() call.
  std::thread thread_;
  volatile bool thread_started_ = false;
  volatile bool thread_stop_ = false;

  void StartThread();
  void StopThread();

  void DumpHeader(const State &state);

 protected:
  std::unique_ptr<std::ofstream> dump_file_ = nullptr;
  mutable std::mutex dump_file_mutex_;
};

std::unique_ptr<PMT> Create(const std::string &name,
                            const std::string &argument = "");

}  // end namespace pmt

std::ostream &operator<<(std::ostream &os, const pmt::Timestamp &timestamp);

#endif

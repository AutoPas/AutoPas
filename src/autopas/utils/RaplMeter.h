/**
 * @file RaplMeter.h
 * @author Vincent Fischer
 * @date 09.11.2021
 */

#pragma once

#include <cstdint>
#include <cstdio>
#include <fstream>

namespace autopas::utils {

/**
 * Measures the energy consumption between calls of reset() and sample()
 * using perf events. This requires root or a paranoid value less than 1.
 */
class RaplMeter {
 private:
  int _type;
  long _psys_raw, _pkg_raw, _cores_raw, _ram_raw;
  double _psys_unit, _pkg_unit, _cores_unit, _ram_unit;
  int _psys_config, _pkg_config, _cores_config, _ram_config;
  int _psys_fd, _pkg_fd, _cores_fd, _ram_fd;

  void open_msr();
  uint64_t read_msr(int msr_offset);
  int open_perf_event(int type, int config);
  long read_perf_event(int fd);

 public:
  RaplMeter();
  ~RaplMeter();
  /** reset perf file descriptors to start new measurement */
  void reset();
  /** measure power consumption since last call of reset
   * the results can be retrieved with the get_<domain>_energy() functions.
   */
  void sample();

  /** returns the energy consumed by the cpu package between the last call of sample() and the preceding
   * call of reset()
   * @return Energy in Joules
   */
  double get_pkg_energy();
  /** returns the energy consumed by the cpu cores between the last call of sample() and the preceding
   * call of reset()
   * @return Energy in Joules
   */
  double get_cores_energy();
  /** returns the energy consumed by ram between the last call of sample() and the preceding call of
   * reset()
   * @return Energy in Joules
   */
  double get_ram_energy();
  /** returns the energy consumed by the entire system between the last call of sample() and the preceding
   * call of reset()
   * @return Energy in Joules
   */
  double get_psys_energy();
  // unscaled for tuning
  /** return energy measurement for tuning purposes. Depending on which perf counters are available this
   * may return psys, or a sum of pkg and ram. Note that the units used here are not consistent across different
   * systems.
   * @return unscaled energy
   */
  long get_total_energy();
};

}  // namespace autopas::utils

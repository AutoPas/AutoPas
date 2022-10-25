/**
 * @file RaplMeter.h
 * @author Vincent Fischer
 * @date 09.11.2021
 */

#pragma once

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <vector>

namespace autopas::utils {

/**
 * Measures the energy consumption between calls of reset() and sample()
 * using perf events. This requires root or a paranoid value less than 1.
 */
class RaplMeter {
 private:
  std::vector<int> _cpus;
  int _type;
  long _psys_raw, _pkg_raw, _cores_raw, _ram_raw;
  double _psys_unit, _pkg_unit, _cores_unit, _ram_unit;
  int _psys_config, _pkg_config, _cores_config, _ram_config;
  std::vector<int> _psys_fd, _pkg_fd, _cores_fd, _ram_fd;

  int open_perf_event(int type, int config, int cpu);
  long read_perf_event(int fd);

 public:
  ~RaplMeter();

  /**
   *initialisation may fail, so moved out of constructor
   */
  void init();

  /**
   * reset perf file descriptors to start new measurement
   */
  void reset();

  /**
   * measure power consumption since last call to reset
   * the results can be retrieved with the get_<domain>_energy() functions.
   */
  void sample();

  /**
   * returns the energy consumed by the cpu package between the last call to sample() and the preceding
   * call to reset()
   * @return Energy in Joules
   */
  double get_pkg_energy();

  /**
   * returns the energy consumed by the cpu cores between the last call to sample() and the preceding
   * call to reset()
   * @return Energy in Joules
   */
  double get_cores_energy();

  /**
   *returns the energy consumed by ram between the last call to sample() and the preceding call to
   * reset()
   * @return Energy in Joules
   */
  double get_ram_energy();

  /**
   * returns the energy consumed by the entire system between the last call to sample() and the preceding
   * call to reset()
   * @return Energy in Joules
   */
  double get_psys_energy();

  /**
   * return energy measurement for tuning purposes. Depending on which perf counters are available this
   * may return psys, or a sum of pkg and ram. Note that the units used here are not consistent across different
   * systems. The scaling is skipped for this purpose, as the AutoTuner expects long values as opposed to doubles
   * and rounding the scaled value seems like the inferior option as well as unnecessary given that the values provided
   * by the powercap framework already have the correct type. The getters for pkg, cores, ram and psys provide scaled
   * values for comparisons between different systems.
   *
   * @return unscaled energy
   */
  long get_total_energy();
};

}  // namespace autopas::utils

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

class RaplMeter {
 private:
  int _type;
  long _psys_raw, _pkg_raw, _cores_raw, _ram_raw;
  double _psys_unit, _pkg_unit, _cores_unit, _ram_unit;
  int _psys_config, _pkg_config, _cores_config, _ram_config;
  int _psys_fd, _pkg_fd, _cores_fd, _ram_fd;

  void open_msr();
  uint64_t read_msr(int msr_offset);

 public:
  RaplMeter();
  ~RaplMeter();
  void reset();
  void sample();

  double get_pkg_energy();
  double get_cores_energy();
  double get_ram_energy();
  double get_psys_energy();
  // unscaled for tuning
  long get_total_energy();
};

}  // namespace autopas::utils

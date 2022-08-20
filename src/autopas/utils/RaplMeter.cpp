/**
 * @file RaplMeter.cpp
 * @author Vincent Fischer
 * @date 09.11.2021
 */

#include "autopas/utils/RaplMeter.h"

#include <linux/perf_event.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <cmath>
#include <cstring>
#include <sstream>

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas::utils {

int RaplMeter::open_perf_event(int type, int config) {
  struct perf_event_attr attr;
  memset(&attr, 0, sizeof(attr));
  attr.type = type;
  attr.config = config;
  if (config == 0) {
    return -1;
  }

  int fd = syscall(__NR_perf_event_open, &attr, -1 /*PID*/, 0 /*CPU*/, -1 /*GROUP FD*/, 0 /*FLAGS*/);
  if (fd < 0) {
    if (errno == EACCES) {
      throw ExceptionHandler::AutoPasException("Failed to open perf event: Permission denied");
    }
    throw ExceptionHandler::AutoPasException("Failed to open perf event");
  }
  return fd;
}

void RaplMeter::init() {
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/type", "r")) {
    if (fscanf(fff, "%d", &this->_type) != 1) {
      throw ExceptionHandler::AutoPasException("Failed to parse /sys/bus/event_source/devices/power/type");
    }
    fclose(fff);
  } else {
    throw ExceptionHandler::AutoPasException("No support for energy measurements detected.");
  }

  int psys_config;
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/events/energy-psys", "r")) {
    if (fscanf(fff, "event=%x", &this->_psys_config) != 1) {
      AutoPasLog(warn, "psys measurement support detected, but failed to parse config file");
      _psys_config = 0;
    }
    fclose(fff);
  } else {
    _psys_fd = -1;
  }
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/events/energy-psys.scale", "r")) {
    if (fscanf(fff, "%lf", &this->_psys_unit) != 1) {
      AutoPasLog(warn,
                 "psys energy measurement support detected, but failed to parse scale file, using 1.0 as fallback");
      _psys_unit = 1.0;
    }
    AutoPasLog(debug, "psys scale={} J", _psys_unit);
    fclose(fff);
  }

  int pkg_config;
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/events/energy-pkg", "r")) {
    if (fscanf(fff, "event=%x", &this->_pkg_config) != 1) {
      AutoPasLog(warn, "pkg energy measurement support detected, but failed to parse config file");
      _pkg_config = 0;
    }
    fclose(fff);
  } else {
    _pkg_fd = -1;
  }
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/events/energy-pkg.scale", "r")) {
    if (fscanf(fff, "%lf", &this->_pkg_unit) != 1) {
      AutoPasLog(warn,
                 "pkg energy measurement support detected, but failed to parse scale file, using 1.0 as fallback");
      _pkg_unit = 1.0;
    }
    AutoPasLog(debug, "pkg scale={} J", _pkg_unit);
    fclose(fff);
  }

  int cores_config;
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/events/energy-cores", "r")) {
    if (fscanf(fff, "event=%x", &this->_cores_config) != 1) {
      AutoPasLog(warn, "cores energy measurement support detected, but failed to parse config file");
      _cores_config = 0;
    }
    fclose(fff);
  } else {
    _cores_fd = -1;
  }
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/events/energy-cores.scale", "r")) {
    if (fscanf(fff, "%lf", &this->_cores_unit) != 1) {
      AutoPasLog(warn,
                 "cores energy measurement support detected, but failed to parse scale file, using 1.0 as fallback");
      _cores_unit = 1.0;
    }
    AutoPasLog(debug, "cores scale={} J", _cores_unit);
    fclose(fff);
  }

  int ram_config;
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/events/energy-ram", "r")) {
    if (fscanf(fff, "event=%x", &this->_ram_config) != 1) {
      AutoPasLog(warn, "ram energy measurement support detected, but failed to parse config file");
      _ram_config = 0;
    }
    fclose(fff);
  } else {
    _ram_fd = -1;
  }
  if (FILE *fff = fopen("/sys/bus/event_source/devices/power/events/energy-ram.scale", "r")) {
    if (fscanf(fff, "%lf", &this->_ram_unit) != 1) {
      AutoPasLog(warn,
                 "ram energy measurement support detected, but failed to parse scale file, using 1.0 as fallback");
      _ram_unit = 1.0;
    }
    AutoPasLog(debug, "ram scale={} J", _ram_unit);
    fclose(fff);
  }
}

RaplMeter::~RaplMeter() {
  if (_psys_fd != -1) {
    close(_psys_fd);
  }
  if (_pkg_fd != -1) {
    close(_pkg_fd);
  }
  if (_cores_fd != -1) {
    close(_cores_fd);
  }
  if (_ram_fd != -1) {
    close(_ram_fd);
  }
}

long RaplMeter::read_perf_event(int fd) {
  if (fd == -1) {
    return 0;
  }
  long value;
  lseek(fd, 0, SEEK_SET);
  if (read(fd, &value, 8) == -1) {
    throw ExceptionHandler::AutoPasException("Failed to read perf event:\n\t");
  }
  return value;
}

void RaplMeter::reset() {
  if (_psys_fd != -1) {
    close(_psys_fd);
  }
  if (_pkg_fd != -1) {
    close(_pkg_fd);
  }
  if (_cores_fd != -1) {
    close(_cores_fd);
  }
  if (_ram_fd != -1) {
    close(_ram_fd);
  }

  this->_psys_fd = open_perf_event(this->_type, this->_psys_config);
  this->_pkg_fd = open_perf_event(this->_type, this->_pkg_config);
  this->_cores_fd = open_perf_event(this->_type, this->_cores_config);
  this->_ram_fd = open_perf_event(this->_type, this->_ram_config);
}

void RaplMeter::sample() {
  this->_pkg_raw = read_perf_event(this->_pkg_fd);
  this->_cores_raw = read_perf_event(this->_cores_fd);
  this->_ram_raw = read_perf_event(this->_ram_fd);
  this->_psys_raw = read_perf_event(this->_psys_fd);
}

double RaplMeter::get_psys_energy() { return this->_psys_unit * this->_psys_raw; }
double RaplMeter::get_pkg_energy() { return this->_pkg_unit * this->_pkg_raw; }
double RaplMeter::get_cores_energy() { return this->_cores_unit * this->_cores_raw; }
double RaplMeter::get_ram_energy() { return this->_ram_unit * this->_ram_raw; }

long RaplMeter::get_total_energy() {
  if (this->_psys_fd != -1) {
    return this->_psys_raw;
  }
  long value = this->_pkg_raw;
  if (this->_ram_fd != -1) {
    value += this->_ram_raw;
  }
  return value;
}

}  // namespace autopas::utils

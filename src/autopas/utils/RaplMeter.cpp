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

int RaplMeter::open_perf_event(int type, int config, int cpu) {
  struct perf_event_attr attr;
  memset(&attr, 0, sizeof(attr));
  attr.type = type;
  attr.config = config;
  if (config == 0) {
    return -1;
  }

  int fd = syscall(__NR_perf_event_open, &attr, -1 /*PID*/, cpu, -1 /*GROUP FD*/, 0 /*FLAGS*/);
  if (fd < 0) {
    if (errno == EACCES) {
      throw ExceptionHandler::AutoPasException("Failed to open perf event: Permission denied");
    }
    throw ExceptionHandler::AutoPasException("Failed to open perf event");
  }
  return fd;
}

void RaplMeter::init() {
  {
    FILE *fd;
    int i = 0;
    int package;
    int last_package = 0;
    do {
      char filename[64];
      snprintf(filename, 63, "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", i);
      fd = fopen(filename, "r");
      if (!fd) {
        break;
      }
      if (fscanf(fd, "%d", &package) == 0) {
        throw ExceptionHandler::AutoPasException("failed to read cpu topology");
      }
      fclose(fd);
      if (this->_cpus.size() == 0 or last_package < package) {
        _cpus.push_back(i);
        last_package = package;
      }
      i++;
    } while (fd);
  }
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
  for (int f : _psys_fd) {
    if (f != -1) {
      close(f);
    }
  }
  for (int f : _pkg_fd) {
    if (f != -1) {
      close(f);
    }
  }
  for (int f : _cores_fd) {
    if (f != -1) {
      close(f);
    }
  }
  for (int f : _ram_fd) {
    if (f != -1) {
      close(f);
    }
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
  for (int f : this->_psys_fd) {
    if (f != -1) {
      close(f);
    }
  }
  for (int f : this->_pkg_fd) {
    if (f != -1) {
      close(f);
    }
  }
  for (int f : this->_cores_fd) {
    if (f != -1) {
      close(f);
    }
  }
  for (int f : this->_ram_fd) {
    if (f != -1) {
      close(f);
    }
  }

  this->_psys_fd.clear();
  this->_pkg_fd.clear();
  this->_cores_fd.clear();
  this->_ram_fd.clear();

  for (int i : this->_cpus) {
    this->_psys_fd.push_back(open_perf_event(this->_type, this->_psys_config, i));
    this->_pkg_fd.push_back(open_perf_event(this->_type, this->_pkg_config, i));
    this->_cores_fd.push_back(open_perf_event(this->_type, this->_cores_config, i));
    this->_ram_fd.push_back(open_perf_event(this->_type, this->_ram_config, i));
  }
}

void RaplMeter::sample() {
  this->_psys_raw = 0;
  for (int f : this->_psys_fd) {
    if (f != -1) {
      this->_psys_raw += read_perf_event(f);
    }
  }

  this->_pkg_raw = 0;
  for (int f : this->_pkg_fd) {
    if (f != -1) {
      this->_pkg_raw += read_perf_event(f);
    }
  }

  this->_cores_raw = 0;
  for (int f : this->_cores_fd) {
    if (f != -1) {
      this->_cores_raw += read_perf_event(f);
    }
  }

  this->_ram_raw = 0;
  for (int f : this->_ram_fd) {
    if (f != -1) {
      this->_ram_raw += read_perf_event(f);
    }
  }
}

double RaplMeter::get_psys_energy() { return this->_psys_unit * this->_psys_raw; }
double RaplMeter::get_pkg_energy() { return this->_pkg_unit * this->_pkg_raw; }
double RaplMeter::get_cores_energy() { return this->_cores_unit * this->_cores_raw; }
double RaplMeter::get_ram_energy() { return this->_ram_unit * this->_ram_raw; }

long RaplMeter::get_total_energy() {
  if (this->_psys_raw > 0) {
    return this->_psys_raw;
  }
  return this->_pkg_raw + this->_ram_raw;
}

}  // namespace autopas::utils

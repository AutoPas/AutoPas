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
  // TODO: Error handling
  return fd;
}

RaplMeter::RaplMeter() {
  FILE *fff = fopen("/sys/bus/event_source/devices/power/type", "r");
  if (fff == NULL) {
    throw ExceptionHandler::AutoPasException("No support for energy measurements detected.");
  } else {
    fscanf(fff, "%d", &this->_type);
    fclose(fff);
  }

  int psys_config;
  fff = fopen("/sys/bus/event_source/devices/power/events/energy-psys", "r");
  if (fff != NULL) {
    fscanf(fff, "event=%x", &this->_psys_config);
    fclose(fff);
  } else {
    _psys_fd = -1;
  }
  fff = fopen("/sys/bus/event_source/devices/power/events/energy-psys.scale", "r");
  if (fff != NULL) {
    fscanf(fff, "%lf", &this->_psys_unit);
    AutoPasLog(debug, "psys scale={} J", _psys_unit);
    fclose(fff);
  }

  int pkg_config;
  fff = fopen("/sys/bus/event_source/devices/power/events/energy-pkg", "r");
  if (fff != NULL) {
    fscanf(fff, "event=%x", &this->_pkg_config);
    fclose(fff);
  } else {
    _pkg_fd = -1;
  }
  fff = fopen("/sys/bus/event_source/devices/power/events/energy-pkg.scale", "r");
  if (fff != NULL) {
    fscanf(fff, "%lf", &this->_pkg_unit);
    AutoPasLog(debug, "pkg scale={} J", _pkg_unit);
    fclose(fff);
  }

  int cores_config;
  fff = fopen("/sys/bus/event_source/devices/power/events/energy-cores", "r");
  if (fff != NULL) {
    fscanf(fff, "event=%x", &this->_cores_config);
    fclose(fff);
  } else {
    _cores_fd = -1;
  }
  fff = fopen("/sys/bus/event_source/devices/power/events/energy-cores.scale", "r");
  if (fff != NULL) {
    fscanf(fff, "%lf", &this->_cores_unit);
    AutoPasLog(debug, "cores scale={} J", _cores_unit);
    fclose(fff);
  }

  int ram_config;
  fff = fopen("/sys/bus/event_source/devices/power/events/energy-ram", "r");
  if (fff != NULL) {
    fscanf(fff, "event=%x", &this->_ram_config);
    fclose(fff);
  } else {
    _ram_fd = -1;
  }
  fff = fopen("/sys/bus/event_source/devices/power/events/energy-ram.scale", "r");
  if (fff != NULL) {
    fscanf(fff, "%lf", &this->_ram_unit);
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
  read(fd, &value, 8);
  return value;
}

void RaplMeter::reset() {
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

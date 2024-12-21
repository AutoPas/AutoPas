/**
 * @file AccumulatorMeasurements.cpp
 * @author L. Llaveshi
 * @date 05.12.2024
 */

#include "AccumulatorMeasurements.h"

#include <sstream>

namespace autopas {

AccumulatorMeasurements::AccumulatorMeasurements()
    : accumulatedEnergyPsys(0.0), accumulatedEnergyPkg(0.0), accumulatedEnergyRam(0.0) {}

void AccumulatorMeasurements::addEnergy(double psys, double pkg, double ram) {
  accumulatedEnergyPsys += psys;
  accumulatedEnergyPkg += pkg;
  accumulatedEnergyRam += ram;
}

double AccumulatorMeasurements::getAccumulatedEnergyPsys() const { return accumulatedEnergyPsys; }

double AccumulatorMeasurements::getAccumulatedEnergyPkg() const { return accumulatedEnergyPkg; }

double AccumulatorMeasurements::getAccumulatedEnergyRam() const { return accumulatedEnergyRam; }

std::string AccumulatorMeasurements::toString() const {
  std::ostringstream oss;
  oss << "Total Accumulated Energy: [Psys = " << accumulatedEnergyPsys
      << ", Pkg = " << accumulatedEnergyPkg
      << ", Ram = " << accumulatedEnergyRam << "]";
  return oss.str();
}

}  // namespace autopas

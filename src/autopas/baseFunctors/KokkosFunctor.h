/**
* @file KokkosFunctor.h
 *
 * @date 25.06.2026
 * @author Luis Gall
 */

#pragma once

namespace autopas {

class KokkosFunctor {

public:

  void setPotentialEnergy(double energy) {
    _potentialEnergy = energy;
  }

  void setVirial(double virial) {
    _virial = virial;
  }

protected:

  double _potentialEnergy = 0.;

  double _virial = 0.;

};

}


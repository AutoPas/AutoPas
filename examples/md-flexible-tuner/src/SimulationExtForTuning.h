#pragma once

#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>

#include "src/TimeDiscretization.h"
#include "autopas/AutoPasDecl.h"
#include "src/ParallelVtkWriter.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/domainDecomposition/DomainDecomposition.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

#include "src/Simulation.h"

class SimulationExtForTuning : public Simulation {
 public:
    SimulationExtForTuning(const MDFlexConfig &configuration, std::shared_ptr<RegularGridDecomposition> &domainDecomposition) : Simulation(configuration, domainDecomposition) {};

    void processCell(int numParticles, size_t sortingThreshold);
    void processCellPair(int numParticles, size_t sortingThreshold);

  template <class T, class F>
  T applyWithChosenFunctor(F f);
};

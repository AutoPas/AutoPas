#include "SimulationExtForTuning.h"

#include "src/TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "autopas/utils/SimilarityFunctions.h"
#include "autopas/utils/WrapMPI.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/cells/ParticleCell.h"

// Declare the main AutoPas class and the iteratePairwise() methods with all used functors as extern template
// instantiation. They are instantiated in the respective cpp file inside the templateInstantiations folder.
//! @cond Doxygen_Suppress
extern template class autopas::AutoPas<ParticleType>;
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
typedef LJFunctorTypeAutove FunctorType;
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
typedef LJFunctorTypeAutovecGlobals FunctorType;
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
typedef LJFunctorTypeAVX FunctorType;
#endif
#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
typedef LJFunctorTypeSVE FunctorType;
#endif
//! @endcond

bool callProcessCell(FunctorType *f , autopas::FullParticleCell<ParticleType> &cell,
                        double interactionLength, autopas::DataLayoutOption dataLayout, bool useNewton3, size_t sortingThreshold) {
    autopas::internal::CellFunctor<autopas::FullParticleCell<ParticleType>, FunctorType, false> cf(f, interactionLength, dataLayout, useNewton3);
    cf.setSortingThreshold(sortingThreshold);
    cf.processCell(cell);
    return true;
}

bool callProcessCellPair(FunctorType *f , autopas::FullParticleCell<ParticleType> &cell1, 
                        autopas::FullParticleCell<ParticleType> &cell2, const std::array<double, 3> &sortingDirection,
                        double interactionLength, autopas::DataLayoutOption dataLayout, bool useNewton3, size_t sortingThreshold) {
    autopas::internal::CellFunctor<autopas::FullParticleCell<ParticleType>, FunctorType, false> cf(f, interactionLength, dataLayout, useNewton3);
    cf.setSortingThreshold(sortingThreshold);
    cf.processCellPair(cell1, cell2, sortingDirection);
    return true;
}

void SimulationExtForTuning::processCell(int numParticles, size_t sortingThreshold) {
    autopas::FullParticleCell<ParticleType> cell;
    cell.reserve(numParticles);

    size_t particleId = 0;
    unsigned long _typeId = 0;
    std::array<double, 3> velocity = {0, 0, 0};
    ParticleType particle{};
    particle.setID(particleId);
    particle.setTypeId(_typeId);
    particle.setOwnershipState(autopas::OwnershipState::owned);
    particle.setV(velocity);
    particle.setF({0.0, 0.0, 0.0});
    particle.setOldF({0.0, 0.0, 0.0});

    //particles are randomly distributed in  a cell
    for (int ii = 0; ii < numParticles; ii++) {
          double fx = static_cast<double>(rand()) / RAND_MAX * cell.getCellLength()[0];
          double fy = static_cast<double>(rand()) / RAND_MAX * cell.getCellLength()[1];
          double fz = static_cast<double>(rand()) / RAND_MAX * cell.getCellLength()[2];
          particle.setR({fx, fy, fz });
          cell.addParticle(particle);
          particle.setID(particle.getID() + 1);
    }

    const double interactionLength = _configuration.cutoff.value + _configuration.verletSkinRadiusPerTimestep.value * _configuration.verletRebuildFrequency.value;
    auto dataLayout = _autoPasContainer->getCurrentConfig().dataLayout; // should be DataLayoutOption::aos but is not
    dataLayout = autopas::DataLayoutOption::aos;
    const auto useNewton3 = _autoPasContainer->getCurrentConfig().newton3;

    applyWithChosenFunctor<bool>([&](auto functor) { return callProcessCell(&functor, cell, interactionLength, dataLayout, useNewton3, sortingThreshold); });
}


void SimulationExtForTuning::processCellPair(int numParticles, size_t sortingThreshold) {
    autopas::FullParticleCell<ParticleType> cell1;
    autopas::FullParticleCell<ParticleType> cell2;
    cell1.reserve(numParticles/2);
    cell2.reserve(numParticles-(numParticles/2));

    size_t particleId = 0;
    unsigned long _typeId = 0;
    std::array<double, 3> velocity = {0, 0, 0};
    ParticleType particle{};
    particle.setID(particleId);
    particle.setTypeId(_typeId);
    particle.setOwnershipState(autopas::OwnershipState::owned);
    particle.setV(velocity);
    particle.setF({0.0, 0.0, 0.0});
    particle.setOldF({0.0, 0.0, 0.0});

    //particles are randomly distributed in cell 1
    for (int ii = 0; ii < numParticles/2; ii++) {
          double fx = static_cast<double>(rand()) / RAND_MAX * cell1.getCellLength()[0];
          double fy = static_cast<double>(rand()) / RAND_MAX * cell1.getCellLength()[1];
          double fz = static_cast<double>(rand()) / RAND_MAX * cell1.getCellLength()[2];
          particle.setR({fx, fy, fz });
          cell1.addParticle(particle);
          particle.setID(particle.getID() + 1);
    }

    //particles are randomly distributed in (neighboring) cell 2
    for (int jj = numParticles/2; jj < numParticles; jj++) {
          double fx = static_cast<double>(rand()) / RAND_MAX * cell2.getCellLength()[0] + 1.0;
          double fy = static_cast<double>(rand()) / RAND_MAX * cell2.getCellLength()[1];
          double fz = static_cast<double>(rand()) / RAND_MAX * cell2.getCellLength()[2];
          particle.setR({fx, fy, fz });
          cell2.addParticle(particle);
          particle.setID(particle.getID() + 1);
    }

    const double interactionLength = _configuration.cutoff.value + _configuration.verletSkinRadiusPerTimestep.value * _configuration.verletRebuildFrequency.value;
    auto dataLayout = _autoPasContainer->getCurrentConfig().dataLayout; // should be DataLayoutOption::aos but is not
    dataLayout = autopas::DataLayoutOption::aos;
    const auto useNewton3 = _autoPasContainer->getCurrentConfig().newton3;

    const std::array<double, 3> &sortingDirection = {1.0, 0.0, 0.0};
    applyWithChosenFunctor<bool>([&](auto functor) { return callProcessCellPair(&functor, cell1, cell2, sortingDirection, interactionLength, dataLayout, useNewton3, sortingThreshold); });
}

template <class T, class F>
T SimulationExtForTuning::applyWithChosenFunctor(F f) {
  const double cutoff = _configuration.cutoff.value;
  auto &particlePropertiesLibrary = *_configuration.getParticlePropertiesLibrary();
  switch (_configuration.functorOption.value) {
    case MDFlexConfig::FunctorOption::lj12_6: {
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC)
      return f(LJFunctorTypeAutovec{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor AutoVec. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AUTOVEC=ON`.");
#endif
    }
    case MDFlexConfig::FunctorOption::lj12_6_Globals: {
#if defined(MD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS)
      return f(LJFunctorTypeAutovecGlobals{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor AutoVec Globals. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AUTOVEC_GLOBALS=ON`.");
#endif
    }
    case MDFlexConfig::FunctorOption::lj12_6_AVX: {
#if defined(MD_FLEXIBLE_FUNCTOR_AVX) && defined(__AVX__)
      return f(LJFunctorTypeAVX{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor AVX. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_AVX=ON`.");
#endif
    }
    case MDFlexConfig::FunctorOption::lj12_6_SVE: {
#if defined(MD_FLEXIBLE_FUNCTOR_SVE) && defined(__ARM_FEATURE_SVE)
      return f(LJFunctorTypeSVE{cutoff, particlePropertiesLibrary});
#else
      throw std::runtime_error(
          "MD-Flexible was not compiled with support for LJFunctor SVE. Activate it via `cmake "
          "-DMD_FLEXIBLE_FUNCTOR_SVE=ON`.");
#endif
    }
  }
  throw std::runtime_error("Unknown functor choice" +
                           std::to_string(static_cast<int>(_configuration.functorOption.value)));
}

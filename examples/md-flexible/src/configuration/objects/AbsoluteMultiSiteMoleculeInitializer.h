//
// Created by johnny on 14.10.23.
//

#pragma once
#include "Object.h"
//#include "../../applicationLibrary/molecularDynamics/molecularDynamicsLibrary/AbsoluteMultiSiteMoleculeLJ.h"

#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
namespace AbsoluteMultiSiteMoleculeInitializer {
/**
 * Finishes the Initialization of MutliSite Particles using the AbsoluteMultiSiteMoleculeLJ class as representation.
 * The problem is that by the time of construction of the particle the absolute position
 * cannot be computed since information laying in the particlePropertiesLibrary is missing.
 * Therefore the Molecules get loaded, then the ppl gets loaded and finally this method finishes the initialization
 * of the molecules
 * @param p Particle (must be of type AbsoluteMultiSiteMoleculeLJ) that gets its absolute Site Positions initialized
 * @param ppl particlePropertiesLibrary
 */
void setAbsoluteSites(ParticleType &p, const std::shared_ptr<const ParticlePropertiesLibraryType> ppl) {
  // determine and initialize absolute Site positions
  std::vector<std::array<double, 3>> relativeSitePositions = ppl->getSitePositions(p.getTypeId());
  std::vector<double> absPositionsX{};
  std::vector<double> absPositionsY{};
  std::vector<double> absPositionsZ{};
  for (int i = 0; i < relativeSitePositions.size(); i++) {
    std::array<double, 3> relativeSitePosition = relativeSitePositions[i];
    absPositionsX.push_back(relativeSitePosition[0] + p.getR()[0]);
    absPositionsY.push_back(relativeSitePosition[1] + p.getR()[1]);
    absPositionsZ.push_back(relativeSitePosition[2] + p.getR()[2]);
  }
  p.setAbsoluteSitePositionsX(absPositionsX);
  p.setAbsoluteSitePositionsY(absPositionsY);
  p.setAbsoluteSitePositionsZ(absPositionsZ);
}
}
#endif
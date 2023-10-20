//
// Created by johnny on 20.10.23.
//

#include "AbsoluteMultiSiteMoleculeInitializer.h"

#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
namespace AbsoluteMultiSiteMoleculeInitializer {
void setAbsoluteSites(ParticleType &p, const std::shared_ptr<const ParticlePropertiesLibraryType> ppl) {
  const ParticlePropertiesLibraryType * ppl_ptr = ppl.get();  //black magic to get const reference out of const std::shared_ptr
  p.setQuaternion(p.getQuaternion(), *ppl_ptr);

  //p.setQuaternion(p.getQuaternion(), ppl->getSitePositions(p.getTypeId()));

  //// determine and initialize absolute Site positions
  //std::vector<std::array<double, 3>> unrotated_site_positions = ppl->getSitePositions(p.getTypeId());
  //const auto q = p.getQuaternion();
  //std::vector<std::array<double, 3>> rotated_site_positions = autopas::utils::quaternion::rotateVectorOfPositions(q, unrotated_site_positions);
  //std::vector<double> absPositionsX{};
  //std::vector<double> absPositionsY{};
  //std::vector<double> absPositionsZ{};
  //for (int i = 0; i < rotated_site_positions.size(); i++) {
  //  std::array<double, 3> relativeSitePosition = rotated_site_positions[i];
  //  absPositionsX.push_back(relativeSitePosition[0] + p.getR()[0]);
  //  absPositionsY.push_back(relativeSitePosition[1] + p.getR()[1]);
  //  absPositionsZ.push_back(relativeSitePosition[2] + p.getR()[2]);
  //}
  //p.setAbsoluteSitePositionsX(absPositionsX);
  //p.setAbsoluteSitePositionsY(absPositionsY);
  //p.setAbsoluteSitePositionsZ(absPositionsZ);
}
}

#endif
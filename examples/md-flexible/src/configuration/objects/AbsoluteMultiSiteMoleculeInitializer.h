//
// Created by johnny on 14.10.23.
//

#pragma once
#include "Object.h"
#include "autopas/utils/Quaternion.h"
// #include "../../applicationLibrary/molecularDynamics/molecularDynamicsLibrary/PositionStoringMultiSiteMolecule.h"

#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
namespace AbsoluteMultiSiteMoleculeInitializer {
/**
 * Finishes the Initialization of MutliSite Particles using the PositionStoringMultiSiteMolecule class as representation.
 * The problem is that by the time of construction of the particle the absolute position
 * cannot be computed since information laying in the particlePropertiesLibrary is missing.
 * Therefore the Molecules get loaded, then the ppl gets loaded and finally this method finishes the initialization
 * of the molecules
 * @param p Particle (must be of type PositionStoringMultiSiteMolecule) that gets its absolute Site Positions initialized
 * @param ppl particlePropertiesLibrary
 */
void setAbsoluteSites(ParticleType &p, const std::shared_ptr<const ParticlePropertiesLibraryType> ppl);
}
#endif
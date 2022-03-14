/**
 * @file SimpleMDAutoPas.h
 * @date 11/03/2022
 * @author S. Newcome
*/

#pragma once

#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include "Particles/MulticenteredMoleculeLJ.h"

class SimpleMDAutoPas : autopas::AutoPas<autopas::MoleculeLJ> {};

class ComplexMDAutoPas : autopas::AutoPas<MulticenteredMoleculeLJ>{};
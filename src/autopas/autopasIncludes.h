/**
 * @file autopasIncludes.h
 * File to handle the includes of autopas. This file should not be included from the outside unless you want to use
 * the lower level routines directly, without the use of the main AutoPas class.
 * In the normal case, please include AutoPas.h.
 * @author tchipevn
 * @date 17.01.2018
 */

#pragma once

/// @todo separate autopas.h and autopasmd.h

// utils
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Logger.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticSelectorMacros.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"

// particles
#include "autopas/particles/MoleculeLJ.h"
#include "autopas/particles/Particle.h"

// cells
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/cells/RMMParticleCell2T.h"
#include "autopas/containers/DirectSum.h"

// iterators
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/iterators/SingleCellIterator.h"

// traversals
#include "autopas/containers/cellPairTraversals/C01Traversal.h"
#include "autopas/containers/cellPairTraversals/C08Traversal.h"
#include "autopas/containers/cellPairTraversals/C18Traversal.h"
#include "autopas/containers/cellPairTraversals/SlicedTraversal.h"

// containers
#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/LinkedCells.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/VerletClusterLists.h"
#include "autopas/containers/VerletLists.h"
#include "autopas/containers/VerletListsCells.h"

// pairwise functors
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"

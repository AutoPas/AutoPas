#pragma once
#include <getopt.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include "autopas/utils/NumberSet.h"

//enspricht genau #include <autopas.h> bloß ohne LJFunctor.h um multiple declarations zu vermeiden
// utils
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Logger.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticSelectors.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"

// particles
#include "autopas/particles/MoleculeLJ.h"
#include "autopas/particles/Particle.h"

// cells
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/cells/RMMParticleCell2T.h"
#include "autopas/containers/directSum/DirectSum.h"

// iterators
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/iterators/SingleCellIterator.h"

// traversals
#include "autopas/containers/directSum/DirectSumTraversal.h"
#include "autopas/containers/linkedCells/traversals/C01Traversal.h"
#include "autopas/containers/linkedCells/traversals/C08Traversal.h"
#include "autopas/containers/linkedCells/traversals/C18Traversal.h"
#include "autopas/containers/linkedCells/traversals/SlicedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/C01TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/C18TraversalVerlet.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/SlicedTraversalVerlet.h"

// containers
#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"

// pairwise functors
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"
#include "autopas/pairwiseFunctors/Functor.h"

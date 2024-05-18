/**
 * @file LCTriTraversalInterface.h
 * @author muehlhaeusser
 * @date 18.05.2024
 */

#pragma once

#include <vector>

#include "LCTraversalInterface.h"
#include "autopas/containers/TriwiseTraversalInterface.h"

namespace autopas {

/**
 * Interface for traversals used by the LinkedCell class.
 *
 * The container only accepts traversals in its iterateTriwise() method that implement this interface.
 */
template <class ParticleCell>
class LCTriTraversalInterface : public TriwiseTraversalInterface, public LCTraversalInterface {};
}  // namespace autopas
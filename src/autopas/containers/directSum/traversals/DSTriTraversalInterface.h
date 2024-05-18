/**
 * @file DSTriTraversalInterface.h
 * @author muehlhaeusser
 * @date 18.05.24
 */

#pragma once

#include <vector>

#include "DSTraversalInterface.h"
#include "autopas/containers/TriwiseTraversalInterface.h"

namespace autopas {

/**
 * Interface for traversals used by the DirectSum container.
 *
 * The container only accepts traversals in its iterateTriwise() method that implement this interface.
 */
class DSTriTraversalInterface : public DSTraversalInterface, public TriwiseTraversalInterface {};
}  // namespace autopas
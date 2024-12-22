
#pragma once
#include "src/zonalMethods/region/RectRegion.h"

/**
 * This class defines a common interface for all zonal methods
 * that use RectRegion for their zones / regions.
 */
class RectRegionMethodInterface {
 public:
  virtual const std::vector<RectRegion> getExportRegions() { return {}; }

  virtual const std::vector<RectRegion> getImportRegions() { return {}; }
};

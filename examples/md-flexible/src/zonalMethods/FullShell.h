
#include "src/zonalMethods/ZonalMethod.h"
#include "src/zonalMethods/region/RectRegion.h"

class FullShell : public ZonalMethod {
 public:
  /**
   * Constructor
   * @param cutoff
   * @param verletSkinWidth
   * */
  FullShell(RectRegion homeBoxRegion, double cutoff, double verletSkinWidth);

  /**
   * Destructor
   * */
  ~FullShell();

  /**
   * Collect particles from the AutoPas container and store them internally.
   * @param autoPasContainer
   */
  void collectParticles(AutoPasType &autoPasContainer) override;

 protected:
  /**
   * The number of export regions
   * */
  inline static const size_t _regionCount = 26;

  /**
   * Stores the export regions of this box
   * */
  std::vector<RectRegion> _exportRegions;

  /**
   * Stores the import regions of this box
   * */
  std::vector<RectRegion> _importRegions;

  /**
   * Stores the particles to be exported in each region
   * */
  std::array<std::vector<ParticleType>, _regionCount> _regionBuffers;
};

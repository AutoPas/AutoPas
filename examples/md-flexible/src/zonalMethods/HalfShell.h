
#include "src/zonalMethods/ZonalMethod.h"
#include "src/zonalMethods/region/RectRegion.h"

/**
 * Class for the HalfShell ZonalMethod.
 * When initialized, it calculates the export and import regions of the 
 * given home box.
 */
class HalfShell : public ZonalMethod {
 public:

  /**
   * Constructor
   * @param cutoff
   * @param verletSkinWidth
   * */
  HalfShell(RectRegion homeBoxRegion, double cutoff, double verletSkinWidth);


  /**
   * Destructor
   * */
  ~HalfShell();

  /**
   * Collect particles from the AutoPas container and store them internally.
   * @param autoPasContainer
   */
  void collectParticles(AutoPasType &autoPasContainer) override;

 protected:

  /**
   * The number of export regions
   * */
  inline static const size_t _regionCount = 13;

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

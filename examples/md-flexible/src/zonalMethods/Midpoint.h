

#pragma once
#include "src/options/BoundaryTypeOption.h"
#include "src/zonalMethods/ZonalMethod.h"
#include "src/zonalMethods/region/RectRegion.h"
#include "src/zonalMethods/RectRegionMethodInterface.h"

/**
 * Class for the Midpoint ZonalMethod.
 * When initialized, it calculates the export and import regions of the
 * given home box.
 */
class Midpoint : public ZonalMethod, public RectRegionMethodInterface {
 public:
  /**
   * Constructor
   * @param cutoff
   * @param verletSkinWidth
   * @param ownRank
   * @param homeBoxRegion
   * @param globalBoxRegion
   * @param comm (optional)
   * @param allNeighbourIndices (optional)
   * @param boundaryType (optional)
   * */
  Midpoint(double cutoff, double verletSkinWidth, int ownRank, RectRegion homeBoxRegion, RectRegion globalBoxRegion,
           autopas::AutoPas_MPI_Comm comm = AUTOPAS_MPI_COMM_WORLD,
           std::array<int, 26> allNeighbourIndices = std::array<int, 26>(),
           std::array<options::BoundaryTypeOption, 3> boundaryType = std::array<options::BoundaryTypeOption, 3>(
               {options::BoundaryTypeOption::periodic, options::BoundaryTypeOption::periodic,
                options::BoundaryTypeOption::periodic}));

  /**
   * Destructor
   * */
  ~Midpoint();

  /**
   * Collect particles from the AutoPas container and store them internally.
   * @param autoPasContainer
   */
  void collectParticles(AutoPasType &autoPasContainer) override;

  /**
   * Send and receive exports.
   * Received particles are stored internally.
   * @param autoPasContainer
   * @param comm
   * @param allNeighbourIndices
   * @param ownRank
   * @param boundaryType
   */
  void SendAndReceiveExports(AutoPasType &autoPasContainer) override;
  /**
   * Send and receive results of the force calculation and
   * store them into the respective particles in the AutoPas container.
   * This function should be called after force calculation, before
   * the velocity calculation.
   * @param autoPasContainer
   * @param comm
   * @param allNeighbourIndices
   */
  void SendAndReceiveResults(AutoPasType &autoPasContainer) override;

  /**
   * Recollect the halo particles from the AutoPas container and
   * save them internally.
   * As this is method specific, this function is pure virtual.
   * This is called before calculateExternalZonalInteractions().
   * @param autoPasContainer
   */
  void recollectResultsFromContainer(AutoPasType &autoPasContainer) override;

  /**
   * Get the export regions
   * @return
   */
  const std::vector<RectRegion> getExportRegions() override;

  /**
   * Get the import regions
   * @return
   */
  const std::vector<RectRegion> getImportRegions() override;

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

  /*
   * Stores the particles imported from neighbours.
   */
  std::array<std::vector<ParticleType>, _regionCount> _importBuffers;

  void calculateZonalInteractionPairwise(std::string zone1, std::string zone2,
                                         std::function<void(ParticleType &, ParticleType &)> aosFunctor) override;

  /**
   * Calculates the _interactionZones and _interactionSchedule
   * @param identifyZone
   */
  void calculateInteractionSchedule(std::function<std::string(const int[3])> identifyZone);
};

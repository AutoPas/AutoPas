
#pragma once
#include <map>

#include "autopas/AutoPas.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/options/InteractionTypeOption.h"
#include "src/TypeDefinitions.h"
#include "src/options/BoundaryTypeOption.h"
#include "src/zonalMethods/region/RectRegion.h"

/**
 * This abstract class represents a zonal method.
 * All specific zonal method classes should extend this class.
 * */
class ZonalMethod {
 public:
  /**
   * Type for the AutoPas container
   */
  using AutoPasType = autopas::AutoPas<ParticleType>;

  /**
   * Constructor
   * @param zoneCount
   * @param ownRank
   * @param homeBoxRegion
   * @param globalBoxRegion
   * @param comm (optional)
   * @param allNeighbourIndices (optional)
   * @param boundaryType (optional)
   */
  ZonalMethod(unsigned int zoneCount, int ownRank, RectRegion homeBoxRegion, RectRegion globalBoxRegion,
              autopas::AutoPas_MPI_Comm comm = AUTOPAS_MPI_COMM_WORLD,
              std::array<int, 26> allNeighbourIndices = std::array<int, 26>(),
              std::array<options::BoundaryTypeOption, 3> boundaryType = std::array<options::BoundaryTypeOption, 3>());

  /**
   * Destructor
   */
  virtual ~ZonalMethod();

  /**
   * Collect particles from the AutoPas container and store them internally.
   * @param autoPasContainer
   */
  virtual void collectParticles(AutoPasType &autoPasContainer) = 0;

  /**
   * Send and receive exports.
   * Received particles are stored internally.
   * @param autoPasContainer
   */
  virtual void SendAndReceiveExports(AutoPasType &autoPasContainer) = 0;

  /**
   * Send and receive results of the force calculation and
   * store them into the respective particles in the AutoPas container.
   * This function should be called after force calculation, before
   * the velocity calculation.
   * @param autoPasContainer
   */
  virtual void SendAndReceiveResults(AutoPasType &autoPasContainer) = 0;

  /**
   * Calculate the interaction between the external zones.
   * This function should be called after recollectResultsFromContainer().
   * @param autoPasContainer
   * @param particleProperties
   * @param cutoff
   */
  virtual void calculateExternalZonalInteractions(AutoPasType &autoPasContainer,
                                                  std::shared_ptr<ParticlePropertiesLibraryType> particleProperties,
                                                  double cutoff,
                                                  std::set<autopas::InteractionTypeOption> interactionTypes);

  /**
   * Recollect the halo particles from the AutoPas container and
   * save them internally.
   * As this is method specific, this function is pure virtual.
   * This is called before calculateExternalZonalInteractions().
   * @param autoPasContainer
   */
  virtual void recollectResultsFromContainer(AutoPasType &autoPasContainer) = 0;


  /**
   * Resizes the home box region and thus recalculates
   * the import and export regions
   * @param homeBoxRegion
   */
  virtual void resizeHomeBoxRegion(RectRegion homeBoxRegion) = 0;

 protected:
  /**
   * Stores the number of zones (except for the home zone)
   */
  unsigned int _zoneCount;

  /**
   * Stores the rank of the home box
   */
  int _ownRank;

  /**
   * Stores the home box region
   */
  RectRegion _homeBoxRegion;

  /**
   * Stores the global box region
   */
  RectRegion _globalBoxRegion;

  /**
   * Stores the MPI communicator
   */
  autopas::AutoPas_MPI_Comm _comm;

  /**
   * Stores the ranks of the neighbour processes
   */
  std::array<int, 26> _allNeighbourIndices;

  /**
   * Stores the boundary types of the global box
   */
  std::array<options::BoundaryTypeOption, 3> _boundaryType;

  /**
   * Stores whether the home box is at the minimal global boundary:
   *  1 - yes
   *  0 - no
   */
  std::array<int, 3> _isAtGlobalBoundaryMin;

  /**
   * Stores whether the home box is at the maximal global boundary
   *  1 - yes
   *  0 - no
   */
  std::array<int, 3> _isAtGlobalBoundaryMax;

  /**
   * Stores the interaction schedule without the home box zone
   */
  std::vector<std::string> _interactionZones;

  /**
   * Stores the interaction zones
   */
  std::map<std::string, std::vector<std::string>> _interactionSchedule;

  /**
   * Calculates import / export regions as RectRegion classes, saving them in the given buffer.
   * This function only considers the 26 neighbours of the home box, as it is assumed that the
   * cutoff radius + verlet skin width is smaller than the box sizes.
   * @param homeBoxRegion: region of the home bex
   * @param cutoffRadius: the cutoff radius
   * @param verletSkinWidth: the verlet skin width
   * @param regions: the export regions (normalized)
   * @param condition: the condition specifying which neighbours to consider
   * @param identifyZone: the zone identification function
   * @param calcImports: whether to calculate import regions (default: true)
   */
  virtual void getRectRegionsConditional(RectRegion &homeBoxRegion, double cutoffRadius, double verletSkinWidth,
                                         std::vector<RectRegion> &regions,
                                         const std::function<bool(const int[3])> &condition,
                                         const std::function<std::string(const int[3])> &identifyZone,
                                         bool calcImports = true);

  /**
   * Calculate the index (for RegularGridDecomposition:_allNeighbours)
   * of the given relative neighbour
   * @param relNeighbour
   */
  virtual size_t convRelNeighboursToIndex(std::array<int, 3> relNeighbour);

  /**
   * Check if there is a need to collect particles from the given neighbour
   * @param relNeighbour
   */
  virtual bool needToCollectParticles(std::array<int, 3> relNeighbour);

  /**
   * If the neighbour is over a global periodic boundary, wrap the particles positions
   * around the global box.
   * @param relNeighbour
   * @param particles
   */
  virtual void wrapAroundPeriodicBoundary(std::array<int, 3> relNeighbour, std::vector<ParticleType> &particles);

  /**
   * Calculate the zonal pairwise interaction between the given zones and the given AoS functor.
   * This is method dependant and should be implemented by each respective subclass.
   * The function is called by calculateExternalZonalInteractions().
   * @param zone1
   * @param zone2
   * @param aosFunctor
   */
  virtual void calculateZonalInteractionPairwise(std::string zone1, std::string zone2,
                                                 std::function<void(ParticleType &, ParticleType &)> aosFunctor) = 0;

  /*
   * Calculate the zonal triwise interactions of a given zone and a AoS functor
   * This is method dependant and should be implemented by each respective subclass.
   * The function is called by calculateExternalZonalInteractions().
   * @param zone
   * @param aosFunctor
   */
  virtual void calculateZonalInteractionTriwise(
      std::string zone, std::function<void(ParticleType &, ParticleType &, ParticleType &, bool)> aosFunctor) = 0;
};

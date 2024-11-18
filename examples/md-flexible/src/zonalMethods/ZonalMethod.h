
#pragma once
#include "autopas/AutoPas.h"
#include "autopas/AutoPasDecl.h"
#include "src/TypeDefinitions.h"
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
   */
  ZonalMethod(unsigned int zoneCount);

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
   * @param comm
   * @param allNeighbourIndices
   */
  virtual void SendAndReceiveExports(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm,
                                     std::array<int, 26> allNeighbourIndices) = 0;

  /**
   * Send and receive results of the force calculation and
   * store them into the respective particles in the AutoPas container.
   * This function should be called after force calculation, before
   * the velocity calculation.
   * @param autoPasContainer
   * @param comm
   * @param allNeighbourIndices
   */
  virtual void SendAndReceiveResults(AutoPasType &autoPasContainer, autopas::AutoPas_MPI_Comm,
                                     std::array<int, 26> allNeighbourIndices) = 0;

 protected:
  /**
   * Stores the number of zones (except for the home zone)
   */
  unsigned int _zoneCount;

  /**
   * Calculates import / export regions as RectRegion classes, saving them in the given buffer.
   * This function only considers the 26 neighbours of the home box, as it is assumed that the
   * cutoff radius + verlet skin width is smaller than the box sizes.
   * @param homeBoxRegion: region of the home bex
   * @param cutoffRadius: the cutoff radius
   * @param verletSkinWidth: the verlet skin width
   * @param regions: the export regions
   * @param condition: the condition specifying which neighbours to consider
   * @param calcImports: whether to calculate import regions (default: true)
   */
  virtual void getRectRegionsConditional(RectRegion &homeBoxRegion, double cutoffRadius, double verletSkinWidth,
                                         std::vector<RectRegion> &regions,
                                         const std::function<bool(const int[3])> &condition, bool calcImports = true);

  /**
   * Calculate the index (for RegularGridDecomposition:_allNeighbours)
   * of the given relative neighbour
   * @param relNeighbour
   */
  virtual size_t convRelNeighboursToIndex(std::array<int, 3> relNeighbour);
};

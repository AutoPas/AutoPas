

#pragma once
#include "src/options/BoundaryTypeOption.h"
#include "src/zonalMethods/RectRegionMethodInterface.h"
#include "src/zonalMethods/ZonalMethod.h"
#include "src/zonalMethods/region/RectRegion.h"

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
   * @param useNewton3
   * @param pairwiseInteraction
   * @param comm (optional)
   * @param allNeighbourIndices (optional)
   * @param boundaryType (optional)
   * */
  Midpoint(double cutoff, double verletSkinWidth, int ownRank, RectRegion homeBoxRegion, RectRegion globalBoxRegion,
           bool useNewton3 = false, bool pairwiseInteraction = true,
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

  void calculateZonalInteractionTriwise(
      std::string zone, std::function<void(ParticleType &, ParticleType &, ParticleType &, bool)> aosFunctor) override;

  /**
   * Calculates the _interactionZones and _interactionSchedule
   * @param identifyZone
   */
  void calculateInteractionSchedule(std::function<std::string(const int[3])> identifyZone);

  /*
   * Class for combining multiple buffers
   */
  class CombinedBuffer {
   public:
    void add_buffer(std::vector<ParticleType> &buffer) { _buffers.push_back(buffer); }

    ParticleType &at(int index) {
      int accSize = _buffers.at(0).get().size();
      size_t bufferIndex = 0;
      while (accSize - 1 < index) {
        ++bufferIndex;
        accSize += _buffers.at(bufferIndex).get().size();
      }
      auto newIndex = index - (accSize - _buffers.at(bufferIndex).get().size());
      if (newIndex >= _buffers.at(bufferIndex).get().size()) {
        throw std::runtime_error("Index out of bounds: " + std::to_string(newIndex) + " from " + std::to_string(index) +
                                 " for buffer of size " + std::to_string(_buffers.at(bufferIndex).get().size()) +
                                 " at buffer index " + std::to_string(bufferIndex) + " of total " +
                                 std::to_string(_buffers.size()) + " with accumulated size " + std::to_string(accSize));
      }
      return _buffers.at(bufferIndex).get().at(newIndex);
    }

    size_t size() {
      return std::accumulate(_buffers.begin(), _buffers.end(), 0, [](auto a, auto b) { return a + b.get().size(); });
    }

   private:
    std::vector<std::reference_wrapper<std::vector<ParticleType>>> _buffers;
  };

  // stores is newton3 is used in the node
  bool _useNewton3;

  // stores if pairwise interaction is used
  bool _pairwiseInteraction;

  // stores cutoff
  double _cutoff;

  // flag to control if brute force schedule is used for triplet interactions
  inline static constexpr bool _bruteForceSchedule3B = false;

  /**
   * Stores the triwise zonal interaction schedule
   */
  std::map<std::string, std::vector<std::set<std::string>>> _interactionScheduleTriwise;

  /**
   * Calculate the interaction schedule for triplet interactions, which considers all
   * zonal interactions which consist of zones that are neighbours.
   * @param identifyZone
   */
  void calculateInteractionScheduleTriwiseBruteForce(std::function<std::string(const int[3])> identifyZone);

  /**
   * Calculate the zonal interaction for triplet interactions, given the bruteforce schedule.
   * @param identifyZone
   */
  void calculateZonalInteractionTriwiseBruteForce(
      std::string zone, std::function<void(ParticleType &, ParticleType &, ParticleType &, bool)> aosFunctor);

  /**
   * Calculate the interaction schedule for triplet interactions, which considers all
   * zonal interactions where at least two zones lie on a distinct "side" or "plane".
   * @param identifyZone
   */
  void calculateInteractionScheduleTriwisePlane(std::function<std::string(const int[3])> identifyZone);

  /**
   * Calculate the zonal interaction for triplet interactions, given the plane schedule.
   * @param identifyZone
   */
  void calculateZonalInteractionTriwisePlane(
      std::string zone, std::function<void(ParticleType &, ParticleType &, ParticleType &, bool)> aosFunctor);

  /**
   * Calculate the relative coordinate of a given zone string
   * @param s
   */
  std::array<int, 3> convZoneStringIntoRelNeighbourIndex(std::string s);
};

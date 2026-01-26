/**
 * @file LiveInfo.h
 * @author humig
 * @date 28.06.2021
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <variant>

#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ParticleBinStructure.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class is able to gather and store important information for a tuning phase from a container and functor.
 * Infos are identified by a string. Their type is InfoType.
 */
class LiveInfo {
  // The class currently needs to be defined in a header only since iteration over a particle container requires to
  // know the Particle type. Actually, no particle specific information can be used here, so only a pointer to
  // ParticleBase would suffice, but iteration doesn't work that way at the moment.

 private:
  /**
   * Returns a particle bin structure that mimics a Linked Cells container with cell size factor 1.
   *
   * @param domainSize Size of (sub)domain.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param boxMin Lower left corner of (sub)domain.
   * @param boxMax Upper right corner of (sub)domain.
   * @param cutoff Cutoff.
   * @return
   */
  static utils::ParticleBinStructure buildCellBinStructure(const std::array<double, 3> &domainSize,
                                                           const double interactionLength,
                                                           const std::array<double, 3> &boxMin,
                                                           const std::array<double, 3> &boxMax, double cutoff) {
    std::array<size_t, 3> cellsPerDim{};
    std::array<double, 3> cellLength{};

    for (int d = 0; d < 3; ++d) {
      // The number of cells is rounded down because the cells will be stretched to fit.
      // std::max to ensure there is at least one cell.
      cellsPerDim[d] = std::max(static_cast<size_t>(std::floor(domainSize[d] / interactionLength)), 1ul);
      cellLength[d] = domainSize[d] / static_cast<double>(cellsPerDim[d]);
    }

    return {cellsPerDim, cellLength, boxMin, boxMax, cutoff};
  }

  /**
   * Returns a bin structure where there are, on average, roughly ten particles per bin, and the bin dimensions are
   * simply a scaling of the domain dimensions. Forces there to be at least one bin, e.g. in the case of no particles.
   *
   * Todo The choice of 10 is arbitrary and probably can be optimized.
   *
   * @param domainSize size of (sub)domain
   * @param numParticles number of particles
   * @param boxMin Lower left corner of (sub)domain.
   * @param boxMax Upper right corner of (sub)domain.
   * @param cutoff Cutoff.
   * @return
   */
  static utils::ParticleBinStructure buildParticleDependentBinStructure(const std::array<double, 3> &domainSize,
                                                                        const size_t numParticles,
                                                                        const std::array<double, 3> &boxMin,
                                                                        const std::array<double, 3> &boxMax,
                                                                        double cutoff) {
    using namespace autopas::utils::ArrayMath::literals;

    const auto targetNumberOfBins = std::max(std::ceil(static_cast<double>(numParticles) / 10.), 1.);
    const auto targetNumberOfBinsPerDim = std::cbrt(targetNumberOfBins);
    // This is probably not an integer, so we floor to get more than 10 particles per bin
    const auto numberOfBinsPerDim = static_cast<size_t>(std::floor(targetNumberOfBinsPerDim));
    const auto binDimensions = domainSize / static_cast<double>(numberOfBinsPerDim);

    return {numberOfBinsPerDim, binDimensions, boxMin, boxMax, cutoff};
  }

  /**
   * Returns a bin structure the domain is divided into 3x3x3 uniform bins.
   *
   * @param domainSize size of (sub)domain
   * @param boxMin Lower left corner of (sub)domain.
   * @param boxMax Upper right corner of (sub)domain.
   * @param cutoff Cutoff.
   * @return
   */
  static utils::ParticleBinStructure buildBlurredBinStructure(const std::array<double, 3> &domainSize,
                                                              const std::array<double, 3> &boxMin,
                                                              const std::array<double, 3> &boxMax, double cutoff) {
    using namespace autopas::utils::ArrayMath::literals;

    const auto binLength = domainSize / 3.;

    return {3, binLength, boxMin, boxMax, cutoff};
  }

 public:
  /**
   * The type of an info.
   */
  using InfoType = std::variant<bool, double, size_t, ContainerOption, TraversalOption, LoadEstimatorOption,
                                DataLayoutOption, Newton3Option>;

  /**
   * Gathers key statistics that define the computational profile of the simulation, in order to provide lower
   * dimensional and human understandable inputs for the tuning strategies.
   *
   * A lot of the information is based on a couple of different spatial bin resolutions:
   * - cells: Bin dimensions are the same as cells in the Linked Cells container with CSF 1.
   * - particleDependentBins: Bin dimensions are designed such that the are approximately 10 particles per bin.
   * - blurredBins: The domain is divided equally into 3x3x3 "blurred" bins
   *
   * @note It is not clear how useful the statistics derived from the particleDependentBins are, as the "resolution"
   * varies depending on the number of particles. E.g. Consider a sparse simulation with a "macroscopic" heterogeneity
   * and a dense simulation with a "microscopic" heterogeneity but at a "macroscopic" level is rather homogeneous.
   * These could have the same particleDependentBin homogeneity (relative std. dev. of bin counts) but would respond
   * to traversals very differently. This could, however, provide a useful metric for homogeneity which is somewhat
   * independent of particle density (e.g. could be useful for determining the best traversal with one statistic
   * independently of what cell-size factor is chosen)
   *
   * Currently, it provides:
   * ---- Bin Independent Statistics ----
   * - numOwnedParticles: The number of particles in the container that are marked as owned.
   * - numHaloParticles: The number of particles in the container that are marked as halo.
   * - cutoff: The configured cutoff radius.
   * - skin: The configured skin radius.
   * - domainSizeX: The size of the domain on the x-axis.
   * - domainSizeY: The size of the domain on the y-axis.
   * - domainSizeZ: The size of the domain on the z-axis.
   * - particleSize: The number of bytes one particle in AoS layout needs.
   * - threadCount: The number of OMP threads that can be used.
   * - rebuildFrequency: The current verlet-rebuild-frequency of the simulation.
   * ---- Cell Statistics ----
   * - numCells: The number of cell-bins in the domain.
   * - numEmptyCells: The number of cell-bins that are empty.
   * - minParticlesPerCell: The minimum number of particles in a cell-bin.
   * - maxParticlesPerCell: The maximum number of particles in a cell-bin.
   * - medianParticlesPerCell: The median number of particles per cell-bin.
   * - lowerQuartileParticlesPerCell: The lower quartile number of particles per cell-bin.
   * - upperQuartileParticlesPerCell: The upper quartile number of particles per cell-bin.
   * - meanParticlesPerCell: The average number of particles per cell-bin.
   * - particlesPerCellStdDev: The standard deviation in the number of particles per cell-bin.
   * - relativeParticlesPerCellStdDev: The, relative to the mean, standard deviation in the number of particles per
   * cell-bin.
   * - estimatedNumNeighborInteractions: Rough estimation of number of neighbor interactions. Assumes that neighboring
   * cell-bins contain roughly the same number of particles. Estimation does not work well if this is not the case.
   * ---- Particle Dependent Bin Statistics ----
   * - particleDependentBinMaxDensity: The maximum density of particles of a particle-dependent bin.
   * - particleDependentBinDensityStdDev: The standard deviation in the density of particles per particle-dependent
   * bin. In earlier versions of AutoPas, this was referred to as "Homogeneity".
   * ---- Blurred Bin Statistics ----
   * - maxParticlesPerBlurredBin: The maximum number of particles in a blurred bin.
   * - minParticlesPerBlurredBin: The minimum number of particles in a blurred bin.
   * - medianParticlesPerBlurredBin: The median number of particles per blurred bin.
   * - lowerQuartileParticlesPerBlurredBin: The lower quartile number of particles per blurred bin.
   * - upperQuartileParticlesPerBlurredBin: The upper quartile number of particles per blurred bin.
   * - meanParticlesPerBlurredBin: The mean number of particles per blurred bin.
   * - particlesPerBlurredBinStdDev: The standard deviation in the number of particles per blurred bin.
   * - relativeParticlesPerBlurredBinStdDev: The, relative to the mean, standard deviation in the number of particles
   * per blurred bin.
   *
   * @tparam Particle_T The type of particle the container stores.
   * @param particleIter Iterator to the first particle in the AutoPas container.
   * @param rebuildFrequency The current verlet rebuild frequency that is used in the simulation.
   * @param numOwnedParticles The number of owned particles in the container.
   * @param boxMin Lower left corner of (sub)domain.
   * @param boxMax Upper right corner of (sub)domain.
   * @param cutoff The cutoff.
   * @param skin The (Verlet) skin.
   */
  template <class Particle_T>
  void gather(ContainerIterator<Particle_T, true, false> particleIter, size_t rebuildFrequency,
              size_t numOwnedParticles, std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff,
              double skin) {
    using namespace utils::ArrayMath::literals;
    using utils::ArrayMath::ceilAndCast;

    // Timer for debugging
    utils::Timer timerGatherLiveInfo;
    timerGatherLiveInfo.start();

    // Aliases and info of particle distribution independent information
    const auto interactionLength = cutoff + skin;

    infos["cutoff"] = cutoff;
    infos["skin"] = skin;
    infos["rebuildFrequency"] = rebuildFrequency;
    infos["particleSize"] = sizeof(Particle_T);
    infos["threadCount"] = static_cast<size_t>(autopas_get_max_threads());

    infos["numOwnedParticles"] = numOwnedParticles;

    const auto domainSize = boxMax - boxMin;
    infos["domainSizeX"] = domainSize[0];
    infos["domainSizeY"] = domainSize[1];
    infos["domainSizeZ"] = domainSize[2];

    // Build bin structures
    auto cellBinStruct = buildCellBinStructure(domainSize, interactionLength, boxMin, boxMax, cutoff);
    auto particleDependentBinStruct =
        buildParticleDependentBinStructure(domainSize, numOwnedParticles, boxMin, boxMax, cutoff);
    auto blurredBinStruct = buildBlurredBinStructure(domainSize, boxMin, boxMax, cutoff);

    infos["numCells"] = cellBinStruct.getNumberOfBins();

    // Count the number of owned particles per bin for each bin structure. Also include total count for halo particles.
    size_t numOwnedParticlesCount = 0;
    size_t numHaloParticlesCount = 0;
    for (; particleIter.isValid(); ++particleIter) {
      if (particleIter->isOwned()) {
        numOwnedParticlesCount++;
        const auto particlePos = particleIter->getR();
        if (utils::inBox(particlePos, boxMin, boxMax)) {
          /*
          TODO: figure out why this fails
          cellBinStruct.countParticle(particlePos);
          particleDependentBinStruct.countParticle(particlePos);
          blurredBinStruct.countParticle(particlePos);
          */
        }
      } else if (particleIter->isHalo()) {
        numHaloParticlesCount++;
      }
    }

    // Sanity Check
    if (numOwnedParticlesCount != numOwnedParticles) {
      //AutoPasLog(ERROR,
      //           "Number of owned particles tracked by AutoPas ({}) does not match number of owned particles "
      //           "counted using the iterator ({}).",
      //           numOwnedParticles, numOwnedParticlesCount);
    }

    infos["numOwnedParticles"] = numOwnedParticlesCount;
    infos["numHaloParticles"] = numHaloParticlesCount;

    // calculate statistics
    cellBinStruct.calculateStatistics();
    particleDependentBinStruct.calculateStatistics();
    blurredBinStruct.calculateStatistics();

    // write cellBinStruct statistics to live info
    infos["numEmptyCells"] = cellBinStruct.getNumEmptyBins();
    infos["maxParticlesPerCell"] = cellBinStruct.getMaxParticlesPerBin();
    infos["minParticlesPerCell"] = cellBinStruct.getMinParticlesPerBin();
    infos["medianParticlesPerCell"] = cellBinStruct.getMedianParticlesPerBin();
    infos["lowerQuartileParticlesPerCell"] = cellBinStruct.getLowerQuartileParticlesPerBin();
    infos["upperQuartileParticlesPerCell"] = cellBinStruct.getUpperQuartileParticlesPerBin();
    infos["meanParticlesPerCell"] = cellBinStruct.getMeanParticlesPerBin();
    infos["relativeParticlesPerCellStdDev"] = cellBinStruct.getRelStdDevParticlesPerBin();
    infos["particlesPerCellStdDev"] = cellBinStruct.getStdDevParticlesPerBin();
    infos["estimatedNumNeighborInteractions"] = cellBinStruct.getEstimatedNumberOfNeighborInteractions();

    // write particle dependent bin statistics to live info
    infos["particleDependentBinMaxDensity"] = particleDependentBinStruct.getMaxDensity();
    infos["particleDependentBinDensityStdDev"] = particleDependentBinStruct.getStdDevDensity();

    // write blurred bin statistics to live info
    infos["maxParticlesPerBlurredBin"] = blurredBinStruct.getMaxParticlesPerBin();
    infos["minParticlesPerBlurredBin"] = blurredBinStruct.getMinParticlesPerBin();
    infos["medianParticlesPerBlurredBin"] = blurredBinStruct.getMedianParticlesPerBin();
    infos["lowerQuartileParticlesPerBlurredBin"] = blurredBinStruct.getLowerQuartileParticlesPerBin();
    infos["upperQuartileParticlesPerBlurredBin"] = blurredBinStruct.getUpperQuartileParticlesPerBin();
    infos["meanParticlesPerBlurredBin"] = blurredBinStruct.getMeanParticlesPerBin();
    infos["relativeParticlesPerBlurredBinStdDev"] = blurredBinStruct.getRelStdDevParticlesPerBin();
    infos["particlesPerBlurredBinStdDev"] = blurredBinStruct.getStdDevParticlesPerBin();

    timerGatherLiveInfo.stop();
    AutoPasLog(DEBUG, "Gathering of LiveInfo took {} ns.", timerGatherLiveInfo.getTotalTime());
  }

  /**
   * Returns a map of all infos.
   * @return A map of all infos.
   */
  [[nodiscard]] const auto &get() const { return infos; }

  /**
   * Gets a single value from the infos that corresponds to the given string.
   * @tparam T return type
   * @param key string key for the infos map
   * @return
   */
  template <typename T>
  T get(const std::string &key) const {
    // Find the key in the map
    const auto it = infos.find(key);

    // If key is not found, log an error
    if (it == infos.end()) {
      AutoPasLog(ERROR, "Key '" + key + "' not found in infos map.");
    }

    // Use std::get<T> to extract the value of the desired type
    try {
      return std::get<T>(it->second);
    } catch (const std::bad_variant_access &e) {
      AutoPasLog(ERROR, "Type mismatch for key '" + key + "'. Requested type does not match the stored type.");
    }
  }

  /**
   * Creates a string containing all live info gathered.
   * @return A string containing all live info gathered.
   */
  [[nodiscard]] std::string toString() const {
    auto typeToString = [](auto type) {
      if constexpr (std::is_same_v<decltype(type), bool> or std::is_same_v<decltype(type), double> or
                    std::is_same_v<decltype(type), size_t>) {
        return std::to_string(type);
      } else if constexpr (std::is_base_of_v<Option<decltype(type)>, decltype(type)>) {
        return type.to_string();
      }
      return std::string{"fail"};
    };
    auto pairToString = [&](const auto &pair) { return pair.first + "=" + std::visit(typeToString, pair.second); };
    std::string res{"Live Info: "};
    if (not infos.empty()) {
      // We initialize the accumulation with the first element. Hence, the accumulation starts at next(begin).
      res += std::accumulate(std::next(infos.begin()), infos.end(), pairToString(*infos.begin()),
                             [&](std::string s, const auto &elem) { return std::move(s) + " " + pairToString(elem); });
    }
    return res;
  }

  /**
   * Stream operator to write the LiveInfo to a stream.
   * @param out
   * @param info
   * @return
   */
  friend std::ostream &operator<<(std::ostream &out, const LiveInfo &info) {
    out << info.infos.size() << ' ';
    for (const auto &[name, val] : info.infos) {
      // val.index here is the index of this value's type in the LiveInfo::InfoType variant type.
      // This is needed for determining the type when reading this file via readIndex.
      out << name << ' ' << val.index() << ' ';
      std::visit([&](const auto &v) { out << v << ' '; }, val);
    }
    return out;
  }

  /**
   * Stream operator to read the LiveInfo in from a stream.
   * @param in
   * @param info
   * @return
   */
  friend std::istream &operator>>(std::istream &in, LiveInfo &info) {
    size_t numElements{0};
    in >> numElements;
    for (size_t i = 0; i < numElements; i++) {
      std::string name;
      in >> name;
      size_t idx{0};
      in >> idx;
      const auto val =
          readIndex<LiveInfo::InfoType>(in, idx, std::make_index_sequence<std::variant_size_v<LiveInfo::InfoType>>());
      info.infos[name] = val;
    }
    return in;
  }

  /**
   * Generate a csv representation containing all values from the toString() method.
   * Since the keys are not necessarily known in advance, this method generates a CSV header and a CSV line.
   * @return A pair of strings in the form of (header, line).
   */
  [[nodiscard]] std::pair<std::string, std::string> getCSVLine() const {
    // match all words (anything that is neither a ' ' or '='), that are followed by a '=',
    // ignoring the 'Live Info: ' prefix
    const auto keyRegex = std::regex("([^= ]+)=[^ ]*");
    // match all words that are preceded by a '='
    const auto valueRegex = std::regex("=([^ ]+)");

    auto searchString = toString();
    // remove leading Live Info:
    searchString = searchString.substr(std::string("Live Info: ").size());

    std::sregex_iterator keyIter(searchString.begin(), searchString.end(), keyRegex);
    std::sregex_iterator valueIter(searchString.begin(), searchString.end(), valueRegex);
    std::sregex_iterator end;

    std::stringstream header;
    std::stringstream line;

    while (keyIter != end) {
      // first submatch is the match of the capture group
      header << keyIter->str(1) << ",";
      ++keyIter;
    }
    while (valueIter != end) {
      // first submatch is the match of the capture group
      line << valueIter->str(1) << ",";
      ++valueIter;
    }

    auto headerStr = header.str();
    auto lineStr = line.str();
    // drop trailing ','
    headerStr.pop_back();
    lineStr.pop_back();

    return std::make_pair(headerStr, lineStr);
  }

 private:
  /**
   * If the template argument Idx is equal to the function argument idx, the method reads in Type from in and stores it
   * in var.
   * @tparam Variant The Variant to return (through var).
   * @tparam Type The type to read in if Idx==idx.
   * @tparam Idx The template argument Idx to compare idx with.
   * @param in The stream to read from.
   * @param idx The index to compare to the template argument Idx.
   * @param var The variant that is filled with the read value if Idx==idx.
   */
  template <class Variant, class Type, size_t Idx>
  static void readIndexHelper(std::istream &in, size_t idx, Variant &var) {
    if (Idx == idx) {
      Type val;
      in >> val;
      var = val;
    }
  }

  /**
   * Reads the value of type std::variant_alternative_t<idx, Variant> from the stream and returns the variant.
   * @tparam Variant The variant that holds the types.
   * @tparam Idx An index sequence for all possible alternatives in the variant.
   * @param in The stream to read from.
   * @param idx The index of the alternative that should be read in.
   * @return A variant that holds the read value of the type indexed by idx in the Variant.
   */
  template <class Variant, size_t... Idx>
  static Variant readIndex(std::istream &in, size_t idx, std::index_sequence<Idx...>) {
    Variant var;
    (readIndexHelper<Variant, std::variant_alternative_t<Idx, Variant>, Idx>(in, idx, var), ...);
    return var;
  }

  /**
   * The map that stores all infos gathered so far.
   */
  std::map<std::string, InfoType> infos;
};

}  // namespace autopas

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

#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/CellBlock3D.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SimilarityFunctions.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
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

 public:
  /**
   * The type of an info.
   */
  using InfoType = std::variant<bool, double, size_t, ContainerOption, TraversalOption, LoadEstimatorOption,
                                DataLayoutOption, Newton3Option>;

  /**
   * Gathers important information from a particle container and functor.
   *
   * The gathered information should allow to estimate the performance of different configurations.
   *
   * Currently, it provides:
   * - numOwnedParticles: The number of particles in the container that are marked as owned.
   * - numHaloParticles: The number of particles in the container that are marked as halo.
   * - numDummyParticles: The number of particles in the container that are marked as dummy.
   * - numTotalParticles: The total number of particles in the container.
   * - cutoff: The configured cutoff radius.
   * - skin: The configured skin radius.
   * - domainSizeX: The size of the domain on the x-axis.
   * - domainSizeY: The size of the domain on the y-axis.
   * - domainSizeZ: The size of the domain on the z-axis.
   * - particleSize: The number of bytes one particle in AoS layout needs.
   * - particleSizeNeededByFunctor: The number of bytes the information needed by the functor from each particle
   * occupies. Important for the SoA data layout, but irrelevant for the AoS data layout.
   * - numCells: The number of cells in the domain if a cell has a side-length equal to the cutoff.
   * - numEmptyCells: The number of empty cells in the domain.
   * - minParticlesPerCell: The minimum number of particles a cell in the domain contains.
   * - maxParticlesPerCell: The maximum number of particles a cell in the domain contains.
   * - avgParticlesPerCell: The average number of particles per cell. (Cells are small so we don't expect outliers
   * that make the average useless).
   * - estimatedNumNeighborInteractions: Rough estimation of number of neighbor interactions. Assumes that neighboring
   * cells contain roughly the same number of particles. Estimation does not work well if this is not the case.
   * - percentParticlesPerCellStdDev: The standard deviation of the number of particles in each cell from the
   * average number of particles per cell, divided by the avgParticlesPerCell.
   * -percentParticlesPerBlurredCellStdDev: The standard deviation of the number of particles in each blurred cell,
   * divided by the average number of particles per blurred cell. A blurred cell is exactly 1/27th of the domain.
   * - threadCount: The number of threads that can be used.
   * - rebuildFrequency: The current verlet-rebuild-frequency of the simulation.
   *
   * @tparam Particle The type of particle the container stores.
   * @tparam PairwiseFunctor The type of functor.
   * @param container The container to gather the infos from.
   * @param functor The functor to gather the infos from.
   * @param rebuildFrequency The current verlet rebuild frequency that is used in the simulation.
   */
  template <class Particle, class PairwiseFunctor>
  void gather(const autopas::ParticleContainerInterface<Particle> &container, const PairwiseFunctor &functor,
              unsigned int rebuildFrequency) {
    using namespace autopas::utils::ArrayMath::literals;
    using autopas::utils::ArrayMath::ceilToInt;
    using autopas::utils::ArrayMath::floorToInt;

    // Some aliases for quicker access
    const auto &boxMin = container.getBoxMin();
    const auto &boxMax = container.getBoxMax();
    const auto cutoff = container.getCutoff();
    const auto skin = container.getVerletSkin();
    const auto interactionLength = cutoff + skin;
    const auto interactionLengthInv = 1. / interactionLength;

    const auto numOwnedParticles = container.getNumberOfParticles(OwnershipState::owned);
    const auto numHaloParticles = container.getNumberOfParticles(OwnershipState::halo);
    const auto numDummyParticles = container.getNumberOfParticles(OwnershipState::dummy);
    const auto numTotalParticles = numOwnedParticles + numHaloParticles + numDummyParticles;

    infos["numOwnedParticles"] = numOwnedParticles;
    infos["numHaloParticles"] = numHaloParticles;
    infos["numDummyParticles"] = numDummyParticles;
    infos["numTotalParticles"] = numTotalParticles;

    infos["cutoff"] = cutoff;
    infos["skin"] = container.getVerletSkin();
    infos["rebuildFrequency"] = static_cast<size_t>(rebuildFrequency);
    const auto domainSize = boxMax - boxMin;

    infos["domainSizeX"] = domainSize[0];
    infos["domainSizeY"] = domainSize[1];
    infos["domainSizeZ"] = domainSize[2];

    infos["particleSize"] = sizeof(Particle);

    // Calculate number of cells for a linked cells container, assuming cell size factor == 1
    const auto [numCells, cellsPerDim, cellLength, cellLengthReciprocal] = [&]() {
      size_t numCellsTmp = 1;
      std::array<size_t, 3> cellsPerDimTmp;
      std::array<double, 3> cellLengthTmp;
      std::array<double, 3> cellLengthReciprocalTmp;
      for (int d = 0; d < 3; ++d) {
        // The number of cells is rounded down because the cells will be stretched to fit.
        // std::max to ensure there is at least one cell.
        cellsPerDimTmp[d] = std::max(static_cast<size_t>(std::floor(domainSize[d] / interactionLength)), 1ul);

//        const auto cellsPerDimWithHalo = cellsPerDimTmp[d] + 2;

        cellLengthTmp[d] = domainSize[d] / static_cast<double>(cellsPerDimTmp[d]);

        cellLengthReciprocalTmp[d] = static_cast<double>(cellsPerDimTmp[d]) / domainSize;
//
//        _haloBoxMin[d] = _boxMin[d] - _cellsPerInteractionLength * _cellLength[d];
//        _haloBoxMax[d] = _boxMax[d] + _cellsPerInteractionLength * _cellLength[d];

        numCellsTmp *= cellsPerDimTmp[d];
      }
      return std::tuple{numCellsTmp, cellsPerDimTmp, cellLengthTmp, cellLengthReciprocalTmp};
    }();

    infos["numCells"] = numCells;

    const auto cellVolume = cellLength[0] * cellLength[1] * cellLength[2];

    // Count how many particles are in each cell via bin counting
    std::vector<size_t> particleCells;
    particleCells.resize(numCells);

    // Blurred bins divide the domain into 3x3x3 equivalent boxes.
    std::vector<size_t> particleBinsBlurred;
    particleBinsBlurred.resize(27);
    const auto blurredBinsDimsReciprocal = std::array<double, 3>{3.0, 3.0, 3.0} / domainSize;

    for (const Particle &particle = container.begin(OwnershipState::owned) ; particle != container.end(); ++particle) {
      if (utils::inBox(particle.getR(), boxMin, boxMax)) {
        // find the "actual" cell
        const auto offsetIntoBox = particle.getR() - boxMin;
        const auto cellIndex3D = floorToInt(offsetIntoBox * cellLengthReciprocal);
        const auto cellIndex1D = utils::ThreeDimensionalMapping::threeToOneD<size_t>(cellIndex3D, cellsPerDim);
        particleCells[cellIndex1D]++;

        // find the blurred bin
        const auto blurredBinIndex3D = floorToInt(offsetIntoBox * blurredBinsDimsReciprocal);
        const auto binIndexBlurred = utils::ThreeDimensionalMapping::threeToOneD<size_t>(blurredBinIndex3D, {3, 3, 3});
        particleBinsBlurred[binIndexBlurred]++;
      }
    }

    // calculate statistics about particle distributions per cell
    const auto avgOwnedParticlesPerCell = numOwnedParticles == 0
                                         ? 0.
                                         : static_cast<double>(numOwnedParticles) / static_cast<double>(numCells);
    const auto avgParticlesPerBlurredCell =
        numOwnedParticles == 0 ? 0.
                          : static_cast<double>(numOwnedParticles) / static_cast<double>(particleBinsBlurred.size());

    // Determine the estimated hit rate if the Linked Cells method with CSF 1 was used
    // This is the ratio of cutoff sphere volume to the volume of 27 cells (number of cells within which one particle
    // could find neighbors given that the cell size factor is 1)
    // This assumes of homogeneous distribution
    const auto volumeOfCutoffSphere = 4. / 3. * M_PI * cutoff * cutoff * cutoff;
    const auto potentialInteractionVolume = cellVolume * 27.;
    const auto estimatedHitRate = volumeOfCutoffSphere / potentialInteractionVolume;

    const auto [estimatedNumNeighborInteractions, maxDiff, sumVariance, numEmptyCells, maxParticlesPerCell,
                minParticlesPerCell] = [&]() {
      double estimatedNumNeighborInteractionsLambda = 0.;
      double maxDiffLambda = 0.;
      double sumVarianceLambda = 0.;
      size_t numEmptyCellsLambda = 0;
      size_t maxParticlesPerCellLambda = std::numeric_limits<size_t>::min();
      size_t minParticlesPerCellLambda = std::numeric_limits<size_t>::max();
      // go over all cells and calculate statistics
      for (size_t i = 0; i < particleCells.size(); i++) {
        const auto particlesInCell = particleCells[i];
        if (particlesInCell == 0) {
          ++numEmptyCellsLambda;
        } else {
          // Assume that the distribution of particles in this cell and the neighboring cells is homogeneous (all cells
          // have the same number of particles as this one and are spread out evenly)
          // At a 3x3x3 cells level, this is maybe okay but the validity of this estimate with this assumption is not tested.

          // Assume no newton3

          // Using the estimateHitRate, the number of particles in this cell, and the above assumption, calculate the
          // number of neighbor interactions.
          // This is
          // - for every particle in this cell: [particlesInCell * ...]
          //   - it interacts with every particle in this cell and neighboring cells [... * (particlesInCell * 27)]
          //   - these interactions have a hit rate of `estimatedHitRate`
          // - remove self interactions [... - particlesInCell]
          // In a very sparse situation, with a large potentialInteractionVolume and small particlesInCell, this could
          // lead to a negative estimate, so just take 0.
          estimatedNumNeighborInteractionsLambda += std::max(
              static_cast<double>(particlesInCell * (particlesInCell * 27)) * estimatedHitRate  - particlesInCell, 0);

        }

        maxParticlesPerCellLambda = std::max(particlesInCell, maxParticlesPerCellLambda);
        minParticlesPerCellLambda = std::min(particlesInCell, minParticlesPerCellLambda);

        const auto diffFromAvg = avgOwnedParticlesPerCell - particlesInCell;
        maxDiffLambda = std::max(diffFromAvg, maxDiffLambda);
        sumVarianceLambda += diffFromAvg * diffFromAvg;
      }

      return std::tuple{estimatedNumNeighborInteractionsLambda,
                        maxDiffLambda,
                        sumVarianceLambda,
                        numEmptyCellsLambda,
                        maxParticlesPerCellLambda,
                        minParticlesPerCellLambda};
    }();

    infos["numEmptyCells"] = numEmptyCells;
    infos["maxParticlesPerCell"] = maxParticlesPerCell;
    infos["minParticlesPerCell"] = minParticlesPerCell;
    infos["particlesPerCellStdDev"] =
        numOwnedParticles == 0 ? 0.
                          : std::sqrt(sumVariance) / static_cast<double>(numCells) / avgOwnedParticlesPerCell;
    infos["avgParticlesPerCell"] = avgOwnedParticlesPerCell;

    const auto [homogeneity, maxDensity] = autopas::utils::calculateHomogeneityAndMaxDensity(container);

    infos["homogeneity"] = homogeneity;
    infos["maxDensity"] = maxDensity;

    infos["estimatedNumNeighborInteractions"] = static_cast<unsigned long>(estimatedNumNeighborInteractions);

    const double sumVarianceBlurred = [&]() {
      double res = 0.;
      for (auto numParticlesInBin : particleBinsBlurred) {
        auto diff = avgParticlesPerBlurredCell - static_cast<int>(numParticlesInBin);
        res += diff * diff;
      }
      return res;
    }();
    infos["particlesPerBlurredCellStdDev"] = numOwnedParticles == 0 ? 0.
                                                               : std::sqrt(sumVarianceBlurred) /
                                                                     static_cast<double>(particleBinsBlurred.size()) /
                                                                     avgParticlesPerBlurredCell;

    infos["threadCount"] = static_cast<size_t>(autopas::autopas_get_max_threads());

    constexpr size_t particleSizeNeededByFunctor = calculateParticleSizeNeededByFunctor<Particle, PairwiseFunctor>(
        std::make_index_sequence<PairwiseFunctor::getNeededAttr().size()>());
    infos["particleSizeNeededByFunctor"] = particleSizeNeededByFunctor;
  }

  /**
   * Returns a map of all infos.
   * @return A map of all infos.
   */
  [[nodiscard]] const auto &get() const { return infos; }

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
   * Private helper to calculate the particle size needed by a functor. This is the sum of the size of the type of all
   * needed attributes.
   * @tparam Particle The type of the Particle
   * @tparam PairwiseFunctor The Functor
   * @tparam Idx An index sequence for all elements of PairwiseFunctor::getNeededAttr() elements.
   * @return The number of bytes needed by the information of a particle the functor needs.
   */
  template <class Particle, class PairwiseFunctor, size_t... Idx>
  constexpr static auto calculateParticleSizeNeededByFunctor(std::index_sequence<Idx...>) {
    return (0 + ... +
            sizeof(typename std::tuple_element<PairwiseFunctor::getNeededAttr()[Idx],
                                               typename Particle::SoAArraysType>::type::value_type));
  }

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

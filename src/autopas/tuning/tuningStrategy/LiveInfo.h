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

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/ArrayMath.h"
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
   * - numParticles: The number of particles in the container.
   * - numHaloParticles: The number of particles in the container that are outside of the domain.
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
    const auto cutoffInv = 1.0 / cutoff;

    infos["numParticles"] = container.getNumberOfParticles();
    infos["cutoff"] = cutoff;
    infos["skin"] = container.getVerletSkin();
    infos["rebuildFrequency"] = static_cast<size_t>(rebuildFrequency);
    const auto domainSize = boxMax - boxMin;

    infos["domainSizeX"] = domainSize[0];
    infos["domainSizeY"] = domainSize[1];
    infos["domainSizeZ"] = domainSize[2];

    infos["particleSize"] = sizeof(Particle);

    // Calculate number of cells for a linked cells container, assuming cell size factor == 1
    const auto cellsPerDim = ceilToInt(domainSize * cutoffInv);
    const auto numCells = static_cast<size_t>(cellsPerDim[0] * cellsPerDim[1] * cellsPerDim[2]);
    infos["numCells"] = numCells;

    // Count how many particles are in each cell via bin counting
    std::vector<size_t> particleBins;
    // +1 because we count all halo particles in the last bin
    particleBins.resize(numCells + 1);
    // Blurred cells divide the domain into 3x3x3 equivalent boxes.
    std::vector<size_t> particleBinsBlurred;
    particleBinsBlurred.resize(27);
    const auto blurredCellDimsReciproc = std::array<double, 3>{3.0, 3.0, 3.0} / domainSize;
    for (const Particle &particle : container) {
      if (utils::inBox(particle.getR(), boxMin, boxMax)) {
        // find the actual cell
        const auto offsetIntoBox = particle.getR() - boxMin;
        const auto cell = floorToInt(offsetIntoBox * cutoffInv);
        const auto binIndex = (cell[2] * cellsPerDim[1] + cell[1]) * cellsPerDim[0] + cell[0];
        particleBins[binIndex]++;

        // find the blurred cell
        const auto cellBlurred = floorToInt(offsetIntoBox * blurredCellDimsReciproc);
        const auto binIndexBlurred = (cellBlurred[2] * 3 + cellBlurred[1]) * 3 + cellBlurred[0];
        particleBinsBlurred[binIndexBlurred]++;
      } else {
        // found a halo particle
        particleBins.back()++;
      }
    }

    infos["numHaloParticles"] = particleBins.back();

    // calculate statistics about particle distributions per cell
    const auto avgParticlesPerCell =
        static_cast<double>(container.getNumberOfParticles() - particleBins.back()) / static_cast<double>(numCells);
    const auto avgParticlesPerBlurredCell =
        static_cast<double>(container.getNumberOfParticles() - particleBins.back()) /
        static_cast<double>(particleBinsBlurred.size());

    const auto [estimatedNumNeighborInteractions, maxDiff, sumStddev, numEmptyCells, maxParticlesPerCell,
                minParticlesPerCell] = [&]() {
      double estimatedNumNeighborInteractionsLambda = 0.;
      double maxDiffLambda = 0.;
      double sumStddevLambda = 0.;
      size_t numEmptyCellsLambda = 0;
      size_t maxParticlesPerCellLambda = std::numeric_limits<size_t>::min();
      size_t minParticlesPerCellLambda = std::numeric_limits<size_t>::max();
      // go over all bins and calculate statistics
      for (size_t i = 0; i < particleBins.size() - 1; i++) {
        const auto particlesInBin = particleBins[i];
        if (particlesInBin == 0) {
          ++numEmptyCellsLambda;
        } else {
          // FIXME: Tobias explain this calculation (:
          estimatedNumNeighborInteractionsLambda +=
              static_cast<double>(particlesInBin * (particlesInBin * 27 - 1)) * 0.155;
        }
        maxParticlesPerCellLambda = std::max(particlesInBin, maxParticlesPerCellLambda);
        minParticlesPerCellLambda = std::min(particlesInBin, minParticlesPerCellLambda);
        const auto diffFromAvg = avgParticlesPerCell - static_cast<int>(particlesInBin);
        maxDiffLambda = std::max(diffFromAvg, maxDiffLambda);
        sumStddevLambda += diffFromAvg * diffFromAvg;
      }
      estimatedNumNeighborInteractionsLambda /= 2;

      return std::tuple{estimatedNumNeighborInteractionsLambda,
                        maxDiffLambda,
                        sumStddevLambda,
                        numEmptyCellsLambda,
                        maxParticlesPerCellLambda,
                        minParticlesPerCellLambda};
    }();

    infos["numEmptyCells"] = numEmptyCells;
    infos["maxParticlesPerCell"] = maxParticlesPerCell;
    infos["minParticlesPerCell"] = minParticlesPerCell;
    infos["particlesPerCellStdDev"] =
        std::sqrt(sumStddev) / static_cast<double>(particleBins.size() - 1) / avgParticlesPerCell;
    infos["avgParticlesPerCell"] = avgParticlesPerCell;
    infos["estimatedNumNeighborInteractions"] = static_cast<unsigned long>(estimatedNumNeighborInteractions);

    const double sumStddevBlurred = [&]() {
      double res = 0.;
      for (auto numParticlesInBin : particleBinsBlurred) {
        auto diff = avgParticlesPerBlurredCell - static_cast<int>(numParticlesInBin);
        res += diff * diff;
      }
      return res;
    }();
    infos["particlesPerBlurredCellStdDev"] =
        std::sqrt(sumStddevBlurred) / static_cast<double>(particleBinsBlurred.size()) / avgParticlesPerBlurredCell;

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
/**
 * @file LiveInfo.h
 * @author humig
 * @date 28.06.2021
 */

#pragma once

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
   * - percentParticlesPerCellStdDev: The standard deviation of the number of particles in each cell from the
   * average number of particles per cell, divided by the avgParticlesPerCell.
   * -percentParticlesPerBlurredCellStdDev: The standard deviation of the number of particles in each blurred cell,
   * divided by the average number of particles per blurred cell. A blurred cell is exactly 1/27th of the domain.
   * - threadCount: The number of threads that can be used.
   *
   * @tparam Particle The type of particle the container stores.
   * @tparam PairwiseFunctor The type of functor.
   * @param container The container to gather the infos from.
   * @param functor The functor to gather the infos from.
   */
  template <class Particle, class PairwiseFunctor>
  void gather(const autopas::ParticleContainerInterface<Particle> &container, const PairwiseFunctor &functor) {
    infos["numParticles"] = container.getNumParticles();
    infos["cutoff"] = container.getCutoff();
    infos["skin"] = container.getSkin();
    auto domainSize = utils::ArrayMath::sub(container.getBoxMax(), container.getBoxMin());

    infos["domainSizeX"] = domainSize[0];
    infos["domainSizeY"] = domainSize[1];
    infos["domainSizeZ"] = domainSize[2];

    infos["particleSize"] = sizeof(Particle);

    using namespace utils::ArrayMath;
    auto cellsPerDim = ceilToInt(mulScalar(domainSize, 1.0 / container.getCutoff()));
    auto numCells = static_cast<size_t>(cellsPerDim[0] * cellsPerDim[1] * cellsPerDim[2]);
    infos["numCells"] = numCells;

    std::vector<size_t> particleBins;
    particleBins.resize(numCells + 1);
    // Each blurred bin is exactly 1/27th of the domain
    std::vector<size_t> particleBinsBlurred;
    particleBinsBlurred.resize(27);
    auto blurredCellDimsReciproc = div({3.0, 3.0, 3.0}, sub(container.getBoxMax(), container.getBoxMin()));
    for (const Particle &particle : container) {
      if (utils::inBox(particle.getR(), container.getBoxMin(), container.getBoxMax())) {
        auto offset = sub(particle.getR(), container.getBoxMin());
        auto cell = floorToInt(mulScalar(offset, 1.0 / container.getCutoff()));
        auto idx = (cell[2] * cellsPerDim[1] + cell[1]) * cellsPerDim[0] + cell[0];
        particleBins.at(idx)++;

        auto cellBlurred = floorToInt(mul(offset, blurredCellDimsReciproc));
        auto idxBlurred = (cellBlurred[2] * 3 + cellBlurred[1]) * 3 + cellBlurred[0];
        particleBinsBlurred.at(idxBlurred)++;
      } else {
        particleBins.back()++;
      }
    }

    infos["numHaloParticles"] = particleBins.back();

    auto avg = static_cast<double>(container.getNumParticles() - particleBins.back()) / static_cast<double>(numCells);
    auto avgBlurred = static_cast<double>(container.getNumParticles() - particleBins.back()) / 27;
    double maxDiff = 0;
    double sumStddev = 0;
    size_t numEmptyCells = 0;
    size_t maxParticlesPerCell = 0;
    size_t minParticlesPerCell = std::numeric_limits<size_t>::max();
    for (size_t i = 0; i < particleBins.size() - 1; i++) {
      size_t particlesInBin = particleBins[i];
      if (particlesInBin == 0) {
        numEmptyCells++;
      }
      if (particlesInBin > maxParticlesPerCell) {
        maxParticlesPerCell = particlesInBin;
      }
      if (particlesInBin < minParticlesPerCell) {
        minParticlesPerCell = particlesInBin;
      }
      auto diff = avg - static_cast<int>(particlesInBin);
      if (diff > maxDiff) {
        maxDiff = diff;
      }
      sumStddev += diff * diff;
    }
    infos["numEmptyCells"] = numEmptyCells;
    infos["maxParticlesPerCell"] = maxParticlesPerCell;
    infos["minParticlesPerCell"] = minParticlesPerCell;
    infos["particlesPerCellStdDev"] = std::sqrt(sumStddev) / static_cast<double>(particleBins.size() - 1) / avg;
    infos["avgParticlesPerCell"] = avg;

    double sumStddevBlurred = 0;
    for (auto numParticlesInBin : particleBinsBlurred) {
      auto diff = avgBlurred - static_cast<int>(numParticlesInBin);
      sumStddevBlurred += diff * diff;
    }

    infos["particlesPerBlurredCellStdDev"] =
        std::sqrt(sumStddevBlurred) / static_cast<double>(particleBinsBlurred.size()) / avgBlurred;

    infos["threadCount"] = static_cast<size_t>(autopas::autopas_get_max_threads());

    constexpr size_t particleSizeNeededByFunctor = calculateParticleSizeNeededByFunctor<Particle, PairwiseFunctor>(
        std::make_index_sequence<PairwiseFunctor::getNeededAttr().size()>());
    infos["particleSizeNeededByFunctor"] = particleSizeNeededByFunctor;
  }

  /**
   * Returns a map of all infos.
   * @return A map of all infos.
   */
  [[nodiscard]] const auto &get() { return infos; }

  /**
   * Creates a string containing all live info gathered.
   * @return A string containing all live info gathered.
   */
  [[nodiscard]] std::string toString() const {
    std::string res{"Live Info: "};
    auto typeToString = [](auto type) {
      if constexpr (std::is_same_v<decltype(type), bool> or std::is_same_v<decltype(type), double> or
                    std::is_same_v<decltype(type), size_t>) {
        return std::to_string(type);
      } else if constexpr (std::is_base_of_v<Option<decltype(type)>, decltype(type)>) {
        return type.to_string();
      }
      return std::string{"fail"};
    };
    auto toString = [&](const auto &pair) { return pair.first + "=" + std::visit(typeToString, pair.second); };
    if (not infos.empty()) {
      res += std::accumulate(std::next(infos.begin()), infos.end(), toString(*infos.begin()),
                             [&](std::string s, const auto &elem) { return std::move(s) + " " + toString(elem); });
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
    size_t numElements;
    in >> numElements;
    for (size_t i = 0; i < numElements; i++) {
      std::string name;
      in >> name;
      size_t idx;
      in >> idx;
      auto val =
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

 private:
  /**
   * The map that stores all infos gathered so far.
   */
  std::map<std::string, InfoType> infos;
};

}  // namespace autopas
/**
* @file LiveInfo.h
 * @author humig
 * @date 28.06.2021
*/

#pragma once

#include <variant>

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

class LiveInfo {
 public:
  LiveInfo() = default;

  using InfoType = std::variant<bool, double, size_t, ContainerOption, TraversalOption, LoadEstimatorOption,
  DataLayoutOption, Newton3Option>;

  template<class Particle, class PairwiseFunctor>
  void gather(const autopas::ParticleContainerInterface<Particle>& container, const PairwiseFunctor& functor) {
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
    for(const Particle& particle : container) {
      if(utils::inBox(particle.getR(), container.getBoxMin(), container.getBoxMax())) {
        auto offset = sub(particle.getR(), container.getBoxMin());
        auto cell = floorToInt(mulScalar(offset, 1.0 / container.getCutoff()));
        auto idx = (cell[2] * cellsPerDim[1] + cell[1]) * cellsPerDim[0] + cell[0];
        particleBins.at(idx)++;
      } else {
        particleBins.back()++;
        std::cout << utils::ArrayUtils::to_string(particle.getR()) << std::endl;
      }
    }

    infos["numHaloParticles"] = particleBins.back();

    auto avg = static_cast<double>(std::accumulate(particleBins.begin(), particleBins.end()-1, 0ul)) /
               static_cast<double>(particleBins.size() - 1);
    double maxDiff = 0;
    double sum = 0;
    size_t numEmptyCells = 0;
    size_t maxParticlesPerCell = 0;
    size_t minParticlesPerCell = std::numeric_limits<size_t>::max();
    for(size_t i = 0; i < particleBins.size() - 1; i++) {
      size_t particlesInBin = particleBins[i];
      if(particlesInBin == 0) {
        numEmptyCells++;
      }
      if(particlesInBin > maxParticlesPerCell) {
        maxParticlesPerCell = particlesInBin;
      }
      if(particlesInBin < minParticlesPerCell) {
        minParticlesPerCell = particlesInBin;
      }
      auto diff = avg - static_cast<int>(particlesInBin);
      if(diff > maxDiff) {
        maxDiff = diff;
      }
      sum += std::sqrt(diff*diff);
    }
    infos["numEmptyCells"] = numEmptyCells;
    infos["maxParticlesPerCell"] = maxParticlesPerCell;
    infos["minParticlesPerCell"] = minParticlesPerCell;
    infos["meanSquaredDiffParticlesPerCell"] = sum / static_cast<double>(particleBins.size() - 1);
    infos["avgParticlesPerCell"] = static_cast<double>(container.getNumParticles() - particleBins.back()) / numCells;

    infos["threadCount"] = static_cast<size_t>(autopas::autopas_get_max_threads());

    constexpr size_t particleSizeNeededByFunctor = calculateParticleSizeNeededByFunctor<Particle, PairwiseFunctor>(
        std::make_index_sequence<PairwiseFunctor::getNeededAttr().size()>()
        );
    infos["particleSizeNeededByFunctor"] = particleSizeNeededByFunctor;
  }

  [[nodiscard]] const auto& get() {
    return infos;
  }

  [[nodiscard]] std::string toString() {
    std::string res{"Live Info: "};
    auto typeToString = [](auto type) {
        if constexpr (std::is_same_v<decltype(type), bool> or
        std::is_same_v<decltype(type), double> or std::is_same_v<decltype(type), size_t>) {
          return std::to_string(type);
        } else if constexpr( std::is_base_of_v<Option<decltype(type)>, decltype(type)>) {
          return type.to_string();
        }
        return std::string{"fail"};
    };
    auto toString = [&](const auto& pair) {
      return pair.first + "=" + std::visit(typeToString, pair.second);};
    if(not infos.empty()) {
      res += std::accumulate(std::next(infos.begin()), infos.end(), toString(*infos.begin()),
                             [&](std::string s, const auto& elem) {
          return std::move(s) + " " + toString(elem);
      });
    }
    return res;
  }

 private:
  std::map<std::string, InfoType> infos;
};

} // namespace autopas
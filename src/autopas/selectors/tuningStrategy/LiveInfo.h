/**
* @file LiveInfo.h
 * @author humig
 * @date 28.06.2021
*/

#pragma once

#include <variant>

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

class LiveInfo {
 public:
  LiveInfo() = default;

  using InfoType = std::variant<bool, double, size_t, ContainerOption, TraversalOption, LoadEstimatorOption,
  DataLayoutOption, Newton3Option>;

  template<class Particle>
  void gather(const autopas::ParticleContainerInterface<Particle>& container) {
    infos["numParticles"] = container.getNumParticles();
    infos["cutoff"] = container.getCutoff();
    auto domainSize = utils::ArrayMath::sub(container.getBoxMax(), container.getBoxMin());

    auto cellsPerDim = utils::ArrayMath::mulScalar(domainSize, 1.0 / container.getCutoff());
    auto numCells = static_cast<unsigned long>(std::ceil(cellsPerDim[0]) * std::ceil(cellsPerDim[1]) *
        std::ceil(cellsPerDim[2]));
    infos["numCells"] = numCells;

    infos["avgParticlesPerCell"] = static_cast<double>(container.getNumParticles() / numCells);

    std::vector<size_t> particleBins;
    particleBins.resize(numCells + 1);
    for(const Particle& particle : container) {
      auto cell = utils::ArrayMath::mulScalar(utils::ArrayMath::add(particle.getR(), container.getBoxMin()),
                                              1.0 / container.getCutoff());
      auto idx = cell[0] * cellsPerDim[1] * cellsPerDim[2] + cell[1] * cellsPerDim[2] + cell[2];
      if(idx > 0 and idx < particleBins.size()) {
        particleBins[idx]++;
      } else {
        particleBins.back()++;
      }
    }
    auto avg = static_cast<double>(std::accumulate(particleBins.begin(), particleBins.end()-1, 0ul)) /
               static_cast<double>(particleBins.size() - 1);
    double maxDiff = 0;
    double sum = 0;
    for(size_t particlesInBin : particleBins) {
      auto diff = avg - static_cast<int>(particlesInBin);
      if(diff > maxDiff) {
        maxDiff = diff;
      }
      sum += std::sqrt(diff*diff);
    }
    infos["meanSquaredDiffParticlesPerCell"] = sum / static_cast<double>(particleBins.size() - 1);
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
          return std::move(s) + toString(elem);
      });
    }
    return res;
  }

 private:
  std::map<std::string, InfoType> infos;
};

} // namespace autopas
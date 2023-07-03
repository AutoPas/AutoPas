/**
 * @file SearchSet.cpp
 * @author F. Gratl
 * @date 26.06.23
 */

#include "SearchSet.h"

#include <algorithm>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/SearchSetIterator.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberInterval.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/NumberSetFinite.h"

namespace autopas {
std::string SearchSet::toString() const {
  using autopas::utils::ArrayUtils::operator<<;
  std::stringstream ss;
  // clang-format off
  ss << "SearchSet:"
     << "\nCellSizeFactors: [" << cellSizeFactors->to_string() << "]"
     << "\nConfigurations: " << configurations;
  // clang-format on
  return ss.str();
}

SearchSetIterator SearchSet::begin(double csfStepSize) { return SearchSetIterator(*this, csfStepSize); }

SearchSetIterator SearchSet::end() { return {configurations.end(), cellSizeFactors->getMax(), *this}; }

bool SearchSet::empty() const { return configurations.empty() or cellSizeFactors->isEmpty(); }

std::vector<SearchSet> SearchSet::deleteConfig(const Configuration &configuration) const {
  // find the config is in this set
  std::vector<SearchSet> returnSets{*this};
  const auto confIter = std::find_if(returnSets.back().configurations.begin(), returnSets.back().configurations.end(),
                                     [&](auto conf) { return configuration.equalsDiscreteOptions(conf); });

  // return if the config is not in this set.
  if (confIter == returnSets.back().configurations.end()) {
    return std::vector<SearchSet>{*this};
  }

  // Create the set that contains the full number set but not the configuration.
  if (configurations.size() != 1) {
    returnSets.back().configurations.erase(confIter);
  } else {
    // Don't keep sets with an empty list of configurations
    returnSets.pop_back();
  }

  // Create the set(s) that contains only the configuration and not the csf.
  if (cellSizeFactors->isInterval()) {
    // special case to avoid: the csf to remove is the only number in the interval. In that case don't add anything.
    // special case: the csf to remove is at the edge of the interval then only add one of the sets.
    if (configuration.cellSizeFactor != cellSizeFactors->getMin()) {
      // set with [n0-csf[
      returnSets.emplace_back<SearchSet>(
          {{configuration},
           std::make_unique<NumberInterval<double>>(
               cellSizeFactors->getMin(), std::nextafter(configuration.cellSizeFactor, cellSizeFactors->getMin()))});
    }
    if (configuration.cellSizeFactor != cellSizeFactors->getMax()) {
      // set with ]csf-nX]
      returnSets.emplace_back<SearchSet>(
          {{configuration},
           std::make_unique<NumberInterval<double>>(
               std::nextafter(configuration.cellSizeFactor, cellSizeFactors->getMax()), cellSizeFactors->getMax())});
    }
  } else {
    auto numbersWithoutCsf = cellSizeFactors->getAll();
    numbersWithoutCsf.erase(configuration.cellSizeFactor);
    returnSets.emplace_back<SearchSet>({{configuration}, std::make_unique<NumberSetFinite<double>>(numbersWithoutCsf)});
  }

  return returnSets;
}  // namespace autopas

void SearchSet::swap(SearchSet &other) {
  std::swap(this->configurations, other.configurations);
  std::swap(this->cellSizeFactors, other.cellSizeFactors);
}

SearchSet &SearchSet::operator=(const SearchSet &rhs) {
  *this = SearchSet(rhs);
  return *this;
}

SearchSet &SearchSet::operator=(SearchSet &&rhs) noexcept {
  this->swap(rhs);
  return *this;
}
}  // namespace autopas

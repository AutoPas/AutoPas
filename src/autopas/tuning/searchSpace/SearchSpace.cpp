/**
 * @file SearchSpace.cpp
 * @author F. Gratl
 * @date 27.06.23
 */

#include "SearchSpace.h"

namespace autopas {

SearchSpace::SearchSpace(const std::vector<SearchSet> &searchSets)
    : _searchSets(searchSets),
      // assume that all sets have either continuous or finite cell size factors
      _csfsAreContinuous(not searchSets.empty() and searchSets.front().cellSizeFactors->isInterval()) {}

void SearchSpace::simplifyInternalSets() {
  // repeat folding of sets until nothing changes anymore
  bool setsChanged = true;
  while (setsChanged) {
    setsChanged = false;
    for (int i = 0; i < _searchSets.size(); ++i) {
      for (int j = 0; j < _searchSets.size(); ++j) {
        // helper function to avoid code duplication
        auto commonActionsAfterMerge = [&]() {
          // remove set j
          _searchSets.erase(_searchSets.begin() + j);
          // decrease i and j, so they stay in place for the next iteration
          --i;
          --j;
          setsChanged = true;
        };

        // Case: Same configs -> merge CSFs if possible
        if (_searchSets[i].configurations == _searchSets[j].configurations) {
          if (_csfsAreContinuous) {
            const auto *csfsI = dynamic_cast<NumberInterval<double> *>(_searchSets[i].cellSizeFactors.get());
            const auto *csfsJ = dynamic_cast<NumberInterval<double> *>(_searchSets[j].cellSizeFactors.get());
            if (csfsI->overlaps(*csfsJ)) {
              // replace set i with the merged set
              _searchSets[i] =
                  SearchSet(_searchSets[i].configurations,
                            std::make_unique<NumberInterval<double>>(std::min(csfsI->getMin(), csfsJ->getMin()),
                                                                     std::max(csfsI->getMax(), csfsJ->getMax())));
            }
          } else {
            const auto csfsI = _searchSets[i].cellSizeFactors->getAll();
            const auto csfsJ = _searchSets[j].cellSizeFactors->getAll();
            std::set<double> combinedCsfs{csfsI};
            combinedCsfs.insert(csfsJ.cbegin(), csfsJ.cend());
            // replace set i with the merged set
            _searchSets[i] =
                SearchSet(_searchSets[i].configurations, std::make_unique<NumberSetFinite<double>>(combinedCsfs));
          }
          commonActionsAfterMerge();
        } else if (_searchSets[i].cellSizeFactors == _searchSets[j].cellSizeFactors) {
          // Case: Same CSFs -> merge Configs
          auto combinedConfigs = _searchSets[i].configurations;
          combinedConfigs.insert(_searchSets[j].configurations.begin(), _searchSets[j].configurations.end());
          // replace set i with the merged set
          _searchSets[i] = SearchSet(combinedConfigs, _searchSets[i].cellSizeFactors->clone());
          commonActionsAfterMerge();
        }
      }
    }
  }
}

}  // namespace autopas
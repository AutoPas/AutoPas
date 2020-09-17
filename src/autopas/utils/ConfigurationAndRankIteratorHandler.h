/**
 * @file ConfigurationAndRankIteratorHandler.h
 * @author W. Thieme
 * @date 21.06.2020
 */

#pragma once

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"

namespace autopas::utils {
/**
 * Functionality needed to iterate through ranks and configurations simultaneously in an evenly distributed manner.
 * If more configurations than ranks exist, each configuration is paired with a single rank.
 * If more ranks exist, it's the other way around.
 * Since it cannot traverse infinite sets of CellSizeFactors, it provides information for splitting those among
 * configurations if more ranks than configurations exist.
 */
class ConfigurationAndRankIteratorHandler {
 public:
  /**
   * Constructor for ConfigurationAndRankIteratorHandler.
   * @param containerOptions
   * @param cellSizeFactors
   * @param traversalOptions
   * @param loadEstimatorOptions
   * @param dataLayoutOptions
   * @param newton3Options
   * @param numConfigs
   * @param commSize
   */
  ConfigurationAndRankIteratorHandler(const std::set<ContainerOption> &containerOptions,
                                      const std::set<double> &cellSizeFactors,
                                      const std::set<TraversalOption> &traversalOptions,
                                      const std::set<LoadEstimatorOption> &loadEstimatorOptions,
                                      const std::set<DataLayoutOption> &dataLayoutOptions,
                                      const std::set<Newton3Option> &newton3Options, const int numConfigs,
                                      const int commSize)
      : _containers(containerOptions),
        _cellSizeFactors(cellSizeFactors),
        _allowedTraversalOptions(traversalOptions),
        _allowedLoadEstimatorOptions(loadEstimatorOptions),
        _dataLayoutOptions(dataLayoutOptions),
        _newton3Options(newton3Options) {
    reset(numConfigs, commSize);
  }

  /**
   * Advances the rankIterator (getRankIterator()) and/or the Option iterators for a single step such that repeated.
   * execution of this function ends up in both reaching their respective ends simultaneously.
   * @param numConfigs
   * @param commSize
   */
  void advanceIterators(int numConfigs, int commSize);

  /**
   * Resets all iterators to the beginning, while keeping the containers.
   * @param numConfigs
   * @param commSize
   */
  void reset(const int numConfigs, const int commSize);

  /**
   * Alternative getter for all Configuration iterators.
   * @param containerIt out
   * @param cellSizeFactorIt out
   * @param loadEstimatorIt out
   * @param traversalIt out
   * @param dataLayoutIt out
   * @param newton3It out
   */
  inline void getConfigIterators(std::set<ContainerOption>::iterator &containerIt,
                                 std::set<double>::iterator &cellSizeFactorIt,
                                 std::set<TraversalOption>::iterator &traversalIt,
                                 std::set<LoadEstimatorOption>::iterator &loadEstimatorIt,
                                 std::set<DataLayoutOption>::iterator &dataLayoutIt,
                                 std::set<Newton3Option>::iterator &newton3It) {
    containerIt = _containerIt;
    cellSizeFactorIt = _cellSizeFactorIt;
    traversalIt = _traversalIt;
    loadEstimatorIt = _loadEstimatorIt;
    dataLayoutIt = _dataLayoutIt;
    newton3It = _newton3It;
  }

  /**
   * Getter for the rankIterator.
   * The value will correspond to the rank that holds the Options that the other iterators
   * point to.
   * @return
   */
  [[nodiscard]] inline int getRankIterator() const { return _rankIterator; }

  /**
   * Getter for the number of ranks smaller than getRankIterator that have the exact same configs assigned to them.
   * Used to reduce infinite CSFs on a rank if more ranks than configurations exist.
   * @return
   */
  [[nodiscard]] inline int getInfiniteCellSizeFactorsOffset() const { return _infiniteCellSizeFactorsOffset; }

  /**
   * Getter for the number of ranks in total that have the exact same configs assigned to them.
   * Used to reduce infinite CSFs on a rank if more ranks than configurations exist.
   * @return
   */
  [[nodiscard]] inline int getInfiniteCellSizeFactorsBlockSize() const { return _infiniteCellSizeFactorsBlockSize; }

  /**
   * Getter for containerIterator.
   * @return
   */
  [[nodiscard]] inline std::set<ContainerOption>::iterator getContainerIterator() const { return _containerIt; }

  /**
   * Getter for the CellSizeFactorIterator.
   * @return
   */
  [[nodiscard]] inline std::set<double>::iterator getCellSizeFactorIterator() const { return _cellSizeFactorIt; }

  /**
   * Getter for the TraversalIterator.
   * @return
   */
  [[nodiscard]] inline std::set<TraversalOption>::iterator getTraversalIterator() const { return _traversalIt; }

  /**
   * Getter for the LoadEstimatorIterator.
   * @return
   */
  [[nodiscard]] inline std::set<LoadEstimatorOption>::iterator getLoadEstimatorIterator() const {
    return _loadEstimatorIt;
  }

  /**
   * Getter for the DataLayoutIterator.
   * @return
   */
  [[nodiscard]] inline std::set<DataLayoutOption>::iterator getDataLayoutIterator() const { return _dataLayoutIt; }

  /**
   * Getter for the Newton3Iterator.
   * @return
   */
  [[nodiscard]] inline std::set<Newton3Option>::iterator getNewton3Iterator() const { return _newton3It; }

 private:
  /**
   * Resets _allowedAndApplicableTraversalOptions to fit with whatever _containerOptions is pointing to.
   */
  void selectTraversalsForCurrentContainer();

  /**
   * Resets _allowedAndApplicableLoadEstimatorOptions to fit with the current container and traversal.
   */
  void selectLoadEstimatorsForCurrentContainerAndTraversal();

  /**
   * Selects the next config from containers X cellSizeFactors X traversals X dataLayouts X newton3Options.
   * Later named options are changed first, because they are easier to switch between simulation steps.
   */
  inline void advanceConfigIterators();

  const std::set<ContainerOption> &_containers;
  const std::set<double> &_cellSizeFactors;
  const std::set<TraversalOption> &_allowedTraversalOptions;
  const std::set<LoadEstimatorOption> &_allowedLoadEstimatorOptions;
  const std::set<DataLayoutOption> &_dataLayoutOptions;
  const std::set<Newton3Option> &_newton3Options;
  std::set<TraversalOption> _allowedAndApplicableTraversalOptions;
  std::set<LoadEstimatorOption> _allowedAndApplicableLoadEstimatorOptions;
  std::set<ContainerOption>::iterator _containerIt;
  std::set<double>::iterator _cellSizeFactorIt;
  std::set<TraversalOption>::iterator _traversalIt;
  std::set<LoadEstimatorOption>::iterator _loadEstimatorIt;
  std::set<DataLayoutOption>::iterator _dataLayoutIt;
  std::set<Newton3Option>::iterator _newton3It;
  int _rankIterator;
  int _remainingBlockSize;
  int _remainder;
  int _infiniteCellSizeFactorsOffset;
  int _infiniteCellSizeFactorsBlockSize;
};

}  // namespace autopas::utils
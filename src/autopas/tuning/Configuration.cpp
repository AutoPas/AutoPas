/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 17. Nov 2022
 */

#include "Configuration.h"

#include <optional>

#include "autopas/utils/StaticCellSelector.h"
#include "autopas/utils/StringUtils.h"
#include "tuning/selectors/TraversalSelector.h"

std::string autopas::Configuration::toString() const {
  return "{Interaction Type: " + interactionType.to_string() + " , Container: " + container.to_string() +
         " , CellSizeFactor: " + std::to_string(cellSizeFactor) + " , Traversal: " + traversal.to_string() +
         " , Load Estimator: " + loadEstimator.to_string() + " , Data Layout: " + dataLayout.to_string() +
         " , Newton 3: " + newton3.to_string() + "}";
}

std::string autopas::Configuration::getCSVHeader() const { return getCSVRepresentation(true); }

std::string autopas::Configuration::getCSVLine() const { return getCSVRepresentation(false); }

bool autopas::Configuration::hasValidValues() const {
  return container != ContainerOption() and cellSizeFactor != -1 and traversal != TraversalOption() and
         loadEstimator != LoadEstimatorOption() and dataLayout != DataLayoutOption() and newton3 != Newton3Option() and
         interactionType != InteractionTypeOption();
}

std::string autopas::Configuration::getCSVRepresentation(bool returnHeaderOnly) const {
  // Escape for '{' and '}' required when using Apple Clang 15.0.0
  auto rgx = returnHeaderOnly ?
                              // match any sequence before a colon and drop any spaces, comma or brackets before it
                 std::regex("[\\{, ]+([^:]+):[^,]*")
                              :
                              // match any sequence after a colon and drop any spaces, comma or brackets around it
                 std::regex(": ([^,]+)(?: ,|\\})");
  auto searchString = toString();
  std::sregex_iterator matchIter(searchString.begin(), searchString.end(), rgx);
  std::sregex_iterator end;
  std::stringstream retStream;

  while (matchIter != end) {
    // first submatch is the match of the capture group
    retStream << matchIter->str(1) << ",";
    ++matchIter;
  }
  auto retString = retStream.str();
  // drop trailing ','
  retString.pop_back();
  return retString;
}

bool autopas::Configuration::hasCompatibleValues(bool silent) const {
  // Check if container and traversal fit together
  const auto &allContainerTraversals = compatibleTraversals::allCompatibleTraversals(container, interactionType);
  if (allContainerTraversals.find(traversal) == allContainerTraversals.end()) {
    if (not silent) {
      AutoPasLog(WARN, "Configuration is incompatible: Container {} does not support the traversal {} with "
                 "interaction type {}.", container, traversal, interactionType);
    }
    return false;
  }

  // Check if the selected load estimator option is applicable.
  const std::set<LoadEstimatorOption> applicableLoadEstimators =
      loadEstimators::getApplicableLoadEstimators(container, traversal, LoadEstimatorOption::getAllOptions());
  if (applicableLoadEstimators.find(loadEstimator) == applicableLoadEstimators.end()) {
    if (not silent) {
      AutoPasLog(WARN, "Configuration is incompatible: Container {} with traversal {} does not support load estimator "
                 "{}.", container, traversal, loadEstimator);
    }
    return false;
  }

  // Check if any of the traversal's newton3 or data layout restrictions are violated.
  if (newton3 == Newton3Option::enabled) {
    const auto newton3DisabledTraversals = compatibleTraversals::allTraversalsSupportingOnlyNewton3Disabled();
    if (newton3DisabledTraversals.find(traversal) != newton3DisabledTraversals.end()) {
      if (not silent) {
        AutoPasLog(WARN, "Configuration is incompatible: Traversal {} does not support newton3 enabled.");
      }
      return false;
    }
  }
  if (newton3 == Newton3Option::disabled) {
    const auto newton3EnabledTraversals = compatibleTraversals::allTraversalsSupportingOnlyNewton3Enabled();
    if (newton3EnabledTraversals.find(traversal) != newton3EnabledTraversals.end()) {
      if (not silent) {
        AutoPasLog(WARN, "Configuration is incompatible: Traversal {} does not support newton3 enabled.");
      }
      return false;
    }
  }
  if (dataLayout == DataLayoutOption::aos) {
    const auto soaTraversals = compatibleTraversals::allTraversalsSupportingOnlySoA();
    if (soaTraversals.find(traversal) != soaTraversals.end()) {
      if (not silent) {
        AutoPasLog(WARN, "Configuration is incompatible: Traversal {} does not support AoS.");
      }
      return false;
    }
  }
  if (dataLayout == DataLayoutOption::soa) {
    const auto soaTraversals = compatibleTraversals::allTraversalsSupportingOnlyAoS();
    if (soaTraversals.find(traversal) != soaTraversals.end()) {
      if (not silent) {
        AutoPasLog(WARN, "Configuration is incompatible: Traversal {} does not support SoA.");
      }
      return false;
    }
  }

  return true;
}

template <class Functor_T, class Particle_T>
bool autopas::Configuration::hasCompatibleValues(Functor_T &functor, ParticleContainerInterface<Particle_T> &builtContainer, bool silent) const {
  // First check using the non-functor version of this function
  if (not hasCompatibleValues(silent)) {
    return false;
  }

  // Check if the required Newton 3 mode is supported by the functor
  if ((newton3 == Newton3Option::enabled and not functor.allowsNewton3()) or
      (newton3 == Newton3Option::disabled and not functor.allowsNonNewton3())) {
    if (not silent) {
      AutoPasLog(WARN, "Configuration is incompatible: The functor does not support Newton3 {}.", newton3);
      return false;
    }
  }

  // Check if the traversal can be generated using the container with all parameters set
  const auto traversalInfo = builtContainer.getTraversalSelectorInfo();

  auto isTraversalApplicable = utils::withStaticCellType<Particle_T>(
      builtContainer.getParticleCellTypeEnum(),
      [&](const auto &particleCellDummy) -> bool {
        // Can't make this unique_ptr const otherwise we can't move it later.
        auto traversalPtr =
            TraversalSelector<std::decay_t<decltype(particleCellDummy)>>::template generateTraversal<Functor>(
                traversal, functor, traversalInfo, dataLayout, newton3);

        if (traversalPtr->isApplicable()) {
          return true;
        } else {
          if (not silent) {
            AutoPasLog(WARN, "Configuration is incompatible: Generation of traversal {} failed. Check traversal's "
                       "isApplicable function for the reason.\n Configuration: [{}]", traversal, this);
          }
          return false;
        }
      });

  return isTraversalApplicable;
}

std::ostream &autopas::operator<<(std::ostream &os, const autopas::Configuration &configuration) {
  return os << configuration.toString();
}

bool autopas::Configuration::equalsDiscreteOptions(const autopas::Configuration &rhs) const {
  return container == rhs.container and traversal == rhs.traversal and loadEstimator == rhs.loadEstimator and
         dataLayout == rhs.dataLayout and newton3 == rhs.newton3 and interactionType == rhs.interactionType;
}

bool autopas::Configuration::equalsContinuousOptions(const autopas::Configuration &rhs, double epsilon) const {
  return std::abs(cellSizeFactor - rhs.cellSizeFactor) < epsilon;
}

bool autopas::operator==(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  return lhs.equalsContinuousOptions(rhs) and lhs.equalsDiscreteOptions(rhs);
}

bool autopas::operator!=(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  return not(lhs == rhs);
}

bool autopas::operator<(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  return std::tie(lhs.container, lhs.cellSizeFactor, lhs.traversal, lhs.loadEstimator, lhs.dataLayout, lhs.newton3,
                  lhs.interactionType) < std::tie(rhs.container, rhs.cellSizeFactor, rhs.traversal, rhs.loadEstimator,
                                                  rhs.dataLayout, rhs.newton3, rhs.interactionType);
}

std::istream &autopas::operator>>(std::istream &in, autopas::Configuration &configuration) {
  constexpr auto max = std::numeric_limits<std::streamsize>::max();
  in.ignore(max, ':');
  in >> configuration.interactionType;
  in.ignore(max, ':');
  in >> configuration.container;
  in.ignore(max, ':');
  in >> configuration.cellSizeFactor;
  in.ignore(max, ':');
  in >> configuration.traversal;
  in.ignore(max, ':');
  in >> configuration.loadEstimator;
  in.ignore(max, ':');
  in >> configuration.dataLayout;
  in.ignore(max, ':');
  in >> configuration.newton3;
  return in;
}

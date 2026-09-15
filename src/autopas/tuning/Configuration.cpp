/**
 * @file Configuration.h
 * @author F. Gratl
 * @date 17. Nov 2022
 */

#include "Configuration.h"

#include <tuple>
#include <type_traits>
#include <utility>

#include "autopas/containers/CompatibleCellSizeFactors.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/StringUtils.h"

namespace {

/**
 * Compare one pair of Configuration components.
 *
 * Floating point components are compared within an absolute tolerance, all others, e.g. the Options, and potential
 * future discrete components are compared exactly.
 *
 * @param lhs
 * @param rhs
 * @param epsilon Maximal allowed absolute difference, only used for floating point components.
 * @return True if the components are considered equal.
 */
bool componentEquals(const auto &lhs, const auto &rhs, double epsilon) {
  if constexpr (std::is_floating_point_v<std::remove_cvref_t<decltype(lhs)>>) {
    return autopas::utils::Math::isNearAbs(lhs, rhs, epsilon);
  } else {
    return lhs == rhs;
  }
}

}  // namespace

std::string autopas::Configuration::toString() const {
  return "{Interaction Type: " + interactionType.to_string() + " , Container: " + container.to_string() +
         " , CellSizeFactor: " + std::to_string(cellSizeFactor) + " , Traversal: " + traversal.to_string() +
         " , Load Estimator: " + loadEstimator.to_string() + " , Data Layout: " + dataLayout.to_string() +
         " , Newton 3: " + newton3.to_string() + " , VectorizationPattern: " + vecPattern.to_string() + "}";
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

bool autopas::Configuration::hasCompatibleValues() const {
  // Check if container and traversal fit together
  const auto allContainerTraversals = compatibleTraversals::allCompatibleTraversals(container, interactionType);
  if (not allContainerTraversals.contains(traversal)) {
    return false;
  }

  // Check if the selected load estimator option is applicable.
  const std::set<LoadEstimatorOption> applicableLoadEstimators =
      loadEstimators::getApplicableLoadEstimators(container, traversal, LoadEstimatorOption::getAllOptions());
  if (not applicableLoadEstimators.contains(loadEstimator)) {
    return false;
  }

  // Check if any of the traversal's newton3 or data layout restrictions are violated.
  if (newton3 == Newton3Option::enabled) {
    const auto newton3DisabledOnlyTraversals = compatibleTraversals::allTraversalsSupportingOnlyNewton3Disabled();
    if (newton3DisabledOnlyTraversals.contains(traversal)) {
      return false;
    }
  }
  if (newton3 == Newton3Option::disabled) {
    const auto newton3EnabledOnlyTraversals = compatibleTraversals::allTraversalsSupportingOnlyNewton3Enabled();
    if (newton3EnabledOnlyTraversals.contains(traversal)) {
      return false;
    }
  }
  if (dataLayout == DataLayoutOption::aos) {
    const auto soaOnlyTraversals = compatibleTraversals::allTraversalsSupportingOnlySoA();
    if (soaOnlyTraversals.contains(traversal)) {
      return false;
    }
  }
  if (dataLayout == DataLayoutOption::soa) {
    const auto aosOnlyTraversals = compatibleTraversals::allTraversalsSupportingOnlyAoS();
    if (aosOnlyTraversals.contains(traversal)) {
      return false;
    }
  }

  const auto allContainersSupportingSuper1CSF = compatibleCSFs::allContainersSupportingSuper1CSF();
  if ((not allContainersSupportingSuper1CSF.contains(container)) and cellSizeFactor > 1.0) {
    return false;
  }

  const auto allContainersSupportingSub1CSF = compatibleCSFs::allContainersSupportingSub1CSF();
  if ((not allContainersSupportingSub1CSF.contains(container)) and cellSizeFactor < 1.0) {
    return false;
  }

  // Check if the container supports the VectorizationPattern
  const auto allowedVecPatterns = compatibleVectorizationPattern::allCompatibleVectorizationPattern(container);
  if (allowedVecPatterns.find(vecPattern) == allowedVecPatterns.end()) {
    return false;
  }

  return true;
}

std::ostream &autopas::operator<<(std::ostream &os, const autopas::Configuration &configuration) {
  return os << configuration.toString();
}

bool autopas::operator==(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  constexpr double epsilon = 1e-12;

  const auto lhsTie = lhs.tie();
  const auto rhsTie = rhs.tie();

  return [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    return (componentEquals(std::get<Is>(lhsTie), std::get<Is>(rhsTie), epsilon) and ...);
  }
  (std::make_index_sequence<std::tuple_size_v<ConfigurationTie>>{});
}

bool autopas::operator!=(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  return not(lhs == rhs);
}

bool autopas::operator<(const autopas::Configuration &lhs, const autopas::Configuration &rhs) {
  return lhs.tie() < rhs.tie();
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
  in.ignore(max, ':');
  in >> configuration.vecPattern;
  return in;
}
